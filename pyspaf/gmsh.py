#!/usr/bin/env python
# encoding: utf-8

"""This file is used for reading a mesh in gmsh format."""

import array, os
from tools import *


def from_gmsh(mesh,file_name):
	"""Read a mesh from a gmsh file.
	
	file_name : a file name"""

	# this dict translates gmsh cell identification to the one used here
	cell_id_to_name = { \
		11 : 'tetrahedron_10', \
		9  : 'triangle_06' }

	# this dict gives, for each cell type, the map from gmsh node ordering to spaf one
	nodes_ordering = {
		'tetrahedron_10' : [0,1,2,3,4,5,6,7,9,8], \
		'triangle_06'    : [0,1,2,3,4,5] }
		
	# in the following sketch, node 3 is in the background
	#             1     
	#           /|\		
	#         /  9 \	
	#       4   |   5	
	#     /     3    \	
	#   /   7`    `8  \	
	#  0-------6-------2


	try:
		file_id = open(file_name,'r')
	except:
		print 'cannot open', file_name
		exit()
	
	print '\nReading gmsh mesh file'

	# check mesh format section
	assert parse_file_line(file_id)[0] == '$MeshFormat', \
		'Missing $MeshFormat section'

	# read format info, this reader reads format 2
	line = parse_file_line(file_id)
	assert line[0] == '2', 'File is not version 2'

	print '\t',

	# read encoding format
	if line[1] == '1':
		print 'binary encoding'
		format = 'binary'
		#  check size of real numbers, should be double
		assert line[2] == '8', 'Bad size of double : %d'%line[2]
		# check endianness
		one = array.array('i')
		one.fromfile(file_id, 1)
		assert one[0] == int(1), 'Bad endianness'
		# finish reading line
		file_id.readline()
	else:
		print 'ASCII encoding'
		format = 'ascii'
	
	# check mesh format section end
	assert parse_file_line(file_id)[0] == '$EndMeshFormat', \
		'Missing $EndMeshFormat section'

	# check physical entity section
	assert parse_file_line(file_id)[0] == '$PhysicalNames', \
		'Missing $PhysicalNames section'

	# read number of zones
	number_of_zones = int(parse_file_line(file_id)[0])
	print '\t', number_of_zones, 'zones'

	# map from zone id to zone name
	zone_id_to_name = {}

	# read zones info (remove double quotes) and
	# create zone instances with default cell type
	for index in range(number_of_zones):
		line = parse_file_line(file_id)
		zone_name = line[1].replace('"','')
		zone_id_to_name[int(line[0])] = zone_name
		mesh.add_zone(zone_name)

	# check physical entity section end
	assert parse_file_line(file_id)[0] == '$EndPhysicalNames', \
		'Missing $EndPhysicalNames section'

	# check nodes section
	assert parse_file_line(file_id)[0] == '$Nodes', \
		'Missing $Nodes section'

	# read number of nodes
	number_of_nodes = int(parse_file_line(file_id)[0])
	print '\t', number_of_nodes, 'nodes'

	# read geometric table

	# loop over the number of nodes, a line format is
	# (node index) (x coordinate) (y coordinate) (z coordinate)
	# we check that the node index matches the loop index
	for index in range(1,number_of_nodes+1):
		if format == 'binary':
			idx = array.array('i')
			idx.fromfile(file_id, 1)
			assert index == idx[0], \
				'Bad node index : %d, %d'%(index,idx[0])
			line = array.array('d')
			line.fromfile(file_id, 3)
		else:
			line = parse_file_line(file_id)
			assert index == int(line[0]), \
				'Bad node index : %d, %d'%(index,int(line[0]))
			line.pop(0)

		mesh.add_point(line)

	# finish reading line
	if format == 'binary': file_id.readline()

	# check nodes section end
	assert parse_file_line(file_id)[0] == '$EndNodes', \
		'Missing $EndNodes section'

	# check element section
	assert parse_file_line(file_id)[0] == '$Elements', \
		'Missing $Elements section'

	# read total number of elements
	total_number_of_elements = int(parse_file_line(file_id)[0])
	print '\t', total_number_of_elements, 'elements'

	if format == 'binary':
		elements_count = 0
		while elements_count < total_number_of_elements:
			# read header
			header = array.array('i')
			header.fromfile(file_id, 3)
			element_type_id = header[0]
			number_of_elements = header[1]
			tags_nb = header[2]
			
			# check type of cell is known
			assert element_type_id in cell_id_to_name, \
				'Bad cell type : %d'%element_type_id

			# get number of nodes per elements
			number_of_nodes = int(cell_id_to_name[element_type_id][-2:])

			# read data
			for element in range(number_of_elements):
				elements_count += 1
				line = array.array('i')
				line.fromfile(file_id,1+tags_nb+number_of_nodes)
				zone_name = zone_id_to_name[line[1]]
				nodes = reorder_nodes(line[1+tags_nb:],nodes_ordering)
				# offset node indexes by 0
				nodes = [ node-1 for node in nodes ]
				mesh.add_cell(zone_name,nodes)
	
		# check the correct number of elements has been read
		assert total_number_of_elements == elements_count, \
			'Element count mis-match : %d != %d'%(total_number_of_elements,elements_count)
	else:
		for index in range(total_number_of_elements):
			line = [ int(x) for x in parse_file_line(file_id) ]
			tags_nb = line[2]
			zone_name = zone_id_to_name[line[3]]
			nodes = reorder_nodes(line[3+tags_nb:],nodes_ordering)
			# offset node indexes by 0
			nodes = [ node-1 for node in nodes ]
			mesh.add_cell(zone_name,nodes)

	# finish reading line
	if format == 'binary': file_id.readline()

	# check element section end
	assert parse_file_line(file_id)[0] == '$EndElements', \
		'Missing $EndElements section'

	# close file
	file_id.close()
	
	print ''



def to_gmsh(mesh,format,file_name):
	"""Read a mesh from a gmsh file.

	file_name : a file name"""

	assert format in ['ascii','binary'], \
		'Bad format: %s'%format

	file_name = check_file_exist(file_name)
	
	try:
		file_id = open(file_name,'w')
	except:
		print 'cannot open', file_name
		exit()
	
	# $MeshFormat section
	file_id.write('$MeshFormat\n')
	
	# print version number
	string = ' 2'
	
	# print file type
	if format == 'ascii':
		string += ' 0'
	else:
		string += ' 1'
		
	# print data size
	string += ' 8\n'
	file_id.write(string)
	
	# print one binary
	if format == 'binary':
		# create array with binary data and write it
		temp_array = array.array('i')
		temp_array[0] = 1
		temp_array.tofile(file_id)
		
	# $EndMeshFormat section
	file_id.write('$EndMeshFormat\n')
	
	# $Nodes section
	file_id.write('$Nodes\n')
	
	# write number of nodes
	file_id.write(str(mesh.number_of_velocity_nodes))
	
	if format == 'ascii':
		for node, point in enumerate(mesh.points):
			string = str(node)+' '+str(point[0])+' '+str(point[1])+' '+str(point[2])+'\n'
			file_id.write(string)
	else:
		pass
	
	# $EndNodes section
	file_id.write('$EndNodes\n')

	# $Elements section
	file_id.write('$Elements\n')
	
	# from now on we will use an alias to the domain zone
	zone = mesh.zones['domain']
	
	# write number of elements
	file_id.write(str(zone.number_of_elements()))
	
	nodes_ordering = [0,1,2,3,4,5,6,7,9,8]
	
	if format == 'ascii':
		for element, connectivity in enumerate(zone.elements):
			# element number, type (only tetrahedron) and tag number
			string = str(element)+' 11 0'
			
			# append connectivity with gmsh order and offset
			for node in connectivity:
				string += ' '+str(nodes_ordering[node]+1)
				
			file_id.write(string+'\n')
	else:
		pass
	
	# $EndElements section
	file_id.write('$EndElements\n')

	# $NodeData section
	file_id.write('$NodeData\n')
	
	# 1 string tag ('velocity'), 0 real tag, 3 integer tag ( timestep 0, 3 components, as many nodes in the mesh)
	file_id.write('1\n'+'velocity\n'+'0\n'+'3\n'+'0\n'+'3\n'+str(mesh.number_of_nodes)+'\n')
	
	if format == 'ascii':
		for node, connectivity in enumerate(mesh.points):
			# element number, type (only tetrahedron) and tag number
			string = str(element)+' 11 0'
			
			# append connectivity with gmsh order and offset
			for node in connectivity:
				string += ' '+str(nodes_ordering[node]+1)
				
			file_id.write(string+'\n')
	else:
		pass

	# $EndNodeData section
	file_id.write('$EndNodeData\n')

	
	
	
	
	