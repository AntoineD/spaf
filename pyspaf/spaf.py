#!/usr/bin/env python
# encoding: utf-8

import sys,os,copy,types,array,time,socket,data_file
from tools import *

def to_spaf(mesh,data_format,dir_name):
	"""Write a mesh to file in spaf format.

	data_format : see to_file
	dir_name : name of the directory to write the mesh to"""

	# make and go to directory
	dir_name = check_mk_ch_dir(dir_name)
	
	print '\n','Writing mesh in spaf format to '+dir_name+':'
	
	# first print some infos to 'info' file : date, host and input file
	infos  = 'date : '+time.asctime()+'\n'
	infos += 'host : '+socket.gethostname()+'\n'
	if mesh.input_file is not None:
		infos += 'input file : '+mesh.input_file

	file_id = open('info','w')
	file_id.write(infos)
	file_id.close()

	# this offset is the index of the first node
	offset = 1
	
	# write points
	data_file.write('points',data_format,mesh.points)
	
	# connectivity for the computational domain
	connectivity = copy.deepcopy(mesh.zones['domain'].connectivity)
	
	# write connectivity
	data_file.write('connectivity',data_format,connectivity,offset)
		
	# open velocity to pressure map file if necessary
	if mesh.v2p_map is not None:
		data_file.write('velocity_pressure_map',data_format,mesh.v2p_map,offset)
	
	# go back to case directory
	os.chdir('..')



def from_spaf(mesh,dir_path):
	"""Read a mesh from file in spaf format.

	dir_path : path to directory containing the mesh
	
	Read the geometric table, the connectivity and the velocity to pressure map.	
	Note : spaf stores with 1-offset, 0-offset conversion is done"""

	print 'Reading spaf mesh from '+dir_path+'\n\t',

	# get the current directory
	caller_directory = os.getcwd()
	
	try:
		os.chdir(dir_path)
	except:
		exit('cannot go to'+dir_path)

	# read points
	print 'reading points,',
	points = data_file.read('points')

	# read connectivity
	print 'connectivity,',
	connectivity = data_file.read('connectivity')

	# read velocity pressure map
	print 'velocity to pressure map\n'
	v2p_map = data_file.read('velocity_pressure_map')
	
	# fill mesh with points
	for i in range(0,len(points),3):
		mesh.add_point(points[i:i+3])

	# add domain zone and fill connectivity
	# make it 0-offset
	mesh.add_zone('domain')
	for i in range(0,len(connectivity),10):
		cell = [ node-1 for node in connectivity[i:i+10] ]
		mesh.add_cell('domain',cell)
	
	# make v2p_map 0-offset
	mesh.v2p_map = [ node-1 for node in v2p_map ]

	# go back to caller directory
	os.chdir(caller_directory)
	