#!/usr/bin/env python
# encoding: utf-8
"""
spaf2viz.py

Converts spaf data to vtk for visualization.
"""

__author__ = "Antoine Dechaume"
__author_email__ = "dechaume@ualberta.ca"
__version__ = "r22"

import sys,os,os.path,data_file
from mesh import *
from tools import *
from optparse import OptionParser


def get_mesh(path_to_mesh):
	"""docstring for get_mesh"""
	
	# check files are present
	for data_file in ['connectivity','points','elements']:
		path_to_file = os.path.join(path_to_mesh,data_file)
		if os.path.isfile(path_to_file) == false:
			exit('No'+path_to_file+'file')
			
	# read connectivity
	path_to_file = os.path.join(path_to_mesh,'connectivity')
	connectivity = data_file.read(path_to_file)

	assert connectivity not None,
		'Cannot read connectivity file : %s'%path_to_file

	# read points
	path_to_file = os.path.join(path_to_mesh,'points')
	points = data_file.read(path_to_file)

	assert points not None,
		'Cannot read points file : %s'%path_to_file

	# read elements
	path_to_file = os.path.join(path_to_mesh,'elements')
	elements = data_file.read(path_to_file)

	assert elements not None,
		'Cannot read elements file : %s'%path_to_file
	





if __name__ == '__main__':
	print
	
	# set up command line parser
	usage = """%prog mesh_file data_container [from:to[:step]]
	
from, to : indexes of timesteps directories to be processed
step     : number of those timesteps to skip

If not provided, spaf2viz will scan the results directory for timesteps directories, if none is found it will process the results directory as the one which contains data to be visualized.
"""
	version = "%prog "+__version__
	parser = OptionParser(usage=usage, version=version)

	# output option
	# it is the name of the output directory
	parser.add_option(
		"-o","--output",
		type="string",
		action="store",
		dest="output",
		default="viz_files",
		help="write to DIR, default is viz_files",
		metavar="DIR")

	# data format option
	# set to ascii when invoqued
	parser.add_option(
		"-a","--ascii",
		action="store_const",
		const="ascii",
		dest="data_format",
		default="binary",
		help="write in ascii format, default is binary" )
		
	# debug option
	# write real values to double when invoqued
	parser.add_option(
		"-g","--debug",
		action="store_true",
		dest="debug",
		default=False,
		help="write real values to double, default is float" )

	# parse command line
	(options,args) = parser.parse_args()

	# process arguments
	# check we have the input mesh file
	if len(args) not in [2,3]:
		parser.error('Bad number of arguments')
    
	# get absolute path to mesh and entries
	mesh_path = os.path.abspath(args[0])
	entries_path = os.path.abspath(args[1])
	
	
	# process list of entries
	entries = get_entries(entries_path)

	
	# Now we create meshes, one for the velocity and one for the pressure,
	# as a dict with 2 keys : 'pressure' and 'velocity'
	mesh = {}
	
	# read velocity mesh
	mesh['velocity'] = Mesh(mesh_path)
	
	# create pressure mesh
	# it has first order elements only
	mesh['pressure'] = Mesh()

	# fill pressure mesh by adding a domain zone
	mesh['pressure'].add_zone('domain','tetrahedron_04')
	
	# just an alias to velocity to pressure map
	v2p_map = mesh['velocity'].v2p_map
	
	# add connectivity to this zone, from the velocity one
	# as pressure nodes are the first 4 velocity nodes
	for cell in mesh['velocity'].zones['domain'].connectivity:
		new_cell = [ v2p_map[node] for node in cell[:4] ]
		mesh['pressure'].add_cell('domain',new_cell)

	# get pressure to velocity map by inverting the velocity to pressure map
	# and filtering -1 values
	p2v_map = []
	for old,new in enumerate(v2p_map):
		if new != -1:
			p2v_map.insert(new,old)
	
	# fill geometric table of the pressure mesh
	# by filtering the velocity one through the pressure to velocity map
	for old in p2v_map:
		mesh['pressure'].add_point(mesh['velocity'].points[old])


	# ask order of velocity elements to write
	# by choosing 1 we save a lot a space
	# and share the mesh with the pressure
	while True:
		answer = raw_input('Store velocity at pressure nodes ? [y/n] ')
		if answer in ['y','n']: break

	# if required, set velocity mesh as an alias to the pressure mesh
	# we'll need the pressure to velocity map too, which we add to mesh class
	if answer == 'y':
		mesh['velocity'] = mesh['pressure']
		mesh['velocity'].p2v_map = p2v_map

	# move to output directory
	options.output = check_mk_ch_dir(options.output)

	# import Vtk now as it takes some time to load
	import Vtk	

	# process timesteps data
	for variable in mesh:
		print mesh[variable]
		Vtk.write_timesteps(options.data_format,mesh[variable],variable,entries,options.debug)
