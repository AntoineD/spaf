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


def get_entries(entries_path):
	"""Process the list of entries."""
	
	# entries dict contains info about data entries as 2 lists
	# a list of names with key 'id' and a list of paths with key 'path'
	# Entries are sorted by increasing values
	entries = {}
	entries['id'] = [ int(entry) for entry in os.listdir(entries_path) if entry.isdigit() ]
	entries['id'].sort()

	# if there is no data entry, check we have data files right under entries directory
	# i.e. pressure, velocity.1, velocity.2, velocity.3
	if len(entries['id']) == 0:
		for variable in ['pressure', 'velocity.1', 'velocity.2', 'velocity.3']:
			variable = os.path.join(entries_path,variable)
			assert os.path.isfile(variable) is True, \
				'Data file %s doesn''t exist in %s'%(variable,entries_path)

		# entries dict is just 
		entries['id'] = [int(0)]
		entries['path'] = [entries_path]
		
		# print some info
		print 'Data entry is'+'\n\t'+entries_path+'\n'	
	
	# otherwise we process entries
	else:
		# print some info
		print '\nFound the following data entries in %s :'%entries_path
		print ' '.join([ str(entry) for entry in entries['id'] ])
	
		# entries slice
		# format is :
		#	_ none : no argument, all entries processed
		#	_ entry_id	: an integer entry index, this one only is processed
		#	_ from:to[:step] : where 'from' and 'to' are entry_id as described above,
		#		which not need to be physically existing. Entries id found between
		#		'from' and 'to' will be processed. If 'from' is omitted, all entries up
		#		to 'to' are considered. If 'to' is omitted all entries from 'from' are
		#		considered. If both are omitted, we get the same behavior as no
		#		argument. The 'step' argument is optional, a value of 1 means all entries
		#		, 2 means half the number of entries.
	
		# get from,to,step
		# we initialize with default values such as
		# all data entries will be processed
		slices = {}
		slices['from'] = entries['id'][0]
		slices['to']   = entries['id'][-1]
		slices['step'] = None
	
		# get slice arguments from command line, if any
		if len(sys.argv) == 4 :
			entries_slice = [ int(entry) for entry in sys.argv[3].split(':') ]
		
			# get slices arguments
			if len(entries_slice) == 1:
				slices['from'] = entries_slice[0]
				slices['to']   = entries_slice[0]

			if len(entries_slice) >= 2:
				if entries_slice[0] != '':
					slices['from'] = entries_slice[0]
				if entries_slice[1] != '':
					slices['to'] = entries_slice[1]

			if len(entries_slice) == 3:
				slices['step'] = entries_slice[2]
	
			# check from is a data entry
			if len(entries_slice) == 1:
				assert slices['from'] in entries['id'], \
					"%d is not a data entry"%slices['from']
			else:
				assert slices['from'] in entries['id'], \
					"Slice argument 'from'=%d is not a data entry"%slices['from']

			# check to is a data entry
			assert slices['to'] in entries['id'], \
				"Slice argument 'to'=%d is not a data entry"%(slices['to'])

			# check 'to >= from'
			assert slices['to'] >= slices['from'], \
				"In slices arguments, 'from'=%d must be >= 'to'=%d"%(slices['from'],slices['to'])

			# check step is lower than 'to - from', if present
			if len(entries_slice) == 3 and slices['to'] - slices['from'] != 0:
				assert slices['step'] <= slices['to'] - slices['from'], \
					"In slices arguments, 'step'=%d is too big"%slices['step']

		# slice the data entries
		# first find the index of the entries
		slice_from = entries['id'].index(slices['from'])
		slice_to   = entries['id'].index(slices['to']) + 1
		entries['id'] = [ entry for entry in entries['id'][slice_from:slice_to:slices['step']] ]

		# print some info
		print 'Will process the following data entries :'
		print ' '.join([ str(entry) for entry in entries['id'] ])
		
		# get absolute paths of all data entries
		entries['path'] = [ os.path.join(entries_path,str(entry)) for entry in entries['id'] ]
		
		# then check they are directories
		for entry in entries['path']:
			assert os.path.isdir(entry), \
				"Data entry %s is not a directory"%entry

	return entries




if __name__ == '__main__':
	print
	
	# set up command line parser
	usage = """%prog mesh_file data_container [from:to[:step]]
	
from, to : indexes of timesteps directories to be processed
step     : number of those timesteps to skip

If not provided, %prog will scan the results directory for timesteps directories, if none is found it will process the results directory as the one which contains data to be visualized.
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
