#!/usr/bin/env python
# encoding: utf-8
"""
tools.py

Created by Antoine Dechaume on 2008-02-17.
Copyright (c) 2008 __MyCompanyName__. All rights reserved.
"""

import math,os

#==============================================================================
# file system stuff
#==============================================================================
def check_mk_ch_dir(dir_name,extension=None):
	"""Check a directory exists, if not make and go to it, otherwise ask new name.
	Returns the eventual new directory name."""
	
	# check directory existence
	while os.path.isdir(dir_name):
		print '\n'
		question = dir_name+' already exists, overwrite it ? [y/n] : '
		while True:
			answer = raw_input(question)
			if answer in ['y','n']: break
		
		if answer == 'y':
			# delete directory
			import shutil
			shutil.rmtree(dir_name) 
		else:
			# get new directory name
			question = 'New file name : '
			while True:
				answer = raw_input(question)
				if answer != dir_name: break
			dir_name = answer

	# make and go to mesh directory
	os.mkdir(dir_name)
	os.chdir(dir_name)
	
	return dir_name

def check_file_exist(file_name,extension=None):
	"""Check a file exists, if not return, otherwise ask new name.
	Returns the eventual new file name."""

	# check directory existence
	while os.path.isfile(file_name):
		print '\n'
		question = file_name+' already exists, overwrite it ? [y/n] : '
		while True:
			answer = raw_input(question)
			if answer in ['y','n']: break

		if answer == 'y':
			# delete directory
			os.remove(file_name) 
		else:
			# get new directory name
			question = 'New file name : '
			while True:
				answer = raw_input(question)
				if answer != file_name: break
			file_name = answer

	return file_name

#==============================================================================
# linear algebra stuff
#==============================================================================

def get_vector(name):
	"""Read a vector from user.
	
	name : a string to print
	return : a list of 3 float"""
	
	# read the components, make sure they are 3 real numbers
	while True:
		vector = raw_input(name+' coordinates : ')
		vector = vector.split(' ')
		if len(vector) == 3:
			test = False
			for coord in vector:
				if not string_is_float(coord):
					test = True
			if test == False: break

	return [ float(coord) for coord in vector ]

def scalar_product(vector1,vector2):
	"""Compute the scalar product of 2 vectors.
	
	vector1,vector2 : 2 vectors with 3 components
	returns the scalar product."""
	
	value = 0.
	for i in range(3):
		value += vector1[i] * vector2[i]
	return value

def norm(vector):
	"""Compute the norm of a vectors.
	
	vector : a vector with 3 components
	returns the norm."""
	
	return math.sqrt(scalar_product(vector,vector))


#==============================================================================
# string to numeric type testing
#==============================================================================

def string_is_int(string):
	"Check whether a string is an integer or not, return bool"
	try:
		x = int(string)
		return True
	except ValueError:
		return False

def string_is_float(string):
	"Check whether a string is a real value or not, return bool"
	try:
		x = float(string)
		return True
	except ValueError:
		return False

#==============================================================================
# Parse tools
#==============================================================================
def parse_line(line):
	"""Parse a line.
	
	line : a string
	return : the list of words in line"""
	
	return line.split()

def parse_file_line(file_id):
	"""Parse the current line of a file.
	
	file_id : a file object
	return : see parse_line return"""
	
	return parse_line(file_id.readline())

def parse_file_line_kv(file_id,key):
	"""Parse the current line of a file in key = value format.

	file_id : a file object
	key : key field string
	return : value"""

	line = parse_file_line(file_id)
	
	assert line[0] == key, \
		'Bad key, should be : %s'%key
	
	assert line[1] == '=', \
		'Bad key value separator, should be ='

	return line[2]

#==============================================================================
# write tools
#==============================================================================
def write_name_value(file_id,name,value):
	"""Write 'name = value' to a file.

	file_id : a file object
	name : a string
	value : anything"""

	file_id.write(name+' = '+str(value)+'\n')

#==============================================================================
# Nodes reordering
#==============================================================================

def reorder_nodes(nodes,local_ordering):
	"""Reorder the node indexes of a cell to spaf ordering.

	nodes : a list of global indexes of nodes that defines a cell
	local_ordering : a dict of cell types to local nodes mapping, see map definition
	return : reordered list of nodes"""

	# select cell type by finding the number of nodes in nodes list
	cell_ordering = None

	for cell in local_ordering:
		if int(cell[-2:]) == len(nodes):
			cell_ordering = local_ordering[cell]
			break

	assert cell_ordering != None, \
		'Number of nodes does not correspond to element in ordering map'

	# re-order nodes to spaf ordering
	return [ nodes[i] for i in cell_ordering ]
