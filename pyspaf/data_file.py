#!/usr/bin/env python
# encoding: utf-8

import sys, types, array
from tools import *

# This file has tools to write and read to and from spaf file formats

def write(file_name,format,data,offset=0,integer_type='int'):
	"""Writes data to a file.

	file_name : a file name
	format : 'ascii' or 'binary'
	data : list of data or list of list of data
	offset : value to offset the data entries
	integer_type : integer data type that will be written
	"""

	# check data format
	assert format in ['binary','ascii'], \
		'Bad data format : %s'%format

	# check data is not empty
	assert len(data) != 0, \
		'Empty data'

	# check integer_type
	assert integer_type in ['int','uint'], \
		'Bad integer type : %s'%integer_type

	# store data in the list data_array
	# if it was a list of list, convert it to a successive list, we check type of first entry to get the entry type, this means we assume that all entries of the list data have the same type, we want a list with no list as entries
	if isinstance(data[0],types.ListType):
		data_array = []
		for datum in data:
			data_array.extend(datum)
	else:
		data_array = data
	
	data_type = None

	# get data type (python FloatType = c double)
	# we assume that all entries have the same type as the first one
	if isinstance(data_array[0],types.FloatType):
		data_type = 'd'
		ctype = 'double'
	elif isinstance(data_array[0],types.IntType):
		if integer_type == 'int':
			data_type = 'i'
			ctype = 'int'
		else:
			data_type = 'I'
			ctype = 'uint'

	assert data_type is not None, \
		'Data must be double or int'

	# open file
	file_id = open(file_name,'w')
	print '\twriting',file_name

	# write header
	write_name_value(file_id,'number',len(data_array))
	write_name_value(file_id,'format',format)
	write_name_value(file_id,'type',ctype)

	# offset
	if offset != 0:
		data_array_offset = [ x+offset for x in data_array ]
	else:
		data_array_offset = data_array

	# write data
	if format == 'ascii':
		# convert entries to strings and write to file
		file_id.writelines([ repr(x)+'\n' for x in data_array_offset ])		
	else:
		# append to header
		write_name_value(file_id,'byte_order',sys.byteorder)		
	
		# create array with binary data and write it
		temp_array = array.array(data_type)
		temp_array.fromlist(data_array_offset)
		temp_array.tofile(file_id)



def read(file_name,verbose=False):
	"""Read a data file in spaf format.

	verbose : will print info if true
	
	return data list or None if a problem occurs"""

	try:
		file_id = open(file_name,'r')
	except:
		if verbose: print 'Can not open',file_name
		return None

	# get the number of points
	number = int(parse_file_line_kv(file_id,'number'))

	# get data format
	format = parse_file_line_kv(file_id,'format')

	assert format in ['ascii','binary'], \
		"Bad format : %s"%format

	# get data type
	ctype = parse_file_line_kv(file_id,'type')

	assert ctype in ['float','double','int','uint'], \
		"Bad ctype : %s"%ctype

	# print some info if required
	if verbose:
		print number,'entries,',format,'format,',ctype,'data type'

	# create the array that will hold the data, according to the ctype
	if ctype == 'float':
		data = array.array('f')
	elif ctype == 'double':
		data = array.array('d')
	elif ctype == 'int':
		data = array.array('i')
	elif ctype == 'uint':
		data = array.array('I')

	# for binary format we need to read the byte order and data type
	if format == 'binary':
		# get byte order
		byte_order = parse_file_line_kv(file_id,'byte_order')

		assert byte_order in ['little','big'], \
			"Bad format : %s"%byte_order

		# only little endian implemented so far
		assert byte_order == 'little', \
			"Big endian not implemented"

		data.fromfile(file_id,number)
	else:
		# create and fill a list of entries
		data_list = []

		while len(data_list) < number:
			data_list.append(float(file_id.readline()))
		
		# fill data from list
		data.fromlist(data_list)

	return data
