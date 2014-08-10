#!/usr/bin/env python
# encoding: utf-8
"""
vtk.py

Created by Antoine Dechaume on 2008-02-17.
Copyright (c) 2008 __MyCompanyName__. All rights reserved.
"""

import os,data_file,vtk

def to_vtk_file(mesh,data_format,file_name):
	"""Write a mesh to vtk xml format.
	
	data_format : see to_file
	file_name : name of the file to write the mesh to
	
	Note that vtk has the same local node ordering as spaf,
	hence no need to use a mapping."""

	# fill in vtk mesh
	writer = mesh.to_vtk(data_format,file_name)
	
	# write to file
	writer.Write()
	

def to_vtk(mesh,data_format,file_name):
	"""Fill in vtk mesh.
	
	data_format : see to_file
	file_name : name of the file to write the mesh to
	
	Note that vtk has the same local node ordering as spaf,
	hence no need to use a mapping."""

	# check data format
	assert data_format in ['binary','ascii'], \
		'Bad data format : %s'%data_format

	# check file extension
	extension = os.path.splitext(file_name)[1]
	assert extension == '.vtu', \
		'Bad extension : %s : should be .vtu'%extension

	print '\t','writing to',file_name

	# map from cell name to identification number
	cell_name_to_id = { \
		'tetrahedron_04' : 10, \
		'triangle_03'    : 5,  \
		'tetrahedron_10' : 24, \
		'triangle_06'    : 22 }

	# fill points
	print '\t','filling points'
	points = vtk.vtkPoints()
	for point in mesh.points:
		points.InsertNextPoint(point[0],point[1],point[2])

	# fill cells
	print '\t','filling cells'
	cells = vtk.vtkCellArray()
	for	cell in mesh.zones['domain'].connectivity:		
		cells.InsertNextCell(len(cell))
		for node in cell:
			cells.InsertCellPoint(node)

	# create grid
	grid = vtk.vtkUnstructuredGrid()
	grid.SetPoints(points)
	
	# get type of cell, make sure it takes 2 char with 0 prefix if necessary
	cell_type_name = mesh.zones['domain'].cell_type
	# cell_type_name = 'tetrahedron_'+str(len(mesh.zones['domain'].connectivity[0])).zfill(2)
	cell_type = cell_name_to_id[cell_type_name]
	grid.SetCells(cell_type,cells)

	# write data to file
	print '\t','writing file'
	writer = vtk.vtkXMLUnstructuredGridWriter()

	if data_format == 'binary':
	# turn base64 off, from vtk doc :
	# Get/Set whether the appended data section is base64 encoded.  If
	# encoded, reading and writing will be slower, but the file will be
	# fully valid XML and text-only.  If not encoded, the XML
	# specification will be violated, but reading and writing will be
	# fast.  The default is to do the encoding.
		writer.EncodeAppendedDataOff()
		# it seems we get no benefit from specifing compression level
		# as compared to default setting
		# compressor = vtk.vtkZLibDataCompressor()
		# compressor.SetCompressionLevel(9)
		# writer.SetCompressor(compressor)
	else:
		writer.SetDataModeToAscii()

	writer.SetFileName(file_name)
	writer.SetInput(grid)

	return writer
	
	
def write_timesteps(data_format,mesh,variable,entries,high_precision=False):
	"""Create vtk file from data entries.

	data_format : ascii or binary
	mesh : mesh of a variable (pressure or velocity)
	files : name of the data files to process
	entries : list of data entries
	high_precision : select whether to write float or double real numbers

	We go through all data entries, then read data files and append them to the vtk file."""

	if variable == 'pressure':
		files = ['pressure']
	elif variable == 'velocity':
		files = ['velocity.1','velocity.2','velocity.3']
	else:
		abort('Bad variable')

	print '\nProcessing',variable

	# get the number of points
	nb_of_points = len(mesh.points)

	# set file name of the file
	file_name = variable+'.vtu'

	# write mesh to vtk file, in binary
	writer = to_vtk(mesh,data_format,file_name)

	# now append time steps to vtk file
	nb_of_timesteps = len(entries['id'])
	writer.SetNumberOfTimeSteps(nb_of_timesteps)

	# start writing time steps
	writer.Start()

	# data array
	data_handle = writer.GetInput().GetPointData()

	# loop over filtered entries
	print '\nProcessing entries :'+'\n\t',

	for index,data_directory in enumerate(entries['path']):
		print entries['id'][index],

		# holds data for each files
		data = {}

		# read time
		time = 0.
		
		# path to time file
		path_to_time = os.path.join(data_directory,'time')
		if os.path.isfile(path_to_time):
			time = data_file.read(path_to_time)[0]

		# read data
		for datafile in files:
			# get path and check file is there
			data_path = os.path.join(data_directory,datafile)
			assert os.path.isfile(data_path), \
				"File %s is missing"%data_path

			# read all data for this iteration
			data[datafile] = data_file.read(data_path)

			# filter velocity if necessary
			if variable == 'velocity' and hasattr(mesh,'p2v_map'):
				filtered_data = []
				for old in mesh.p2v_map:
					filtered_data.append(data[datafile][old])
				data[datafile] = filtered_data

			# check the number of points
			# if not hasattr(mesh,'p2v_map'):
			# 	assert len(data[datafile]) == nb_of_points, \
			# 		"Data set %s has a wrong number of points : %d != %d" \
			# 		%(datafile,len(data[datafile]),nb_of_points)

		# write vtk
		if high_precision:
			field = vtk.vtkDoubleArray()
		else:
			field = vtk.vtkFloatArray()

		field.SetName(variable)
		field.SetNumberOfComponents(len(files))
		field.SetNumberOfTuples(nb_of_points)

		if variable == 'velocity':
			for i in range(nb_of_points):
				field.SetTuple3(i,data['velocity.1'][i], data['velocity.2'][i], data['velocity.3'][i])

			data_handle.SetVectors(field)

		elif variable == 'pressure':
			for i in range(nb_of_points):
				field.SetTuple1(i,data['pressure'][i])

			data_handle.SetScalars(field)

		writer.WriteNextTime(time)

	print

	writer.Stop()			
