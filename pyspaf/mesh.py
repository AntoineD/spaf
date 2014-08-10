#!/usr/bin/env python
# encoding: utf-8

import os, spaf, gmsh
from zone import *
from zone_definitions import *


class Mesh:
	"""Mesh class.
	
	A mesh is defined by following items:
		geometric table: a list of points, each point is defined by the list of 3 coordinates.
		zones: at least the global computational domain, and eventually the boundary conditions. The key of each zone of the dictionary is the name of the zone. A zone name must be valid, see add_zone.
		map from velocity node indexes to pressure nodes indexes
	"""
	
	def __init__(self,input_file=None):
		"""input_file : name of the file that contain the mesh, optional, for reading a mesh from file."""
		
		# store input file for informations
		self.input_file = None
		
		if input_file is not None:
			self.input_file = os.path.abspath(input_file)
		
		# geometric table, each element of this list is an list of 3 coordinates.
		# points = [ ... , [x,y,z] , ... ]
		self.points = []
		
		self.number_of_velocity_nodes = None

		# dict of zones.
		# a key is the name of the zone, which must be valid, see add_zone
		# zones = { ... , 'name of the zone' : instance of zone class , ... }
		self.zones = {}
		
		# mapping from velocity node numbering to pressure one
		# usefull because velocity usually is 2nd order and pressure is 1st
		# also because of renumbering according to pressure bc
		self.v2p_map = None
		
		# read mesh from file if required
		if input_file is not None:
			self.from_file(input_file)
			self.number_of_velocity_nodes = len(self.points)
			
			# create the velocity to pressure map if not read from file
			if self.v2p_map is None:
				self.get_v2p_map()
		
			# count number of pressure nodes
			self.number_of_pressure_nodes = len(self.v2p_map) - self.v2p_map.count(-1)
			
		

	def __str__(self):
		"""Print mesh informations."""

		string = 'Mesh info :\n'
		string += '\t'+str(self.number_of_velocity_nodes)+' nodes \n'

		if self.v2p_map is not None:
			string += '\t'+str(self.number_of_pressure_nodes)
			string += ' pressure nodes\n'

		string += '\n'

		# print entities elements info
		for name,zone in self.zones.items():
			string += '\t'+name+' : '+zone.__str__()

		return string

		
	def add_zone(self,zone_name,cell_type=None):
		"""Add a zone to the mesh.
		
		zone_name : a string which is the name of the zone
		cell_type : type of the cells in the zone, optional.
		
		The name must be valid, as defined in zone_definitions.
		This last argument is for unusual mesh such as the pressure one for spaf2viz."""
		
		# check validity of the name of the zone
		assert zone_name_is_valid(zone_name)

		# set cell_type
		# if no type has been given, we get it from zones definitions
		if cell_type is None:
			cell_type = get_zone_cell_type(zone_name)

		# add a zone
		self.zones[zone_name] = Zone(cell_type)


	def add_point(self,coordinates):
		"""Add a point to the geometric table.
		
		coordinates : a list of the 3 coordinates of a point"""	
			
		# check number of coordinates
		assert len(coordinates) == 3, \
			'Bad number of coordinates : %d'%len(coordinates)
		
		# add the list of a point coordinates to the list of points
		# converted to floats
		self.points.append([ float(x) for x in coordinates ])
	
	
	def add_cell(self,zone_name,nodes):
		"""Add a cell to a zone.
		
		zone_name : a string which is the name of the zone
		nodes : a list of nodes defining the cell, must be in spaf ordering"""
		
		# check zone name exists in zones dict
		assert zone_name in self.zones, \
			'Element has invalid zone name : %d'%zone_name

		# add cell to zone
		self.zones[zone_name].append(nodes)
		
		
	def get_v2p_map(self):
		"""Create the map from velocity nodes indexes to pressure nodes indexes.
		
		A velocity node which is not a pressure is associated to -1 through the
		mapping, otherwise it is associated to a unique pressure node index."""
		
		# initialize to -1
		self.v2p_map = [-1] * self.number_of_velocity_nodes

		# We loop over all cells of the domain zone connectivity, then for each cell
		# we set velocity nodes that are vertex of the cell to be pressure nodes
		# if they haven't been already processed
		pressure_index = 0
		for cell in self.zones['domain'].connectivity:
			for node in cell[0:4]:
				if self.v2p_map[node] == -1:
					self.v2p_map[node] = pressure_index
					pressure_index += 1


	def from_file(self,input_file):
		"""Dispatches to file reader according to file extension.
		
		input_file : name of the input file"""
		 
		# get file extension or base name
		extension = os.path.splitext(input_file)[1]
		base_name = os.path.basename(input_file)

		if extension == '.msh':
			reader = gmsh.from_gmsh
		elif extension == '.mesh' or base_name == 'mesh':
			reader = spaf.from_spaf
		else:
			exit('Bad file name : %s : should be .msh or .case or just mesh'%input_file)
		
		reader(self,input_file)


	def to_file(self,data_format,output_file):
		"""Dispatch to file writer according to file format.

		output_file : format of the file, spaf or vtk
		data_format : format of the data, ascii or binary"""

		# dispatch to reader according to extension
		extension = os.path.splitext(output_file)[1]

		# dispatch to writer
		if extension == '.vtu':
			import Vtk
			writer = Vtk.to_vtk_file
		elif extension == '.msh':
			import gmsh
			writer = gmsh.to_gmsh
		else:
			writer = spaf.to_spaf

		# call selected writer and return path to file
		return writer(self,data_format,output_file)
	

	def to_vtk(self,data_format,file_name):
		"""Write in vtk format, leave file open and returns vtk file handle.

		data_format : format of the data, ascii or binary
		file_name : name of the output file
		returns a writer handle for future usage."""

		import Vtk
		return Vtk.to_vtk(self,data_format,file_name)
