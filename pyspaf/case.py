#!/usr/bin/env python
# encoding: utf-8

import os,math,meshreorder,spaf,data_file
from tools import *
from opening import *
from fields import *


class Case(meshreorder.Mesh):
	"""Mesh superclass with initial, boundary conditions and openings data."""
	
	def __init__(self,file_name=None):
		# read mesh from base class
		meshreorder.Mesh.__init__(self,file_name)
		
		# print mesh info
		print self
		
		# list of openings
		self.openings = []
		
		# read openings
		self.get_openings()
		
		# dict of boundary conditions
		self.ibc = {}
		
		# fill the dict with bc zones names
		# each zone has the list of nodes of the corresponding mesh zone
		# plus variables fields
		for zone_name in self.zones:
			list_of_nodes = self.zones[zone_name].list_of_nodes()
			self.ibc[zone_name] = Fields(list_of_nodes,self.points)

		# initial condition
		# fields from the domain zone
		self.ic = self.ibc['domain']
		
		# boundary condition
		# same as ibc but without domain zone
		self.bc = self.ibc.copy()
		del self.bc['domain']
		
		# process ic
		self.get_ic()

		# process bc
		self.get_bc()
		
		# consistency
		self.get_ic_and_bc_consistent()


	def get_openings(self):
		"""Get openings data from user.
		
		If the mesh contains a zone with inlet or outlet, we ask for the opening data. For an outlet zone, we can have several physical outlets."""

		for zone_name in self.zones:
			if zone_name in ['outlet','inlet']:
				# read number
				while True:
					number = raw_input('Number of '+zone_name+' : ')
					if number in ['1','2']: break
				
				for outlet in range(int(number)):
					self.openings.append(Opening(zone_name))


	def to_file(self,data_format,dir_name):
		"""Write mesh, openings, initial conditions and boundary conditions to file.
		
		file_format,data_format : see base class doc."""
		
		# check extension
		extension = os.path.splitext(dir_name)[1]
		assert extension == '.case', \
			'Bad extension : %s : should be .case'%extension

		# create the case directory
		# check directory existence
		dir_name = check_mk_ch_dir(dir_name)
		
		print '\n','Writing case to',dir_name

		# write reordered mesh to mesh directory
		meshreorder.Mesh.to_file(self,data_format,'mesh')

		# write free nodes to mesh directory, this is for submatrices
		os.chdir('mesh')
		
		free_nodes_velocity = [ self.number_of_velocity_nodes - self.number_of_prescribed_velocity_nodes ]
		free_nodes_pressure = [ self.number_of_pressure_nodes - self.number_of_prescribed_pressure_nodes ]
		
		data_file.write('free_nodes_velocity', 'ascii', free_nodes_velocity)
		data_file.write('free_nodes_pressure', 'ascii', free_nodes_pressure)

		# for version z, write surface connectivity
		if 'z' in self.zones:
			offset = 1
			data_file.write('connectivity.z',data_format,self.zones['z'].connectivity,offset)	
		
		# write openings if any
		if len(self.openings) != 0:
			print '\twriting openings'
			file_id = open('openings','w')
			for opening in self.openings:
				opening.to_file(file_id)
			file_id.close()

		# go back to case directory
		os.chdir('..')

		# make and move to ic directory
		os.mkdir('initial_conditions')
		os.chdir('initial_conditions')

		# for each variables, write the ic
		print '\n','Writing initial conditions:'
		for variable,field in self.ic.fields.iteritems():
			data_file.write(variable,data_format,field)
		
		# go back to case directory
		os.chdir('..')
		
		# make and move to bc directory
		os.mkdir('boundary_conditions')
		os.chdir('boundary_conditions')

		# write bc
		# return if no bc
		if len(self.bc) == 0: return
		
		# we use the fact that bc could only be all same Dirichlet (all known or all unknown at pre-processing time), that nodes have been reordered according to them. So we only write the bc values, which is enough to figure out the prescribed nodes, and get the values from ic since those have been made consistent.
		print '\n','Writing boundary conditions:'
		for variable,field in self.ic.fields.iteritems():
			if variable == 'pressure':
				number_of_nodes = self.number_of_pressure_nodes
				number_of_prescribed_nodes = self.number_of_prescribed_pressure_nodes
			else:
				number_of_nodes = self.number_of_velocity_nodes
				number_of_prescribed_nodes = self.number_of_prescribed_velocity_nodes

			first_prescribed_node = number_of_nodes - number_of_prescribed_nodes
				
			# make a list with bc values
			values = field[first_prescribed_node:number_of_nodes]
			
			# check we have bc
			data_file.write(variable,data_format,values)
		
		# go back to case directory
		os.chdir('..')

		

	def get_bc(self):
		"""Create boundary conditions."""
		
		# this function is used for getting the bc from the ic fields
		def from_ic(bc_name):
			question = 'Fill '+bc_name+' boundary condition from initial condition ? [y/n] : '
			while True:
				answer = raw_input(question)
				if answer in ['y','n']: break
			return answer == 'y'
		
		for zone_name,bc_zone in self.bc.iteritems():
			if zone_name == 'wall':
				bc_zone.wall()
			
			elif zone_name == 'outlet':
				bc_zone.outlet(self.v2p_map)
			
			elif zone_name == 'shear':
				if from_ic(zone_name) == True:
					bc_zone.from_field(self.ic)
				else:
					bc_zone.linear_shear()
			
			elif zone_name == 'inlet':
				if from_ic(zone_name) == True:
					bc_zone.from_field(self.ic)
				else:					
					bc_zone.inlet(self.openings)
	

	def get_ic(self):
		"""Create initial conditions."""

		# for version z, ic and bc do not matter, we take fluid at rest
		if 'z' in self.bc and len(self.bc) == 1:
			print 'Setting fluid at rest for version z','\n'
			self.ic.rest(self.v2p_map)
			return

		string = 'Choose an initial condition :' + \
		'(take 0 if you do not know it beforehand)\n' + \
		'\t 0 : fluid at rest \n' + \
		'\t 1 : poiseuille flow in straight tube \n' + \
		'\t 2 : linear shear flow \n' + \
		'\t 3 : free surface flow on inclined plane \n' + \
		'\t 4 : constant fields \n'

		# ask initial condition
		while True:
			type_of_ic = raw_input(string)
			if type_of_ic in ['0','1','2','3','4']: break

		print ''
		
		# call ic processor
		if type_of_ic == '0':
			self.ic.rest(self.v2p_map)
	
		elif type_of_ic == '1':
			self.ic.poiseuille_straight_tube(self.openings,self.v2p_map)
	
		elif type_of_ic == '2':
			self.ic.linear_shear(self.v2p_map)
		
		elif type_of_ic == '3':
			self.ic.free_surface_on_inclined_plane(self.openings,self.v2p_map)
		
		elif type_of_ic == '4':
			self.ic.constant_fields(self.v2p_map)


	def get_ic_and_bc_consistent(self):
		"""Explain here"""

		print 'Making fields consistent'+'\n'

		# explain order	here !!!!	
		for zone in self.bc.values():
			if 'wall' in self.bc:
				zone.impose(self.bc['wall'],self.v2p_map)
			elif 'outlet' in self.bc:
				zone.impose(self.bc['outlet'],self.v2p_map)
					
		# process domain, impose all bc
		for zone in self.bc.values():
			self.ic.impose(zone,self.v2p_map,False)					
					
