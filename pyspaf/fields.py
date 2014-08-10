#!/usr/bin/env python
# encoding: utf-8
"""
fields.py

Created by Antoine Dechaume on 2008-04-28.
Copyright (c) 2008 __MyCompanyName__. All rights reserved.

This file contains code for managing velocity and pressure values in a set of nodes. It is used for prescribing initial and boundary conditions.
"""

import copy
from tools import *



def select_opening(openings,kind):
	"""Search for a specified opening kind and return it."""

	opening = None

	# search for outlet
	for op in openings:
		if op.kind == kind:
			opening = op
			break

	assert opening != None, \
		'no '+kind+' in openings'

	return opening


def get_non_zero_float(question):
	"""Ask a question, get a float, check it is non zero and return it.
	
	question : string to display
	"""

	while True:
		number = raw_input(question)
		if string_is_float(number) and float(number) != 0: break

	# convert to float
	return float(number)


def get_float(question):
	"""Ask a question, get a float and return it.
	
	question : string to display
	"""

	while True:
		number = raw_input(question)
		if string_is_float(number): break

	# convert to float
	return float(number)


def get_direction(string):
	"""Ask direction and return it."""

	while True:
		direction = raw_input('direction of '+string+' : ')
		if direction in ['1','2','3']: break

	return direction
	
	
	
	

class Fields:
	"""Class for field data.
	
	A Field contains velocity and pressure values at certain nodes. It is defined by a list of nodes, with corresponding coordinates in the point list, and a dict of variables with pressure and 3 velocities components."""

	def __init__(self,list_of_nodes,points):
		"""
		list_of_nodes : nodes that define the field
		points : map from nodes indexes to nodes coordinates
		"""
		
		# list of the variables of a field
		field_names = ['pressure','velocity.1','velocity.2','velocity.3']

		# initialize the dict of variables
		self.fields = {}
		
		# each element of the dict is a list of values
		# each element of this list is the value of the variable
		# at the node defined in list_of_nodes
		for field in field_names:
			self.fields[field] = []
		
		# make an alias points and list of nodes
		self.points = points
		self.nodes = list_of_nodes
		

	def wall(self):
		"""Wall condition.
		
		Velocity : Dirichlet homogeneous.
		Pressure : None.
		"""

		for variable in ['velocity.1','velocity.2','velocity.3']:
			self.fields[variable] = [0.] * len(self.nodes)


	def outlet(self,v2p_map):
		"""Outlet condition.
		
		v2p_map : velocity to pressure map
		
		Velocity : Neumann homogeneous, but None actually.
		Pressure : Dirichlet homogeneous.
		"""

		for node in self.nodes:
			if v2p_map[node] != -1:
				self.fields['pressure'].append(0.)


	def constant_fields(self,v2p_map=None):
		"""Constant fields.
		"""
		
		# ask constant vector
		vector = get_vector('constant velocity vector')
		
		for i,variable in enumerate(['velocity.1','velocity.2','velocity.3']):
			self.fields[variable] = [ vector[i] ] * len(self.nodes)
			
		# pressure same as outlet
		if v2p_map is not None:
			# ask constant value
			value = get_float('constant pressure value')
			
			for node in self.nodes:
				if v2p_map[node] != -1:
					self.fields['pressure'].append(value)


	def rest(self,v2p_map=None):
		"""Fluid at rest condition.
		
		v2p_map : velocity to pressure map, optional.
		
		Velocity : homogeneous.
		Pressure : homogeneous.
		"""

		# velocity same as wall
		self.wall()

		# pressure same as outlet
		if v2p_map is not None:
			self.outlet(v2p_map)


	def linear_shear(self,v2p_map):
		"""Linear shear flow condition.
		
		v2p_map : velocity to pressure map, optional.
		
		Velocity : linear shear.
		Pressure : homogeneous.
		"""

		# get shear direction and 0-offset it
		shear_direction = int(get_direction('shear')) - 1

		# get flow direction and set relevant variable
		variable = 'velocity.'+get_direction('flow')

		# get shear rate and velocity norm at origin
		rate = get_non_zero_float('shear rate : ')
		velocity_at_origin = get_float('velocity norm at origin : ')

		# initialize with rest flow
		self.rest(v2p_map)

		# set non zero component to linear shear flow
		for node in self.nodes:
			point = self.points[node]
			self.fields[variable][node] = velocity_at_origin + rate * point[shear_direction]
		

	def free_surface_on_inclined_plane(self,openings,v2p_map=None):
		""".

		openings : list of openings, should have one inlet and one outlet.
		v2p_map : velocity to pressure map, optional.
		
		Velocity : parabolic.
		Pressure : linear.
		
		Requires an outlet opening data, the maximum fluid velocity and the Reynolds number based on the maximum velocity and tube radius (equivalent to specifying viscosity).
		"""

		# get an outlet
		opening = select_opening(openings,'outlet')
		
		# get maximum velocity
		max_vel = get_non_zero_float('Maximum velocity : ')

		# ask shear direction and 0-offset it
		shear_direction = int(get_direction('shear')) - 1

		# get angle of inclinaison
		# angle = math.acos(opening.normal[shear_direction]) - math.pi / 2
		# print angle

		# get height in shear direction
		height = get_non_zero_float('Height : ')
		
		# get Reynolds number
		Reynolds = get_non_zero_float('Reynolds number : ')

		# pressure gradient coefficient
		pressure_gradient = - 2 / Reynolds

		# initialize to fluid at rest
		self.rest(v2p_map)

		# loop over nodes and prescribe poiseuille flow
		for node in self.nodes:
			# coordinates of the current point
			point = self.points[node]
			
			# position relative to opening center
			CM = point[:]
			for i in range(3):
				CM[i] -= opening.center[i]

			# algebraic distance from outlet
			distance_from_outlet = scalar_product(CM,opening.normal)

			# distance from bottom
			distance = opening.center[shear_direction] + CM[shear_direction]
			distance /= height
			
			# algebraic velocity magnitude
			vel_magn = distance * ( 2. - distance )

			# compute velocity
			for direction in [1,2,3]:
				variable = 'velocity.'+str(direction)
				value = vel_magn * opening.normal[direction-1]
				self.fields[variable][node] = value * max_vel

			# get pressure node index
			pressure_node = v2p_map[node]

			# is it a pressure node ?
			if pressure_node != -1:
				self.fields['pressure'][pressure_node] = 0. + pressure_gradient * distance_from_outlet

		
	def poiseuille(self,opening,v2p_map=None,pressure_gradient=None):
		"""Set a poiseuille flow."""

		# initialize to fluid at rest
		self.rest(v2p_map)
		
		# get maximum velocity
		max_vel = get_non_zero_float('Maximum velocity : ')

		# loop over nodes and prescribe poiseuille flow
		for node_id,node in enumerate(self.nodes):
			# node's position relative to opening center
			CM = [0]*3
			for i in range(3):
				CM[i] = self.points[node][i] - opening.center[i]

			# distance to opening plane
			oc = scalar_product(CM,opening.normal)

			# radius square normed
			radius2 = ( scalar_product(CM,CM) - pow(oc,2) ) / pow(opening.radius,2)

			# compute velocity components
			for direction in [1,2,3]:
				variable = 'velocity.'+str(direction)
				value = ( 1. - radius2 ) * opening.normal[direction-1]
				assert value <= 1., \
					"non dimensional node radius > 1 : %g"%value
				self.fields[variable][node_id] = value * max_vel


			# this part is only for initial condition
			# i.e. for domain zone wich has all nodes
			# so we don't care about array index and node index
			if v2p_map is not None:
				# get pressure node index
				pressure_node = v2p_map[node]

				# is it a pressure node ?
				if pressure_node != -1:
					self.fields['pressure'][pressure_node] = 0. + pressure_gradient * oc


	def poiseuille_straight_tube(self,openings,v2p_map):
		"""Poiseuille in a circular 3d straight tube condition.

		openings : list of openings, should have one inlet and one outlet.
		v2p_map : velocity to pressure map, optional.
		
		Velocity : parabolic.
		Pressure : linear.
		
		Requires an opening data, the maximum fluid velocity and the Reynolds number based on the maximum velocity and tube radius (equivalent to specifying viscosity).
		"""

		# get an outlet
		opening = select_opening(openings,'outlet')

		# get Reynolds number
		Reynolds = get_non_zero_float('Reynolds number : ')

		# set pressure gradient coefficient
		pressure_gradient = - 4. / Reynolds

		self.poiseuille(opening,v2p_map,pressure_gradient)


	def inlet(self,openings):
		"""Set inlet bc.
		
		openings : openings data, from which we need an inlet"""
		
		# get an inlet and copy it since we'll modify it
		opening = copy.deepcopy(select_opening(openings,'inlet'))
		
		# flow direction is opposed to inlet normal
		opening.normal = [ -component for component in opening.normal ]

		self.poiseuille(opening)
	
	
	def from_field(self,field,v2p_map=None):
		"""docstring for from_fields"""
		
		for variable in ['velocity.1','velocity.2','velocity.3']:
			for node in self.nodes:
				value = field.fields[variable][node]
				self.fields[variable].append(value) 
		
		if v2p_map is not None:
			variable = 'pressure'
			for node in self.nodes:
				pressure_node = v2p_map[node]
				if pressure_node != -1:
					value = field.fields[variable][pressure_node]
					self.fields[variable].append(value) 


	def impose(self,fields_from,v2p_map,check=True):			
		"""Impose a bc to all other zones on nodes they share.
		
		fields_from : the field that will overwrite self field
		v2p_map : velocity to pressure map
		check : hack to speed up for domain"""

		if self == fields_from: return

		for variable,field_from in fields_from.fields.iteritems():
			# just an alias
			field_to = self.fields[variable]			

			# skip if no bc for this variable
			# or for corresponding variable of imposing bc
			if len(field_from) == 0 or len(field_to) == 0: continue
			
			nodes_from = fields_from.nodes
			
			# for pressure we need to map
			if variable == 'pressure':
				nodes_from = [v2p_map[node] for node in nodes_from if v2p_map[node] != -1 ]
			
			# loop over nodes in self zone
			for idx_from,node in enumerate(nodes_from):
				# is this node in the imposing bc zone ?
				if check == True:
					if self.nodes.count(node) == 0: continue
					# get the corresponding index of the value
					idx_to = self.nodes.index(node)
				else:
					idx_to = node
				
				# impose value
				field_to[idx_to] = field_from[idx_from]
		