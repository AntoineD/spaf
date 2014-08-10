#!/usr/bin/env python
# encoding: utf-8

import os,mesh
from zone_definitions import *

import data_file

# Free nodes come first, this order has to be matched by spaf !

class Mesh(mesh.Mesh):
	"""Mesh superclass with nodes reordering according to bc.
	
	This class inherits the what's defined in mesh module. It adds 2 attributes to base class.
	"""

	def __init__(self,file_name=None):
		"""file_name : name of the input mesh file."""
		
		# create an instance of the base class
		mesh.Mesh.__init__(self,file_name)

		# new attributes to mesh class :
		self.number_of_prescribed_velocity_nodes = 0
		self.number_of_prescribed_pressure_nodes = 0
		
		# now do the reordering
		self.reorder_velocity()
		self.reorder_pressure()
		

	def __str__(self):
		"""Print mesh informations.
		
		Appends info from reordering to base class info."""

		string  = '\n\tvelocity prescribed nodes : '
		string += str(self.number_of_prescribed_velocity_nodes)+'\n'
		string += '\tpressure prescribed nodes : '
		string += str(self.number_of_prescribed_pressure_nodes)+'\n'

		return mesh.Mesh.__str__(self)+string


	def reorder_velocity(self):
		"""Reorder velocity nodes numbering.

		Here the nodes are reordered according to whether they are prescribed or not. The nodes are then reordered such that free nodes come first."""


		print 'Reordering velocity nodes\n'

		# make a set of prescribed velocity nodes as we do not want duplicates
		prescribed_nodes = set()

		# create it
		for zone in self.zones:
			if zone_is_prescribed(zone,'velocity'):
				for cell in self.zones[zone].connectivity:
					for node in cell:
						prescribed_nodes.add(node)
		
		# make it a list and sort it, this is just for having a well defined order, for debugging purpose only
		prescribed_nodes = list(prescribed_nodes)
		prescribed_nodes.sort()
		
		self.number_of_prescribed_velocity_nodes = len(prescribed_nodes)

		# initialize old to new index permutation mapping list to None, in order to know which nodes have been processed by each of the renumbering steps
		old_to_new_node_index = [None] * self.number_of_velocity_nodes

		# iterator for new node indexes, here we start from the first prescribed node
		new_node_id = self.number_of_velocity_nodes - self.number_of_prescribed_velocity_nodes

		# process prescribed nodes
		for old_node_id in prescribed_nodes:
			old_to_new_node_index[old_node_id] = new_node_id
			new_node_id += 1

		# then free nodes
		new_node_id = 0
		for node in range(self.number_of_velocity_nodes):
			if old_to_new_node_index[node] is None:
				old_to_new_node_index[node] = new_node_id
				new_node_id += 1

		# copy the list of points with old order
		old_order = self.points[:]

		# reorder geometric table
		for node_id,node in enumerate(old_to_new_node_index):
			self.points[node] = old_order[node_id]

		# reorder connectivity tables for all zones
		for zone in self.zones.values():
			for cell in zone.connectivity:
				for local_node,global_node in enumerate(cell):
					cell[local_node] = old_to_new_node_index[global_node]



	def reorder_pressure(self):
		"""Reorder pressure nodes numbering.
		
		The velocity to pressure map is recomputed, firstly to account for velocity nodes renumbering, and secondly to have free pressure nodes indexes first, with increasing index as for velocity, i.e. if index v_node_1 > index v_node_2 then index p_node_1 > index p_node_2. This way the pressure operators will have successive rows. """
		
		print 'Reordering pressure nodes\n'
		
		# make a set of velocity nodes that are prescribed pressure nodes
		prescribed_nodes = set()

		# build this set, by selecting velocity nodes that are pressure nodes (vertex)
		for zone in self.zones:
			if zone_is_prescribed(zone,'pressure'):
				for cell in self.zones[zone].connectivity:
					for node in cell[0:3]:
						prescribed_nodes.add(node)

		# in case there is no dirichlet bc for the pressure
		# we have to set one pressure node to avoid ill-posedness
		# this choise is not important, we choose the first pressure node
		# of the first cell of the domain zone
		if len(prescribed_nodes) == 0:
			prescribed_nodes.add(self.zones['domain'].connectivity[0][0])

		self.number_of_prescribed_pressure_nodes = len(prescribed_nodes)

		# re-create the velocity to pressure map
		# as velocity node numbering has changed if we have called reorder_velocity
		self.get_v2p_map()

		# Renumber the map
		# free nodes come first, the map has to be in increasing order.
		# The 2 following variables are counters.
		free_node_id = -1
		prescribed_node_id = self.number_of_pressure_nodes - self.number_of_prescribed_pressure_nodes - 1

		# we loop over all velocity nodes, then check whether the current one
		# is a pressure node, then check whether it is prescribed for pressure
		# and set index accordingly
		for v_node,p_node in enumerate(self.v2p_map):
			if p_node != -1:
				if v_node in prescribed_nodes:
					prescribed_node_id += 1
					self.v2p_map[v_node] = prescribed_node_id
				else:
					free_node_id += 1
					self.v2p_map[v_node] = free_node_id
