#!/usr/bin/env python
# encoding: utf-8

import types

class Zone:
	"""Zone class.
	
	A zone is a set of cells, it is defined by the type of cells and a connectivity. The connectivity is a list of cells local connectivity, each cell local connectivity is a list of the global indexes of the nodes that define the cell.	
	 Numbering is 0-offset.
	"""
	
	def __init__(self,cell_type):
		"""cell_type : a string which defines the type of cell in the zone, the last 2 characters are the number of nodes per cell."""
		
		# type of cell
		self.cell_type = cell_type
		
		# number of nodes per cell, for checking purposes
		self.nb_of_nodes = int(self.cell_type[-2:])
		
		# connectivity, a list of list, each cell being a list of nodes
		self.connectivity = []


	def __str__(self):
		"""Print zone info."""
		return str(len(self.connectivity)) + ' ' + self.cell_type + '\n'
		
	
	def number_of_elements(self):
		"""Returns the number of elements in the zone"""
		return len(self.connectivity)


	def append(self,nodes):
		"""Add a cell to the zone.
		
		nodes : list of connectivity of a cell"""
		
		# check it is a list
		assert isinstance(nodes,types.ListType), \
			'nodes must be a list'
		
		# check number of nodes
		assert len(nodes) == self.nb_of_nodes, \
			'Bad number of nodes : %d, should be %d' \
			%(len(nodes),self.nb_of_nodes)
		
		# add cell's nodes list to connectivity
		self.connectivity.append(nodes)
	
	
	def list_of_nodes(self):
		"""Create and return the list of nodes that belong to a zone."""
		
		# nodes are first stored in a set to avoid duplicates
		nodes = set()
		
		# build nodes set
		for cell in self.connectivity:
			for node in cell:
				nodes.add(node)
		
		# convert to list and return it
		return list(nodes)
