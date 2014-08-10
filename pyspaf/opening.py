#!/usr/bin/env python
# encoding: utf-8
"""
opening.py

Created by Antoine Dechaume on 2008-04-28.
Copyright (c) 2008 __MyCompanyName__. All rights reserved.

This file contains code for managing openings geometric data, i.e. inlet and outlet.
"""

from tools import *

class Opening:
	"""Class for opening data.
	
	Openings are inlet or outlet of a domain. The Lagrange characteristic method for computing convection term requires to check if a characteristic crosses one of the opening. To do so we need to know the geometry of the disc which includes the corresponding inlet or outlet, such as its radius, center and normal vector pointing outside the physical domain."""
	
	def __init__(self,kind):
		"""kind : kind of opening, 'inlet' or 'outlet'"""
		
		# check kind
		assert kind in ['inlet','outlet'], \
			"Bad opening kind, should be 'inlet' or 'outlet' : %s"%kind

		self.kind = kind
		
		# initialize other attribute
		# radius of the disc
		self.radius = 0.
		# coordinates of the center of the disc
		self.center = []
		# components of the normal of the disc, pointing outside the domain
		self.normal = []
		
		# read opening data from user
		self.get_opening()


	def get_opening(self):
		"""Read an opening data from user."""

		# get radius
		# make sure it a strictly positive real number
		while True:
			radius = raw_input('\tradius : ')
			if string_is_float(radius):
				radius = float(radius)
				if radius > 0: break

		self.radius = radius

		# get center point
		self.center = get_vector('center')
		
		# get normal and make sure it is not 0
		while True:
			self.normal = get_vector('normal')

			if norm(self.normal) != 0: break
			print 'normal is null, try again !'
		
		# make normal normed
		self.normal = [ coord / norm(self.normal) for coord in self.normal ]
		
		print ''


	def to_file(self,file_id):
		"""Write an opening to a file.
		
		file_id : file index to write to"""
		
		file_id.write('section = '+self.kind+'\n')
		file_id.write('\tradius = '+str(self.radius)+'\n')
		file_id.write('\tcenter = '+str(self.center[0])+' '+ \
			str(self.center[1])+' '+str(self.center[2])+'\n')
		file_id.write('\tnormal = '+str(self.normal[0])+' '+ \
			str(self.normal[1])+' '+str(self.normal[2])+'\n')
		file_id.write('section = end'+'\n')
		