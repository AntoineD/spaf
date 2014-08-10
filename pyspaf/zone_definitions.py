#!/usr/bin/env python
# encoding: utf-8
"""
zone_definitions.py

Created by Antoine Dechaume on 2008-05-22.
Copyright (c) 2008 __MyCompanyName__. All rights reserved.
"""

# This dict defines zones, by a name (fisrt key), then each key is associated to another dict, wich defines the the type of cell in the zone, the kind of boundary conditions for the velocity and the pressure. A bool defines whether a zone has prescirbed nodes or not.

# domain bc are set to free for convenience but will have no impact as a non free bc overrides a free one.

zone_definitions = { \
'domain' : { 'cell_type' : 'tetrahedron_10', 'velocity' : False , 'pressure' : False }, \
'wall'   : { 'cell_type' : 'triangle_06'   , 'velocity' : True  , 'pressure' : False }, \
'inlet'  : { 'cell_type' : 'triangle_06'   , 'velocity' : True  , 'pressure' : False }, \
'outlet' : { 'cell_type' : 'triangle_06'   , 'velocity' : False , 'pressure' : True  }, \
'shear'  : { 'cell_type' : 'triangle_06'   , 'velocity' : True  , 'pressure' : False }, \
'z'      : { 'cell_type' : 'triangle_06'   , 'velocity' : True  , 'pressure' : False }, \
}


def zone_name_is_valid(zone_name):
	"""Check a zone name is in the definitions, return True if it is, False otherwise."""
	
	if zone_name not in zone_definitions:
		print 'Bad zone name : ',zone_name
		print 'It should be one of the following :\n\t'
		for name in zone_definitions:
			print name,
		print ''
		return False
	
	return True


def get_zone_cell_type(zone_name):
	"""Return the type of cell.
	
	zone_name : name of the zone"""
	
	assert zone_name_is_valid(zone_name)
	
	return zone_definitions[zone_name]['cell_type']


def zone_is_prescribed(zone_name,variable):
	"""Return true if a given variable and zone is prescribed.
	
	zone_name : name of the zone
	variable : 'velocity' or 'pressure'"""

	assert zone_name_is_valid(zone_name)
	assert variable in ['velocity','pressure']

	return zone_definitions[zone_name][variable]
	