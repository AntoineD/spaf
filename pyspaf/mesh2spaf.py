#!/usr/bin/env python

"""
Mesh converter.

Converts mesh from a file format to another. Currently it can read mesh from gmsh and spaf formats, and write to spaf format.

"""

__author__ = "Antoine Dechaume"
__author_email__ = "dechaume@ualberta.ca"
__version__ = "r22"

import sys,os
from case import *	
from optparse import OptionParser

if __name__ == '__main__':

	# set up command line parser
	usage = "usage: %prog mesh_file [options]"
	version = "%prog "+__version__
	parser = OptionParser(usage=usage, version=version)

	# output option
	# it is the name without extension of the output file
	parser.add_option(
		"-o",
		type="string",
		action="store",
		dest="output",
		help="write to FILE.case",
		metavar="FILE")

	# data format option
	# set to ascii when invoqued
	parser.add_option(
		"-a","--ascii",
		action="store_const",
		const="ascii",
		dest="data_format",
		default="binary",
		help="write in ascii format, default is binary" )
		
	# parse command line
	(options,args) = parser.parse_args()


	# process arguments
	# check we have the input mesh file
	if len(args) != 1:
		parser.error("You have to give a mesh file as argument.")
    
	input_file = args[0]
	
	# output file
	if options.output:
		output_file = options.output
	else:
		output_file = os.path.basename(os.path.splitext(input_file)[0])
		
	# we'll write to spaf format
	output_file += '.case'


	# create case instance
	case = Case(input_file)

	# write to file
	case.to_file(options.data_format,output_file)
	