#!/usr/bin/env python3

# This script simulates SHAPE data for an RNA input in fasta format
# Import system-specific functions, allowing passing arguments to python

import sys

# import random method

import random

infile = sys.argv[1]              

fasta_file = open(infile, "r")    								# "r" indicates read-only
fasta = fasta_file.readlines()       							# the readlines function separates each line
fasta_file.close()                   							# close gene_file to save memory

seq = fasta[1]
counter = 1														# Keep track of nt number
	
	
for nt in seq:													# Loop through each nt from seq 
	simulated_reactivity = (random.randint(0,1000))/float(1000)	# Simulate SHAPE data 0-1 value
	output = "%s\t%s\t%s"%(counter,nt,simulated_reactivity)
	print(output)												# Print constraint format data for RNAfold
	counter += 1												# Increment counter 
	
