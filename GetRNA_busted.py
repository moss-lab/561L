#!/usr/bin/env python3

# GetRNA is a program for randomly selecting an RNA sequence from a file.
# This Python script requires an input text file.
# This file should contain a set of RNA sequences, one per line


# Import system-specific functions, allowing passing arguments to python

import sys

# import shuffle from random method
   
from random import shuffle


# Create a variable called file_name and get this informatio from sys.argv.
# sys.argv is a list of arguments that are passed to the script via command-line.
# FYI argv[0] is the script name  and the subsequent argv's are the arguments

file_name = sys.argv[1]              
user_name = sys.argv[2]                  

# Create a file object "gene_file" and place each line/gene into a list called genes

rna_file = open(file_name, "r")    # "r" indicates read-only
rnas = rna_file.readlines()       # the readlines function separates each line
rna_file.close()                   # close gene_file to save memory

# Randomize the rnas by shuffling the array with shuffle function

shuffle(rnas)

# Pick RNA from randomized list using the length of the user's name

write (">",user_name, "RNA\n", rnas[len(user_name)])
