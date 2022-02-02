#!/usr/bin/env python3

# GenePicker is a program for randomly selecting a gene.
# This Python script requires an input text file.
# This file should contain a list of gene names (one per line)


# Import system-specific functions, allowing passing arguments to python

import sys

# import shuffle from random method
   
from random import shuffle


# Create a variable called file_name and get this informatio from sys.argv.
# sys.argv is a list of arguments that are passed to the script via command-line. 

file_name = sys.argv[1]             # FYI argv[0] is the script name    

# Create a file object "gene_file" and place each line/gene into a list called genes

gene_file = open(file_name, "r")    # "r" indicates read-only
genes = gene_file.readlines()       # the readlines function separates each line
gene_file.close()                   # close gene_file to save memory

# Randomize the genes by shuffling the array with shuffle function

shuffle(genes)

# Prompt user to enter names

user_name = input("Please enter your first name: ")           

# Pick gene from randomized list using the length of the user's name

print (user_name, "your gene is:", genes[len(user_name)])



