#!/usr/bin/python
#
# This is a script to convert RegRNA prediction output files into BED and GFF3 tracks for JBrowse.
# RegRNA2.0 - http://regrna2.mbc.nctu.edu.tw/ can be used to predict regulatory elements
# of a given RNA sequence.
# This program was written based on an regRNA prediction file in the following format:
#
#   Motif Type      Motif Name      Position        Length  Sequence
#   open reading frame (ORF)        ORF_0   132 ~ 398       267     atgca
#
# For each entry in the input file then, there should be five columns
#
# Usage:
#
# $ python3.6  561L_REGRNA_to_BED_GFF3.py inputfile outputfile
#

import sys      # this allows for sys.argv to be called
import re       # this is used to implement a split command

filename = sys.argv[1]
output = sys.argv[2]
gene_coordinate = input('What are the coordinates of your gene (e.g., chr1:2555639..2565382): ')
strand = input('What strand is your gene on (type + or -): ')
bedfile = open(output+'.bed', 'w')
gff3file = open(output+'.gff3', 'w')
with open(filename, newline='') as f:

#this is a function to check for numbers
    #def hasNumbers (inputString):
    #    return bool(re.search(r'\d', inputString))

# f.read()  <-- This would read the entire file and print it as a string
    #header = f.readlines()[0:1] # this will focus on just the header portion of the input file
    #for row in header:
    #    print(header)
    lines = f.readlines()[2:] #this will skip the header
    for row in lines:
        data = row.split('\t')
        motif_type = data[0]
        #print(motif+'\t') #this is for testing
        motif_name = data[1]
        #print(motif_name+'\t') #this is for testing
        motif_position = re.split('\~', data[2])
        motif_i = int(motif_position[0])
        #print(motif_i+'\t') #this is for testing
        motif_j = int(motif_position[1])
        #print(motif_j+'\t') #this is for testing
        motif_length = int(data [3])
        #print(motif_length+'\t') #this is for testing   chr2:27370315..27370486
        sequence = data [4]

        if "+" in strand:
            chromosome_data = re.split('\:', gene_coordinate)
            chromosome = chromosome_data[0]
            print(chromosome) #for testing
            fasta_coordinates = re.split('\.\.', chromosome_data[1])
            gene_start = int(fasta_coordinates[0])-1
            print(gene_start) #for testing
            icoordinate = str(motif_i + gene_start)
            jcoordinate = str(motif_j + (gene_start))
            print(chromosome, icoordinate, jcoordinate, motif_name, "0", strand)
            bedfile.write(chromosome+'\t'+str(int(icoordinate)-1)+'\t'+jcoordinate+'\t'+motif_name+'\t'+"0"+'\t'+strand+'\n')
            gff3file.write(chromosome+'\t'+'.'+'\t'+'sequence_attribute'+'\t'+icoordinate+'\t'+jcoordinate+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+'Motif_type='+motif_type+';'+'Name='+motif_name+';'+'Length='+str(motif_length)+';'+'Sequence='+sequence+'\n')
        if "-" in strand:
            chromosome_data = re.split('\:', gene_coordinate)
            print(chromosome_data)
            chromosome = chromosome_data[0]
            print(chromosome) #for testing
            fasta_coordinates =  re.split('\.\.', chromosome_data[1])
            print(fasta_coordinates)
            gene_start = int(fasta_coordinates[1])+1
            print(gene_start) #for testing
            jcoordinate = str(gene_start - motif_i)
            icoordinate = str(gene_start - motif_j)
            print(chromosome, icoordinate, jcoordinate, motif_name, "0", strand)
            bedfile.write(chromosome+'\t'+str(int(icoordinate)-1)+'\t'+jcoordinate+'\t'+motif_name+'\t'+"0"+'\t'+strand+'\n')
            gff3file.write(chromosome+'\t'+'.'+'\t'+'sequence_attribute'+'\t'+icoordinate+'\t'+jcoordinate+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+'Motif_type='+motif_type+';'+'Name='+motif_name+';'+'Length='+str(motif_length)+';'+'Sequence='+sequence+'\n')
