#!/usr/bin/python
#
# This program converts RBPmap prediction files into BED and GFF3 tracks for JBrowse
# This assumes that the RBPMap prediction file contains accurate genomic coordinates
# for each RNA binding site. The input file should be in the following format:
#
#   Protein: BRUNOL4(Hs/Mm)
#   Sequence Position	Genomic Coordinate	Motif     	K-mer     	Z-score	P-value
#   2745             	chr1:2558383      	kgugukk   	agugugu   	1.784
#
# For each entry in the input file then, there should be a title containing the
# name of the binding protein followed by a six-column entry containing the
# internal coordinates of the binding site, the genomic coordinates, the "consensus"
# motif sequence for the protein, the matching sequence on the given RNA sequences,
# followed by the z-score and p-value of match.
#
# Usage:
#
# $ python3.6  561L_REGRNA_to_BED_GFF3.py inputfile outputfile


import sys
import re

filename = sys.argv[1]
output = sys.argv[2]
strand = input('What strand is your gene on (type + or -): ')
bedfile = open(output+'.bed', 'w')
gff3file = open(output+'.gff3', 'w')
with open(filename, newline='') as f:

    lines = f.readlines()[9:] #this will skip the input file header
    for row in lines:
        if 'Protein' in row: #this singles out lines with "Protein" in them
                protein_header = row.split()
                protein_data = re.split('\(', protein_header[1]) # this splits the protein line based on this character "("
                protein = protein_data[0] # this defines the protein name as the first portion of the protein
                print(protein)


        if 'chr' in row:

                #print(row)
                data = row.split()
                motif = data[2]
                match = data[3]
                zscore = data[4]
                pvalue = data [5]
                chromosome_data = re.split('\:', data[1]) # the second column actually has two elements we need to split at the colon chr1:2558383
                chromosome = chromosome_data[0]
                icoordinate = int(chromosome_data[1])
                jcoordinate = str(icoordinate+((len(str(motif)))-1))
                #print(icoordinate)
                #print(jcoordinate)
                #print(len(motif))
                #print(motif)
                #print(match)
                if "+" in strand:
                    print("POSITIVE"+chromosome, icoordinate, jcoordinate, protein, "0", strand)
                    bedfile.write(chromosome+'\t'+str(icoordinate-1)+'\t'+jcoordinate+'\t'+protein+'\t'+"0"+'\t'+strand+'\n')
                    gff3file.write(chromosome+'\t'+'.'+'\t'+'sequence_attribute'+'\t'+str(icoordinate)+'\t'+jcoordinate+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+'Name='+protein+';'+'Binding_Motif='+motif+';'+'Binding_Match='+match+';'+'Zscore='+zscore+';'+'Pvalue='+pvalue+'\n')
                if "-" in strand:
                    print("NEGATIVE"+chromosome, icoordinate, jcoordinate, protein, "0", strand)
                    bedfile.write(chromosome+'\t'+str(icoordinate-(len(str(motif))))+'\t'+str(icoordinate)+'\t'+protein+'\t'+"0"+'\t'+strand+'\n')
                    gff3file.write(chromosome+'\t'+'.'+'\t'+'sequence_attribute'+'\t'+str(icoordinate-(len(str(motif))-1))+'\t'+str(icoordinate)+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+'Name='+protein+';'+'Binding_Motif='+motif+';'+'Binding_Match='+match+';'+'Zscore='+zscore+';'+'Pvalue='+pvalue+'\n')
