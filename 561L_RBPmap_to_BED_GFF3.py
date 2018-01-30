#this is a script to convert RBPmap prediction files into bed tracks for JBrowse
# This assumes that the RBPMap prediction file contains the genomic coordinate descriptions
# for each RNA binding site, in the following format:
# Protein: BRUNOL4(Hs/Mm)
# Sequence Position	Genomic Coordinate	Motif     	K-mer     	Z-score	P-value
# 2745             	chr1:2558383      	kgugukk   	agugugu   	1.784

import sys
import re

filename = sys.argv[1]
output = sys.argv[2]
strand = input('What strand is your gene on (type + or -): ')
bedfile = open(output+'.bed', 'w')
gff3file = open(output+'.gff3', 'w')
with open(filename, newline='') as f:

#this is a function to check for numbers
    #def hasNumbers (inputString):
    #    return bool(re.search(r'\d', inputString))

# f.read()  <-- This would read the entire file and print it as a string
    #header = f.readlines()[0:9] # this will focus on just the header portion of the input file
    #for row in header:
    #    print(header)
    lines = f.readlines()[9:] #this will skip the header
    for row in lines:
        if 'Protein' in row:
                #print("true")
                protein_header = row.split()
                #print (protein_header)
                protein_data = re.split('\(', protein_header[1])
                protein = protein_data[0]
                print(protein)
        if 'chr' in row:
                #print(row)
                data = row.split()
                motif = data[2]
                match = data[3]
                zscore = data[4]
                pvalue = data [5]
                chromosome_data = re.split('\:', data[1])
                chromosome = chromosome_data[0]
                icoordinate = int(chromosome_data[1])
                jcoordinate = str(icoordinate+((len(str(motif)))-1))
                #print(icoordinate)
                #print(jcoordinate)
                #print(len(motif))
                #print(motif)
                #print(match)
                print(chromosome, icoordinate, jcoordinate, protein, "0", strand)
                bedfile.write(chromosome+'\t'+str(icoordinate-1)+'\t'+jcoordinate+'\t'+protein+'\t'+"0"+'\t'+strand+'\n')
                gff3file.write(chromosome+'\t'+'.'+'\t'+'sequence_attribute'+'\t'+str(icoordinate)+'\t'+jcoordinate+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+'Name='+protein+';'+'Binding_Motif='+motif+';'+'Binding_Match='+match+';'+'Zscore='+zscore+';'+'Pvalue='+pvalue+'\n')
