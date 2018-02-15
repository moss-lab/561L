#!usr/bin/python3.6 -w

#This program will take a the output from the ZScore_calculation.pl program append
# convert it into bigwig tracks.
#
# Usage:
#
# $ python thisScrip.py inputfile outfile
#

import sys # this will allow for the use of system argument inputs
import re
import numpy as np
import pyBigWig
import RNA

filename = sys.argv[1] # this should be the output of a z-score analysis in tab-delimited format
#output = sys.argv[2]
MFE_wig = pyBigWig.open('MFE_'+filename+'.bw', 'w')
zscore_wig = pyBigWig.open('zscore_'+filename+'.bw', 'w')
pscore_wig = pyBigWig.open('pscore_'+filename+'.bw', 'w')
ED_wig = pyBigWig.open('ED_'+filename+'.bw', 'w')
fMFE_wig = pyBigWig.open('fMFE_'+filename+'.bw', 'w')

chromosome = str(raw_input("What chromsome is your gene on? (e.g., chr1): "))
genomic_start = int(raw_input("What is the starting coordinate of your gene sequence? (e.g., 2555639; without commas): " ))
step_size = int(raw_input("What is the step size for your data? (e.g. 1): " ))

MFE_wig.addHeader([("chr1",	248956422),("chr2", 242193529),("chr3", 198295559),("chr4", 190214555),("chr5", 181538259),("chr6", 170805979),("chr7", 159345973),("chr8", 145138636),("chr9",138394717),("chr10",133797422),("chr10",133797422),("chr11",135086622)])
zscore_wig.addHeader([("chr1",	248956422),("chr2", 242193529),("chr3", 198295559),("chr4", 190214555),("chr5", 181538259),("chr6", 170805979),("chr7", 159345973),("chr8", 145138636),("chr9",138394717),("chr10",133797422),("chr10",133797422),("chr11",135086622)])
pscore_wig.addHeader([("chr1",	248956422),("chr2", 242193529),("chr3", 198295559),("chr4", 190214555),("chr5", 181538259),("chr6", 170805979),("chr7", 159345973),("chr8", 145138636),("chr9",138394717),("chr10",133797422),("chr10",133797422),("chr11",135086622)])
ED_wig.addHeader([("chr1",	248956422),("chr2", 242193529),("chr3", 198295559),("chr4", 190214555),("chr5", 181538259),("chr6", 170805979),("chr7", 159345973),("chr8", 145138636),("chr9",138394717),("chr10",133797422),("chr10",133797422),("chr11",135086622)])
fMFE_wig.addHeader([("chr1",	248956422),("chr2", 242193529),("chr3", 198295559),("chr4", 190214555),("chr5", 181538259),("chr6", 170805979),("chr7", 159345973),("chr8", 145138636),("chr9",138394717),("chr10",133797422),("chr10",133797422),("chr11",135086622)])
#print(chromosome)


MFE_list = []
zscore_list = []
pscore_list = []
ED_list = []
fMFE_list = []

with open(filename, 'r') as f:
    lines = f.readlines()[1:]
    for row in lines:
        data = row.split('\t') # this splits each row based on "tab
        print(len(data))
        if (len(data) == 14) and ('\-' in data[5]):
                print(data)
                data.remove(data[5])

                #corrected_row = ('\t'.join(data))
                print(data)
                #print(row)
                i = data[0]
                j = data[1]
                genomic_end = int(genomic_start)+int(step_size)
                MFE = float(data[2])
                MFE_list.append(MFE)
                zscore = float(data[3])
                zscore_list.append(zscore)
                pscore = float(data[4])
                pscore_list.append(pscore)
                ED = float(data[5])
                ED_list.append(ED)
                fMFE = float(data[6])
                fMFE_list.append(fMFE)

        if len(data) == 14:
            #print(row)
            i = data[0]
            j = data[1]
            genomic_end = int(genomic_start)+int(step_size)
            MFE = float(data[2])
            MFE_list.append(MFE)
            zscore = float(data[3])
            zscore_list.append(zscore)
            pscore = float(data[4])
            pscore_list.append(pscore)
            ED = float(data[5])
            ED_list.append(ED)
            fMFE = float(data[6])
            fMFE_list.append(fMFE)




#print(MFE_list)
#print(chromosome)
#print(step_size)
MFE_wig.addEntries(chromosome, genomic_start,  values=MFE_list, span=step_size, step=step_size)
zscore_wig.addEntries(chromosome, genomic_start,  values=zscore_list, span=step_size, step=step_size)
pscore_wig.addEntries(chromosome, genomic_start,  values=pscore_list, span=step_size, step=step_size)
ED_wig.addEntries(chromosome, genomic_start,  values=ED_list, span=step_size, step=step_size)
fMFE_wig.addEntries(chromosome, genomic_start,  values=fMFE_list, span=step_size, step=step_size)

MFE_wig.close()
zscore_wig.close()
pscore_wig.close()
ED_wig.close()
fMFE_wig.close()
