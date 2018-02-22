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
import pyBigWig



filename = sys.argv[1] # this should be the output of a z-score analysis in tab-delimited format
#output = sys.argv[2]

gene_coordinate = input('What are the coordinates of your gene (e.g., chr1:2555639..2565382): ')
#chromosome = str(input("What chromsome is your gene on? (e.g., chr1): "))
#genomic_start = int(input("What is the starting coordinate of your gene sequence? (e.g., 2555639; without commas): " ))
step_size = int(input("What is the step size for your data? (e.g., 1):  "))
window_size = int(input("What is the window size for your data? (e.g., 120): "))
strand = input("What strand is your gene on (+ or -)?: ")

chromosome_data = re.split('\:', gene_coordinate)
chromosome = chromosome_data[0]
genomic_coordinates = re.split('\.\.', chromosome_data[1])
genomic_start = genomic_coordinates[0]
genomic_end = genomic_coordinates[1]

MFE_wig = pyBigWig.open(filename+'.strand'+strand+'_MFE_'+'.bw', 'w')
zscore_wig = pyBigWig.open(filename+'.strand'+strand+'_zscore_'+'.bw', 'w')
pscore_wig = pyBigWig.open(filename+'.strand'+strand+'_pscore_'+'.bw', 'w')
ED_wig = pyBigWig.open(filename+'.strand'+strand+'_Ed_'+'.bw', 'w')
fMFE_wig = pyBigWig.open(filename+'.strand'+strand+'_fMFE_'+'.bw', 'w')
gff3file = open(filename+'.strand'+strand+'.gff3', 'w')

MFE_wig.addHeader([("chr1",248956422),("chr2",242193529),("chr3",198295559),("chr4",190214555),("chr5",181538259),("chr6",170805979),("chr7",159345973),("chr8",145138636),("chr9",138394717),("chr10",133797422),("chr11",135086622),("chr12",133851895),("chr13",115169878),("chr14",107349540),("chr15",102531392),("chr16",90354753),("chr17",107349540),("chr18",78077248),("chr19",59128983),("chr20",63025520),("chr21",48129895),("chr22",51304566),("chrX",155270560),("chrY",59373566)])
zscore_wig.addHeader([("chr1",248956422),("chr2",242193529),("chr3",198295559),("chr4",190214555),("chr5",181538259),("chr6",170805979),("chr7",159345973),("chr8",145138636),("chr9",138394717),("chr10",133797422),("chr11",135086622),("chr12",133851895),("chr13",115169878),("chr14",107349540),("chr15",102531392),("chr16",90354753),("chr17",107349540),("chr18",78077248),("chr19",59128983),("chr20",63025520),("chr21",48129895),("chr22",51304566),("chrX",155270560),("chrY",59373566)])
pscore_wig.addHeader([("chr1",248956422),("chr2",242193529),("chr3",198295559),("chr4",190214555),("chr5",181538259),("chr6",170805979),("chr7",159345973),("chr8",145138636),("chr9",138394717),("chr10",133797422),("chr11",135086622),("chr12",133851895),("chr13",115169878),("chr14",107349540),("chr15",102531392),("chr16",90354753),("chr17",107349540),("chr18",78077248),("chr19",59128983),("chr20",63025520),("chr21",48129895),("chr22",51304566),("chrX",155270560),("chrY",59373566)])
ED_wig.addHeader([("chr1",248956422),("chr2",242193529),("chr3",198295559),("chr4",190214555),("chr5",181538259),("chr6",170805979),("chr7",159345973),("chr8",145138636),("chr9",138394717),("chr10",133797422),("chr11",135086622),("chr12",133851895),("chr13",115169878),("chr14",107349540),("chr15",102531392),("chr16",90354753),("chr17",107349540),("chr18",78077248),("chr19",59128983),("chr20",63025520),("chr21",48129895),("chr22",51304566),("chrX",155270560),("chrY",59373566)])
fMFE_wig.addHeader([("chr1",248956422),("chr2",242193529),("chr3",198295559),("chr4",190214555),("chr5",181538259),("chr6",170805979),("chr7",159345973),("chr8",145138636),("chr9",138394717),("chr10",133797422),("chr11",135086622),("chr12",133851895),("chr13",115169878),("chr14",107349540),("chr15",102531392),("chr16",90354753),("chr17",107349540),("chr18",78077248),("chr19",59128983),("chr20",63025520),("chr21",48129895),("chr22",51304566),("chrX",155270560),("chrY",59373566)])

MFE_list = []
zscore_list = []
pscore_list = []
ED_list = []
fMFE_list = []
length = (int(genomic_end) - int(genomic_start)) + 1
#print(length)


if strand == "+":
    with open(filename, 'r') as g:
        glines = g.readlines()[1:]
        #print(row)
        for row in glines:
            gdata = row.split('\t') # this splits each row based on "tab"
            #print(data)
            #print(gdata)
            #print(data)
            icoordinate = int(int(gdata[0])+int(genomic_start))
            jcoordinate = int(int(gdata[1])+int(genomic_start))
            gMFE = gdata[2]
            gzscore = gdata[3]
            if gzscore == "Undef":
                gzscore = 00000
            gpvalue = gdata[4]
            gED = gdata[5]
            gfMFE = gdata[6]
            gsequence = gdata[7]
            gfold = gdata[8]
            gcentroid = gdata[9]
            fA = gdata[10]
            fG = gdata[11]
            fC = gdata[12]
            fU = gdata[13]
            #icoordinate = int(chromosome_data[1])
            #jcoordinate = str(icoordinate+((len(str(motif)))-1))
            #print(icoordinate)
            #print(jcoordinate)
            #print(len(motif))
            #print(motif)
            #print(match)
            #print(chromosome, icoordinate, jcoordinate, protein, "0", strand)
            #bedfile.write(chromosome+'\t'+str(icoordinate-1)+'\t'+jcoordinate+'\t'+protein+'\t'+"0"+'\t'+strand+'\n')
            gff3file.write(chromosome+'\t'+'.'+'\t'+'sequence_attribute'+'\t'+str(icoordinate)+'\t'+str(jcoordinate)+'\t'+'.'+'\t'+strand+'\t'+'.\t'+'MFE='+str(gMFE)+';'+'Z-score='+str(gzscore)+';'+'P-value='+str(gpvalue)+';'+'EnsDiv='+str(gED)+';'+'fMFE='+str(gfMFE)+';'+'Sequence='+gsequence+';'+'MFE_Fold='+gfold+';'+'Centroid='+gcentroid+'\n')

if strand == "-":
    with open(filename, 'r') as g:
        glines = g.readlines()[1:]
        for row in glines:
            gdata = row.split('\t') # this splits each row based on "tab"
            #print(data)
            #print(gdata)
            #print(data)
            icoordinate = int(int(genomic_start)+(int(length)-int(gdata[1])))
            jcoordinate = int(int(genomic_start)+(int(length)-int(gdata[0])))
            gMFE = gdata[2]
            gzscore = gdata[3]
            if gzscore == "Undef":
                gzscore = 00000
            gpvalue = gdata[4]
            gED = gdata[5]
            gfMFE = gdata[6]
            gsequence = gdata[7]
            gfold = gdata[8]
            gcentroid = gdata[9]
            fA = gdata[10]
            fG = gdata[11]
            fC = gdata[12]
            fU = gdata[13]
            #icoordinate = int(chromosome_data[1])
            #jcoordinate = str(icoordinate+((len(str(motif)))-1))
            #print(icoordinate)
            #print(jcoordinate)
            #print(len(motif))
            #print(motif)
            #print(match)
            #print(chromosome, icoordinate, jcoordinate, protein, "0", strand)
            #bedfile.write(chromosome+'\t'+str(icoordinate-1)+'\t'+jcoordinate+'\t'+protein+'\t'+"0"+'\t'+strand+'\n')
            gff3file.write(chromosome+'\t'+'.'+'\t'+'sequence_attribute'+'\t'+str(icoordinate)+'\t'+str(jcoordinate)+'\t'+'.'+'\t'+strand+'\t'+'.\t'+'MFE='+str(gMFE)+';'+'Z-score='+str(gzscore)+';'+'P-value='+str(gpvalue)+';'+'EnsDiv='+str(gED)+';'+'fMFE='+str(gfMFE)+';'+'Sequence='+gsequence+';'+'MFE_Fold='+gfold+';'+'Centroid='+gcentroid+'\n')



with open(filename, 'r') as f:
    if strand == "+":
        genomic_start = int(genomic_start)
        lines = f.readlines()[1:]
        for row in lines:
            data = row.split('\t') # this splits each row based on "tab
            #print(len(data))
            i = data[0]
            j = data[1]
            genomic_end = int(genomic_start)+int(step_size)
            MFE = float(data[2])
            MFE_list.append(MFE)
            if data[3] == "Undef":
                zscore = float(00000)
            else:
                zscore = float(data[3])
            zscore_list.append(zscore)
            pscore = float(data[4])
            pscore_list.append(pscore)
            ED = float(data[5])
            ED_list.append(ED)
            fMFE = float(data[6])
            fMFE_list.append(fMFE)

        MFE_wig.addEntries(chromosome, genomic_start,  values=MFE_list, span=step_size, step=step_size)
        zscore_wig.addEntries(chromosome, genomic_start,  values=zscore_list, span=step_size, step=step_size)
        pscore_wig.addEntries(chromosome, genomic_start,  values=pscore_list, span=step_size, step=step_size)
        ED_wig.addEntries(chromosome, genomic_start,  values=ED_list, span=step_size, step=step_size)
        fMFE_wig.addEntries(chromosome, genomic_start,  values=fMFE_list, span=step_size, step=step_size)

    if strand == "-":
        lines = reversed(open(filename).readlines()[1:])
        start = genomic_start
        genomic_start = int(start) + int(window_size)
        for row in lines:
                data = row.split('\t') # this splits each row based on "tab
                #print(len(data))
                #print(row)
                i = data[0]
                j = data[1]
                #genomic_start = int(genomic_start)+int(window_size)
                MFE = float(data[2])
                MFE_list.append(MFE)
                if data[3] == "Undef":
                    zscore = float(00000)
                else:
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
