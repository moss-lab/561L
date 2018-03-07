#!usr/bin/python3.6 -w

#This program will take a the output from the ZScore_calculation.pl program append
# convert it into bigwig tracks.
#
# Usage:
#
# $ python3.6 thisScrip.py inputfile
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

# Create and open output files for writing
MFE_wig = pyBigWig.open(filename+'.strand('+strand+')_MFE.bw', 'w')
zscore_wig = pyBigWig.open(filename+'.strand('+strand+')_zscore.bw', 'w')
pscore_wig = pyBigWig.open(filename+'.strand('+strand+')_pvalue.bw', 'w')
ED_wig = pyBigWig.open(filename+'.strand('+strand+')_Ed.bw', 'w')
fMFE_wig = pyBigWig.open(filename+'.strand('+strand+')_fMFE.bw', 'w')
gff3file = open(filename+'.strand('+strand+').gff3', 'w')
#corrected_file = open(filename+'.corrected.txt', 'w')

# Write header for corrected file:
#corrected_file.write("i\tj\tMFE\trandomMFE\tZscore\tPscore\tED\tfMFE\tSequence\tFold\tCentroid\t#A\t#G\t#C\t#U\n")


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

with open(filename, 'r') as g:
    if strand == "+": # Generating GFF3 file for forward strand
        glines = g.readlines()[1:]
        #print(row)
        for row in glines:
            if not row.strip():
                continue
            else:
                gdata = row.split('\t') # this splits each row based on "tab"
                if ((len(str(gdata[8])) == window_size) and
                (float(gdata[3]) != 0) and
                (("A" or "G" or "C" or "U") in gdata[8])) or (len(str(gdata[9])) == window_size and
                (float(gdata[4]) == 0) and
                (("A" or "G" or "C" or "U") in gdata[9])):
                    print("Sequence in column 9")
                    #float(gdata[7])
                    gdata = row.split('\t') # this splits each row based on "tab"
                    if len(gdata) > 15:
                        print("Errors found in file")
                        print("Error in column six:", gdata)
                        gdata.remove(gdata[5])
                        #corrected_row = ('\t'.join(data))
                        print("Error removed:", gdata)
                        #print(row)
                        #icoordinate = int((int(data[0])+int(genomic_start)-1))
                        #jcoordinate = int(int(data[1])+int(genomic_start))
                        icoordinate = (int(gdata[0])-1)+int(genomic_start)
                        jcoordinate = (int(gdata[1]))+int(genomic_start)
                        window_size = jcoordinate - icoordinate
                    #       if strand == '-1':
                    #           icoordinate = int(int(genomic_start)+(int(length)-int(data[1])))
                    #           jcoordinate = int(int(genomic_start)+(int(length)-int(data[0])))
                        MFE = float(gdata[2])
                        rand_MFE = float(gdata[3])
                        zscore = gdata[4]
                        if zscore == "Undef":
                            zscore = float(00000)
                        else:
                            zscore = float(gdata[4])
                        pvalue = gdata[5]
                        #try:
                        #    pvalue = float(gdata[5])
                        #except ValueError:
                        #    print(str(gdata[5]))
                        #    pvalue =float(0)
                        gED = float(gdata[6])
                        gfMFE = float(gdata[7])
                        gwindow_sequence = gdata[8]
                        gfold = gdata[9]
                        gcentroid = gdata[10]
                        gfA = gdata[11]
                        gfG = gdata[12]
                        gfC = gdata[13]
                        gfU = gdata[14]
                        #corrected_file.write(gdata[0]+'\t'+gdata[1]+'\t'+gdata[2]+'\t'+gdata[3]+'\t'+gdata[4]+'\t'+gdata[5]+'\t'+gdata[6]+'\t'+gdata[7]+'\t'+gdata[8]+'\t'+gdata[9]+'\t'+gdata[10]+'\t'+gdata[11]+'\t'+gdata[12]+'\t'+gdata[13]+'\t'+gdata[14])
                        gff3file.write(chromosome+'\t'+'.'+'\t'+'sequence_attribute'+'\t'+str(icoordinate)+'\t'+str(jcoordinate)+'\t'+'.'+'\t'+strand+'\t'+'.\t'+'MFE='+str(gMFE)+';'+'Z-score='+str(gzscore)+';'+'P-value='+str(gpvalue)+';'+'EnsDiv='+str(gED)+';'+'fMFE='+str(gfMFE)+';'+'Sequence='+gsequence+';'+'MFE_Fold='+gfold+';'+'Centroid='+gcentroid+'\n')

                    else:
                        print("passed length test")
                        icoordinate = (int(gdata[0])-1)+int(genomic_start)
                        jcoordinate = (int(gdata[1]))+int(genomic_start)
                        window_size = jcoordinate - icoordinate
                    #        if strand == '-1':
                    #            icoordinate = int(int(genomic_start)+(int(length)-int(data[1])))
                    #            jcoordinate = int(int(genomic_start)+(int(length)-int(data[0])))
                        gMFE = float(gdata[2])
                        grand_MFE = float(gdata[3])
                        gzscore = gdata[4]
                        if gzscore == "Undef":
                            gzscore = float(00000)
                        else:
                            gzscore = float(gdata[4])
                        gpvalue = gdata[5]
                        try:
                            gpvalue = float(gdata[5])
                        except ValueError:
                            print(str(gdata[5]))
                            gpvalue =float(0)
                        gED = float(gdata[6])
                        gfMFE = float(gdata[7])
                        gsequence = gdata[8]
                        gfold = gdata[9]
                        gcentroid = gdata[10]
                        gfA = gdata[11]
                        gfG = gdata[12]
                        gfC = gdata[13]
                        gfU = gdata[14]

                        #corrected_file.write(gdata[0]+'\t'+gdata[1]+'\t'+gdata[2]+'\t'+gdata[3]+'\t'+gdata[4]+'\t'+gdata[5]+'\t'+gdata[6]+'\t'+gdata[7]+'\t'+gdata[8]+'\t'+gdata[9]+'\t'+gdata[10]+'\t'+gdata[11]+'\t'+gdata[12]+'\t'+gdata[13]+'\t'+gdata[14])
                        gff3file.write(chromosome+'\t'+'.'+'\t'+'sequence_attribute'+'\t'+str(icoordinate)+'\t'+str(jcoordinate)+'\t'+'.'+'\t'+strand+'\t'+'.\t'+'MFE='+str(gMFE)+';'+'Z-score='+str(gzscore)+';'+'P-value='+str(gpvalue)+';'+'EnsDiv='+str(gED)+';'+'fMFE='+str(gfMFE)+';'+'Sequence='+gsequence+';'+'MFE_Fold='+gfold+';'+'Centroid='+gcentroid+'\n')

                else:
                    print("else")
                    if len(gdata) > 14:
                        print("Errors found in file")
                        print("Error in column six:", gdata)
                        gdata.remove(gdata[5])
                        #corrected_row = ('\t'.join(data))
                        print("Error removed:", gdata)
                        #print(row)
                        #icoordinate = int((int(data[0])+int(genomic_start)-1))
                        #jcoordinate = int(int(data[1])+int(genomic_start))
                        icoordinate = (int(gdata[0])-1)+int(genomic_start)
                        jcoordinate = (int(gdata[1]))+int(genomic_start)
                        window_size = jcoordinate - icoordinate
                #        if strand == '-1':
                #            icoordinate = int(int(genomic_start)+(int(length)-int(data[1])))
                #            jcoordinate = int(int(genomic_start)+(int(length)-int(data[0])))
                        MFE = float(gdata[2])
                        #rand_MFE = float(gdata[3])
                        zscore = gdata[3]
                        if zscore == "Undef":
                            zscore = float(00000)
                        else:
                            zscore = float(gdata[3])
                        pvalue = gdata[4]
                        #try:
                        #    pvalue = float(gdata[5])
                        #except ValueError:
                        #    print(str(gdata[5]))
                        #    pvalue =float(0)
                        gED = float(gdata[5])
                        gfMFE = float(gdata[6])
                        gwindow_sequence = gdata[7]
                        gfold = gdata[8]
                        gcentroid = gdata[9]
                        gfA = gdata[10]
                        gfG = gdata[11]
                        gfC = gdata[12]
                        gfU = gdata[13]
                        #corrected_file.write(gdata[0]+'\t'+gdata[1]+'\t'+gdata[2]+'\t'+gdata[3]+'\t'+gdata[4]+'\t'+gdata[5]+'\t'+gdata[6]+'\t'+gdata[7]+'\t'+gdata[8]+'\t'+gdata[9]+'\t'+gdata[10]+'\t'+gdata[11]+'\t'+gdata[12]+'\t'+gdata[13]+'\t'+gdata[14])
                        gff3file.write(chromosome+'\t'+'.'+'\t'+'sequence_attribute'+'\t'+str(icoordinate)+'\t'+str(jcoordinate)+'\t'+'.'+'\t'+strand+'\t'+'.\t'+'MFE='+str(gMFE)+';'+'Z-score='+str(gzscore)+';'+'P-value='+str(gpvalue)+';'+'EnsDiv='+str(gED)+';'+'fMFE='+str(gfMFE)+';'+'Sequence='+gsequence+';'+'MFE_Fold='+gfold+';'+'Centroid='+gcentroid+'\n')

                    else:
                        #print(len(gdata))
                        icoordinate = (int(gdata[0])-1)+int(genomic_start)
                        jcoordinate = (int(gdata[1]))+int(genomic_start)
                        window_size = jcoordinate - icoordinate
                #        if strand == '-1':
                #            icoordinate = int(int(genomic_start)+(int(length)-int(data[1])))
                #            jcoordinate = int(int(genomic_start)+(int(length)-int(data[0])))
                        gMFE = float(gdata[2])
                        #grand_MFE = float(gdata[3])
                        gzscore = gdata[3]
                        if gzscore == "Undef":
                            gzscore = float(00000)
                        else:
                            gzscore = float(gdata[3])
                        gpvalue = gdata[3]
                        try:
                            gpvalue = float(gdata[4])
                        except ValueError:
                            print(str(gdata[4]))
                            gpvalue =float(0)
                        gED = float(gdata[5])
                        gfMFE = float(gdata[6])
                        gsequence = gdata[7]
                        gfold = gdata[8]
                        gcentroid = gdata[9]
                        gfA = gdata[10]
                        gfG = gdata[11]
                        gfC = gdata[12]
                        gfU = gdata[13]

                        #corrected_file.write(gdata[0]+'\t'+gdata[1]+'\t'+gdata[2]+'\t'+gdata[3]+'\t'+gdata[4]+'\t'+gdata[5]+'\t'+gdata[6]+'\t'+gdata[7]+'\t'+gdata[8]+'\t'+gdata[9]+'\t'+gdata[10]+'\t'+gdata[11]+'\t'+gdata[12]+'\t'+gdata[13]+'\t'+gdata[14])
                        gff3file.write(chromosome+'\t'+'.'+'\t'+'sequence_attribute'+'\t'+str(icoordinate)+'\t'+str(jcoordinate)+'\t'+'.'+'\t'+strand+'\t'+'.\t'+'MFE='+str(gMFE)+';'+'Z-score='+str(gzscore)+';'+'P-value='+str(gpvalue)+';'+'EnsDiv='+str(gED)+';'+'fMFE='+str(gfMFE)+';'+'Sequence='+gsequence+';'+'MFE_Fold='+gfold+';'+'Centroid='+gcentroid+'\n')

    if strand == "-": #Generating GFF3 file for reverse strand
        glines = g.readlines()[1:]
        for row in glines:
            if not row.strip():
                continue
            else:
                gdata = row.split('\t') # this splits each row based on "tab"
                if ((len(str(gdata[8])) == window_size) and
                (float(gdata[3]) != 0) and
                (("A" or "G" or "C" or "U") in gdata[8])) or (len(str(gdata[9])) == window_size and
                (float(gdata[4]) == 0) and
                (("A" or "G" or "C" or "U") in gdata[9])):
                    #float(gdata[7])
                    if len(gdata) > 15:
                        print("Errors found in file")
                        print("Error in column six:", gdata)
                        gdata.remove(gdata[5])
                        #corrected_row = ('\t'.join(data))
                        print("Error removed:", gdata)
                        #print(row)
                        #icoordinate = int((int(data[0])+int(genomic_start)-1))
                        #jcoordinate = int(int(data[1])+int(genomic_start))
                        #print(data)
                        #print(gdata)
                        #print(data)
                        icoordinate = int(int(genomic_start)+(int(length)-int(gdata[1])))
                        jcoordinate = int(int(genomic_start)+(int(length)-int(gdata[0])))
                        gMFE = gdata[2]
                        g_rand_MFE = gdata[3]
                        gzscore = gdata[4]
                        if gzscore == "Undef":
                            gzscore = 00000
                        gpvalue = gdata[5]
                        gED = gdata[6]
                        gfMFE = gdata[7]
                        gsequence = gdata[8]
                        gfold = gdata[9]
                        gcentroid = gdata[10]
                        fA = gdata[11]
                        fG = gdata[12]
                        fC = gdata[13]
                        fU = gdata[14]

                        #corrected_file.write(gdata[0]+'\t'+gdata[1]+'\t'+gdata[2]+'\t'+gdata[3]+'\t'+gdata[4]+'\t'+gdata[5]+'\t'+gdata[6]+'\t'+gdata[7]+'\t'+gdata[8]+'\t'+gdata[9]+'\t'+gdata[10]+'\t'+gdata[11]+'\t'+gdata[12]+'\t'+gdata[13]+'\t'+gdata[14])
                        gff3file.write(chromosome+'\t'+'.'+'\t'+'sequence_attribute'+'\t'+str(icoordinate)+'\t'+str(jcoordinate)+'\t'+'.'+'\t'+strand+'\t'+'.\t'+'MFE='+str(gMFE)+';'+'Z-score='+str(gzscore)+';'+'P-value='+str(gpvalue)+';'+'EnsDiv='+str(gED)+';'+'fMFE='+str(gfMFE)+';'+'Sequence='+gsequence+';'+'MFE_Fold='+gfold+';'+'Centroid='+gcentroid+'\n')

                    else:
                        icoordinate = int(int(genomic_start)+(int(length)-int(gdata[1])))
                        jcoordinate = int(int(genomic_start)+(int(length)-int(gdata[0])))
                        gMFE = gdata[2]
                        g_rand_MFE = gdata[3]
                        gzscore = gdata[4]
                        if gzscore == "Undef":
                            gzscore = 00000
                        gpvalue = gdata[5]
                        gED = gdata[6]
                        gfMFE = gdata[7]
                        gsequence = gdata[8]
                        gfold = gdata[9]
                        gcentroid = gdata[10]
                        fA = gdata[11]
                        fG = gdata[12]
                        fC = gdata[13]
                        fU = gdata[14]

                        #corrected_file.write(gdata[0]+'\t'+gdata[1]+'\t'+gdata[2]+'\t'+gdata[3]+'\t'+gdata[4]+'\t'+gdata[5]+'\t'+gdata[6]+'\t'+gdata[7]+'\t'+gdata[8]+'\t'+gdata[9]+'\t'+gdata[10]+'\t'+gdata[11]+'\t'+gdata[12]+'\t'+gdata[13]+'\t'+gdata[14])
                        gff3file.write(chromosome+'\t'+'.'+'\t'+'sequence_attribute'+'\t'+str(icoordinate)+'\t'+str(jcoordinate)+'\t'+'.'+'\t'+strand+'\t'+'.\t'+'MFE='+str(gMFE)+';'+'Z-score='+str(gzscore)+';'+'P-value='+str(gpvalue)+';'+'EnsDiv='+str(gED)+';'+'fMFE='+str(gfMFE)+';'+'Sequence='+gsequence+';'+'MFE_Fold='+gfold+';'+'Centroid='+gcentroid+'\n')

                else:
                    if len(gdata) > 14:
                        print("Errors found in file")
                        print("Error in column six:", gdata)
                        gdata.remove(gdata[5])
                        #corrected_row = ('\t'.join(data))
                        print("Error removed:", gdata)
                        #print(row)
                        #icoordinate = int((int(data[0])+int(genomic_start)-1))
                        #jcoordinate = int(int(data[1])+int(genomic_start))
                        #print(data)
                        #print(gdata)
                        #print(data)
                        icoordinate = int(int(genomic_start)+(int(length)-int(gdata[1])))
                        jcoordinate = int(int(genomic_start)+(int(length)-int(gdata[0])))
                        gMFE = gdata[2]
                        #g_rand_MFE = gdata[3]
                        gzscore = gdata[3]
                        if gzscore == "Undef":
                            gzscore = 00000
                        gpvalue = gdata[4]
                        gED = gdata[5]
                        gfMFE = gdata[7]
                        gsequence = gdata[8]
                        gfold = gdata[9]
                        gcentroid = gdata[10]
                        fA = gdata[11]
                        fG = gdata[12]
                        fC = gdata[13]
                        fU = gdata[14]

                        #corrected_file.write(gdata[0]+'\t'+gdata[1]+'\t'+gdata[2]+'\t'+gdata[3]+'\t'+gdata[4]+'\t'+gdata[5]+'\t'+gdata[6]+'\t'+gdata[7]+'\t'+gdata[8]+'\t'+gdata[9]+'\t'+gdata[10]+'\t'+gdata[11]+'\t'+gdata[12]+'\t'+gdata[13]+'\t'+gdata[14])
                        gff3file.write(chromosome+'\t'+'.'+'\t'+'sequence_attribute'+'\t'+str(icoordinate)+'\t'+str(jcoordinate)+'\t'+'.'+'\t'+strand+'\t'+'.\t'+'MFE='+str(gMFE)+';'+'Z-score='+str(gzscore)+';'+'P-value='+str(gpvalue)+';'+'EnsDiv='+str(gED)+';'+'fMFE='+str(gfMFE)+';'+'Sequence='+gsequence+';'+'MFE_Fold='+gfold+';'+'Centroid='+gcentroid+'\n')

                    else:
                        icoordinate = int(int(genomic_start)+(int(length)-int(gdata[1])))
                        jcoordinate = int(int(genomic_start)+(int(length)-int(gdata[0])))
                        gMFE = gdata[2]
                        #g_rand_MFE = gdata[3]
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

                        #corrected_file.write(gdata[0]+'\t'+gdata[1]+'\t'+gdata[2]+'\t'+gdata[3]+'\t'+gdata[4]+'\t'+gdata[5]+'\t'+gdata[6]+'\t'+gdata[7]+'\t'+gdata[8]+'\t'+gdata[9]+'\t'+gdata[10]+'\t'+gdata[11]+'\t'+gdata[12]+'\t'+gdata[13]+'\t'+gdata[14])
                        gff3file.write(chromosome+'\t'+'.'+'\t'+'sequence_attribute'+'\t'+str(icoordinate)+'\t'+str(jcoordinate)+'\t'+'.'+'\t'+strand+'\t'+'.\t'+'MFE='+str(gMFE)+';'+'Z-score='+str(gzscore)+';'+'P-value='+str(gpvalue)+';'+'EnsDiv='+str(gED)+';'+'fMFE='+str(gfMFE)+';'+'Sequence='+gsequence+';'+'MFE_Fold='+gfold+';'+'Centroid='+gcentroid+'\n')

with open(filename, 'r') as f:
    if strand == "+": #Generating BW tracks for forward strand.
        genomic_start = int(genomic_start)
        lines = f.readlines()[1:]
        for row in lines:
            if not row.strip():
                continue
            else:
                data = row.split('\t') # this splits each row based on "tab"
                if ((len(str(data[8])) == window_size) and
                (float(data[3]) != 0) and
                (("A" or "G" or "C" or "U") in data[8])) or (len(str(data[9])) == window_size and
                (float(data[4]) == 0) and
                (("A" or "G" or "C" or "U") in data[9])):
                    if len(data) > 15:
                        print("Errors found in file")
                        print("Error in column six:", data)
                        data.remove(data[5])
                        #corrected_row = ('\t'.join(data))
                        print("Error removed:", data)

                        i = data[0]
                        j = data[1]
                        genomic_end = int(genomic_start)+int(step_size)
                        MFE = float(data[2])
                        MFE_list.append(MFE)
                        if data[4] == "Undef":
                            zscore = float(00000)
                        else:
                            zscore = float(data[4])
                        zscore_list.append(zscore)
                        pscore = float(data[5])
                        pscore_list.append(pscore)
                        ED = float(data[6])
                        ED_list.append(ED)
                        fMFE = float(data[7])
                        fMFE_list.append(fMFE)

                    #print(len(data))
                    else:
                        i = data[0]
                        j = data[1]
                        genomic_end = int(genomic_start)+int(step_size)
                        MFE = float(data[2])
                        MFE_list.append(MFE)
                        if data[4] == "Undef":
                            zscore = float(00000)
                        else:
                            zscore = float(data[4])
                        zscore_list.append(zscore)
                        pscore = float(data[5])
                        pscore_list.append(pscore)
                        ED = float(data[6])
                        ED_list.append(ED)
                        fMFE = float(data[7])
                        fMFE_list.append(fMFE)

                else:
                    if len(data) > 14:
                        print("Errors found in file")
                        print("Error in column six:", data)
                        data.remove(data[5])
                        #corrected_row = ('\t'.join(data))
                        print("Error removed:", data)

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

                    #print(len(data))
                    else:
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

    if strand == "-":
        lines = reversed(open(filename).readlines()[1:])
        start = genomic_start
        genomic_start = int(start) + int(window_size)
        for row in lines:
            if not row.strip():
                continue
            else:
                data = row.split('\t') # this splits each row based on "tab"
                if ((len(str(data[8])) == window_size) and
                (float(data[3]) != 0) and
                (("A" or "G" or "C" or "U") in data[8])) or (len(str(data[9])) == window_size and
                (float(data[4]) == 0) and
                (("A" or "G" or "C" or "U") in data[9])):
                    if len(data) > 15:
                        print("Errors found in file")
                        print("Error in column six:", data)
                        data.remove(data[5])
                        #corrected_row = ('\t'.join(data))
                        print("Error removed:", data)


                        #print(row)
                        #print(len(data))
                        #print(row)
                        i = data[0]
                        j = data[1]
                        #genomic_start = int(genomic_start)+int(window_size)
                        MFE = float(data[2])
                        MFE_list.append(MFE)
                        if data[4] == "Undef":
                            zscore = float(00000)
                        else:
                            zscore = float(data[4])
                        zscore_list.append(zscore)
                        pscore = float(data[5])
                        pscore_list.append(pscore)
                        ED = float(data[6])
                        ED_list.append(ED)
                        fMFE = float(data[7])
                        fMFE_list.append(fMFE)

                    else:
                        i = data[0]
                        j = data[1]
                        #genomic_start = int(genomic_start)+int(window_size)
                        MFE = float(data[2])
                        MFE_list.append(MFE)
                        if data[4] == "Undef":
                                zscore = float(00000)
                        else:
                            zscore = float(data[4])
                        zscore_list.append(zscore)
                        pscore = float(data[5])
                        pscore_list.append(pscore)
                        ED = float(data[6])
                        ED_list.append(ED)
                        fMFE = float(data[7])
                        fMFE_list.append(fMFE)

                else:
                    if len(data) > 14:
                        print("Errors found in file")
                        print("Error in column six:", data)
                        data.remove(data[5])
                        #corrected_row = ('\t'.join(data))
                        print("Error removed:", data)

                        #print(row)
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
                    else:
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
