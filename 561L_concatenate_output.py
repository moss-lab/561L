#This program will take a the output from the ZScore_calculation.pl program and
# output the sequences of regions which are below the user input zscore
#
# Usage:
#
# $ python thisScrip.py inputfile fastafile
#

import sys # This will allow for the use of system argument inputs
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np

# Definition for an interval.
class Interval:
    def __init__(self, s=0, e=0):
        self.start = s
        self.end = e

    def __repr__(self):
        return "[{}, {}]".format(self.start, self.end)

# This will create an object based on the merging of the intervals
class Solution(object):
    def merge(self, intervals):
        """
        :type intervals: List[Interval]
        :rtype: List[Interval]
        """
        if not intervals:
            return intervals
        intervals.sort(key=lambda x: x.start)
        result = [intervals[0]]
        for i in range(1, len(intervals)):
            prev, current = result[-1], intervals[i]
            if current.start <= prev.end:
                prev.end = max(prev.end, current.end)
            else:
                result.append(current)
        return result
###

filename = sys.argv[1] # this should be the output of a z-score analysis in tab-delimited format
fastafile = sys.argv[2]
#output = sys.argv[2]

threshold = input('Enter Z-score cutoff value (e.g., -1): ')

hotspot_regions = []
MFE_list = []
z_score_list = []
p_score_list = []
ED_list = []
frags = []


# Open fasta with SeqIO to grab sequence
for cur_record in SeqIO.parse(fastafile, "fasta"):
    gene_sequence = cur_record.seq
    #print(cur_record.name)
    #print(gene_sequence)

#corrected_file = open(filename+'.corrected.txt', 'w')
#corrected_file.write("i\tj\tMFE\trandomMFE\tZscore\tPscore\tED\tfMFE\tSequence\tFold\tCentroid\t#A\t#G\t#C\t#U")

# Open fasta with python to grab the header
with open(fastafile, 'r') as f:
    header = f.readlines()[:1]
    #print(header)
    header_data = re.split('\s', str(header))
    chromosome_data = header_data[2]
    split_chromosome_data = re.split('\:', chromosome_data)
    chromosome = split_chromosome_data[2]
    genomic_start = split_chromosome_data[3]
    genomic_end = split_chromosome_data[4]
    strand_data = split_chromosome_data[5]
    split_strand_data = re.split('\n', strand_data)
    strand = strand_data[0]

length = (int(genomic_end) - int(genomic_start)) + 1
#print(strand)
# Open the z-score output file, and split each row into its metrics
if strand == "1":
    #print("positive_strand")
    with open(filename, 'r') as g:
        lines = g.readlines()[1:]
        for row in lines:
            data = row.split('\t')
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
                    #icoordinate = int((int(data[0])+int(genomic_start)-1))
                    #jcoordinate = int(int(data[1])+int(genomic_start))
                    icoordinate = int(data[0])-1
                    jcoordinate = int(data[1])
                    window_size = jcoordinate - icoordinate
            #        if strand == '-1':
            #            icoordinate = int(int(genomic_start)+(int(length)-int(data[1])))
            #            jcoordinate = int(int(genomic_start)+(int(length)-int(data[0])))
                    MFE = float(data[2])
                    rand_MFE = float(data[3])
                    zscore = data[4]
                    if zscore == "Undef":
                        zscore = float(00000)
                    else:
                        zscore = float(data[4])
                    pvalue = data[5]
                    try:
                        pvalue = float(data[5])
                    except ValueError:
                        print(str(data[5]))
                        pvalue =float(0)
                    ED = float(data[6])
                    fMFE = float(data[7])
                    window_sequence = data[8]
                    fold = data[9]
                    centroid = data[10]
                    fA = data[11]
                    fG = data[12]
                    fC = data[13]
                    fU = data[14]
                    #corrected_file.write(data[0]+'\t'+data[1]+'\t'+data[2]+'\t'+data[3]+'\t'+data[4]+'\t'+data[5]+'\t'+data[6]+'\t'+data[7]+'\t'+data[8]+'\t'+data[9]+'\t'+data[10]+'\t'+data[11]+'\t'+data[12]+'\t'+data[13]+'\t'+data[14])
            # With metrics defined, we can place each window with
                    if float(zscore) < float(threshold):
                        window = Interval(icoordinate, jcoordinate)
                        hotspot_regions.append(window)
                        MFE_list.append(MFE)
                        z_score_list.append(zscore)
                        p_score_list.append(pvalue)
                        ED_list.append(ED)
                else:
                    pvalue = data[5]
                    icoordinate = int(data[0])-1
                    jcoordinate = int(data[1])
                    window_size = jcoordinate - icoordinate
                    MFE = float(data[2])
                    rand_MFE = float(data[3])
                    zscore = data[4]
                    if zscore == "Undef":
                        zscore = float(00000)
                    else:
                        zscore = float(data[4])
                    if data[5] == '':
                        print(HELP)
                    try:
                        pvalue = float(data[5])
                    except ValueError:
                        print(str(data[5]))
                        pvalue = float(0)
                    ED = float(data[6])
                    fMFE = float(data[7])
                    window_sequence = data[8]
                    fold = data[9]
                    centroid = data[10]
                    fA = data[11]
                    fG = data[12]
                    fC = data[13]
                    fU = data[14]
                    if float(zscore) < float(threshold):
                        window = Interval(icoordinate, jcoordinate)
                        hotspot_regions.append(window)
                        MFE_list.append(MFE)
                        z_score_list.append(zscore)
                        p_score_list.append(pvalue)
                        ED_list.append(ED)

            else:
                if len(data) > 14:
                    print("Errors found in file")
                    print("Error in column five:", data)
                    data.remove(data[5])
                    #corrected_row = ('\t'.join(data))
                    print("Error removed:", data)
                    #print(row)
                    #icoordinate = int((int(data[0])+int(genomic_start)-1))
                    #jcoordinate = int(int(data[1])+int(genomic_start))
                    icoordinate = int(data[0])-1
                    jcoordinate = int(data[1])
                    window_size = jcoordinate - icoordinate
            #        if strand == '-1':
            #            icoordinate = int(int(genomic_start)+(int(length)-int(data[1])))
            #            jcoordinate = int(int(genomic_start)+(int(length)-int(data[0])))
                    MFE = float(data[2])
                    #rand_MFE = float(data[3])
                    zscore = data[4]
                    if zscore == "Undef":
                        zscore = float(00000)
                    else:
                        zscore = float(data[4])
                    pvalue = data[5]
                    try:
                        pvalue = float(data[5])
                    except ValueError:
                        print(str(data[4]))
                        pvalue =float(0)
                    ED = float(data[5])
                    fMFE = float(data[6])
                    window_sequence = data[7]
                    fold = data[8]
                    centroid = data[9]
                    fA = data[10]
                    fG = data[11]
                    fC = data[12]
                    fU = data[13]
                    #corrected_file.write(data[0]+'\t'+data[1]+'\t'+data[2]+'\t'+data[3]+'\t'+data[4]+'\t'+data[5]+'\t'+data[6]+'\t'+data[7]+'\t'+data[8]+'\t'+data[9]+'\t'+data[10]+'\t'+data[11]+'\t'+data[12]+'\t'+data[13]+'\t'+data[14])
            # With metrics defined, we can place each window with
                    if float(zscore) < float(threshold):
                        window = Interval(icoordinate, jcoordinate)
                        hotspot_regions.append(window)
                        MFE_list.append(MFE)
                        z_score_list.append(zscore)
                        p_score_list.append(pvalue)
                        ED_list.append(ED)
                else:
                    pvalue = data[5]
                    icoordinate = int(data[0])-1
                    jcoordinate = int(data[1])
                    window_size = jcoordinate - icoordinate
                    MFE = float(data[2])
                    #rand_MFE = float(data[3])
                    zscore = data[3]
                    if zscore == "Undef":
                        zscore = float(00000)
                    else:
                        zscore = float(data[3])
                    if data[4] == '':
                        print(HELP)
                    try:
                        pvalue = float(data[4])
                    except ValueError:
                        print(str(data[4]))
                        pvalue = float(0)
                    ED = float(data[5])
                    fMFE = float(data[6])
                    window_sequence = data[7]
                    fold = data[8]
                    centroid = data[9]
                    fA = data[10]
                    fG = data[11]
                    fC = data[12]
                    fU = data[13]
                    if float(zscore) < float(threshold):
                        window = Interval(icoordinate, jcoordinate)
                        hotspot_regions.append(window)
                        MFE_list.append(MFE)
                        z_score_list.append(zscore)
                        p_score_list.append(pvalue)
                        ED_list.append(ED)    # This will merge all of the overlapping windows which were below threshold
    concatenated_regions = Solution().merge(hotspot_regions)

    # Now we are going to be working directly on the merged regions, opening output filename
    with open(fastafile+".seq_cutoff_("+threshold+").fa", "w") as out_handle:

        #print(type(concatenated_regions))
        for x in concatenated_regions:
            print("Region coordinates:", x)
            #print(type(x))
            i = x.start
            j = x.end
            #i = Interval.start()
            #j = Interval.end()
            frag = cur_record.seq[i:j]
            record = SeqRecord(frag, '%s dna:chromosome chromosome:GRCh38:%s:%i:%i:%s' % (chromosome, chromosome, i+int(genomic_start), j+int(genomic_start), strand), '', '')
            frags.append(record)

    SeqIO.write(frags, fastafile+".seq_cutoff_("+threshold+").fa", "fasta")
    #print(MFE_list)
    print("Mean MFE of concatenated windows = ", np.mean(MFE_list))
    print("Mean z-score of concatenated windows = ", np.mean(z_score_list))
    print("Mean p-value of concatenated windows = ", np.mean(p_score_list))
    print("Mean ED of concatenated windows = ", np.mean(ED_list))

if strand == "-":
    #print("Negative_strand")
    with open(filename, 'r') as g:
        lines = g.readlines()[1:]
        for row in lines:
            data = row.split('\t')
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
                    #icoordinate = int((int(data[0])+int(genomic_start)-1))
                    #jcoordinate = int(int(data[1])+int(genomic_start))
                    icoordinate = int(data[0])-1
                    jcoordinate = int(data[1])
                    window_size = jcoordinate - icoordinate
                    MFE = float(data[2])
                    rand_MFE = float(data[3])
                    zscore = data[4]
                    if zscore == "Undef":
                        zscore = float(00000)
                    else:
                        zscore = float(data[4])
                    pvalue = data[5]
                    try:
                        pvalue = float(data[5])
                    except ValueError:
                        print(str(data[5]))
                        pvalue =float(0)
                    ED = float(data[6])
                    fMFE = float(data[7])
                    window_sequence = data[8]
                    fold = data[9]
                    centroid = data[10]
                    fA = data[11]
                    fG = data[12]
                    fC = data[13]
                    fU = data[14]
                    #corrected_file.write(data[0]+'\t'+data[1]+'\t'+data[2]+'\t'+data[3]+'\t'+data[4]+'\t'+data[5]+'\t'+data[6]+'\t'+data[7]+'\t'+data[8]+'\t'+data[9]+'\t'+data[10]+'\t'+data[11]+'\t'+data[12]+'\t'+data[13]+'\t'+data[14])
            # With metrics defined, we can place each window with
                    if float(zscore) < float(threshold):
                        window = Interval(icoordinate, jcoordinate)
                        hotspot_regions.append(window)
                        MFE_list.append(MFE)
                        z_score_list.append(zscore)
                        p_score_list.append(pvalue)
                        ED_list.append(ED)
                else:
                    pvalue = data[5]
                    icoordinate = int(data[0])-1
                    jcoordinate = int(data[1])
                    window_size = jcoordinate - icoordinate
                    MFE = float(data[2])
                    rand_MFE = float(data[3])
                    zscore = data[4]
                    if zscore == "Undef":
                        zscore = float(00000)
                    else:
                        zscore = float(data[4])
                    if data[5] == '':
                        print(HELP)
                    try:
                        pvalue = float(data[5])
                    except ValueError:
                        print(str(data[5]))
                        pvalue = float(0)

                    ED = float(data[6])
                    fMFE = float(data[7])
                    window_sequence = data[8]
                    fold = data[9]
                    centroid = data[10]
                    fA = data[11]
                    fG = data[12]
                    fC = data[13]
                    fU = data[14]
                    #corrected_file.write(data[0]+'\t'+data[1]+'\t'+data[2]+'\t'+data[3]+'\t'+data[4]+'\t'+data[5]+'\t'+data[6]+'\t'+data[7]+'\t'+data[8]+'\t'+data[9]+'\t'+data[10]+'\t'+data[11]+'\t'+data[12]+'\t'+data[13]+'\t'+data[14])
            # With metrics defined, we can place each window with
                    if float(zscore) < float(threshold):
                        window = Interval(icoordinate, jcoordinate)
                        hotspot_regions.append(window)
                        MFE_list.append(MFE)
                        z_score_list.append(zscore)
                        p_score_list.append(pvalue)
                        ED_list.append(ED)

            else:
                if len(data) > 14:
                    print("Errors found in file")
                    print("Error in column five:", data)
                    data.remove(data[4])
                    #corrected_row = ('\t'.join(data))
                    print("Error removed:", data)
                    #print(row)
                    #icoordinate = int((int(data[0])+int(genomic_start)-1))
                    #jcoordinate = int(int(data[1])+int(genomic_start))
                    icoordinate = int(data[0])-1
                    jcoordinate = int(data[1])
                    window_size = jcoordinate - icoordinate
                    MFE = float(data[2])
                    #rand_MFE = float(data[3])
                    zscore = data[3]
                    if zscore == "Undef":
                        zscore = float(00000)
                    else:
                        zscore = float(data[3])
                    pvalue = data[4]
                    try:
                        pvalue = float(data[4])
                    except ValueError:
                        print(str(data[4]))
                        pvalue =float(0)
                    ED = float(data[5])
                    fMFE = float(data[6])
                    window_sequence = data[7]
                    fold = data[8]
                    centroid = data[9]
                    fA = data[10]
                    fG = data[11]
                    fC = data[12]
                    fU = data[13]
                    #corrected_file.write(data[0]+'\t'+data[1]+'\t'+data[2]+'\t'+data[3]+'\t'+data[4]+'\t'+data[5]+'\t'+data[6]+'\t'+data[7]+'\t'+data[8]+'\t'+data[9]+'\t'+data[10]+'\t'+data[11]+'\t'+data[12]+'\t'+data[13]+'\t'+data[14])
            # With metrics defined, we can place each window with
                    if float(zscore) < float(threshold):
                        window = Interval(icoordinate, jcoordinate)
                        hotspot_regions.append(window)
                        MFE_list.append(MFE)
                        z_score_list.append(zscore)
                        p_score_list.append(pvalue)
                        ED_list.append(ED)
                else:
                    pvalue = data[4]
                    icoordinate = int(data[0])-1
                    jcoordinate = int(data[1])
                    window_size = jcoordinate - icoordinate
                    MFE = float(data[2])
                    #rand_MFE = float(data[3])
                    zscore = data[3]
                    if zscore == "Undef":
                        zscore = float(00000)
                    else:
                        zscore = float(data[3])
                    if data[4] == '':
                        print(HELP)
                    try:
                        pvalue = float(data[4])
                    except ValueError:
                        print(str(data[4]))
                        pvalue = float(0)
                    ED = float(data[5])
                    fMFE = float(data[6])
                    window_sequence = data[7]
                    fold = data[8]
                    centroid = data[9]
                    fA = data[10]
                    fG = data[11]
                    fC = data[12]
                    fU = data[13]
                    #corrected_file.write(data[0]+'\t'+data[1]+'\t'+data[2]+'\t'+data[3]+'\t'+data[4]+'\t'+data[5]+'\t'+data[6]+'\t'+data[7]+'\t'+data[8]+'\t'+data[9]+'\t'+data[10]+'\t'+data[11]+'\t'+data[12]+'\t'+data[13]+'\t'+data[14])
            # With metrics defined, we can place each window with
                    if float(zscore) < float(threshold):
                        window = Interval(icoordinate, jcoordinate)
                        hotspot_regions.append(window)
                        MFE_list.append(MFE)
                        z_score_list.append(zscore)
                        p_score_list.append(pvalue)
                        ED_list.append(ED)

    # This will merge all of the overlapping windows which were below threshold
    concatenated_regions = Solution().merge(hotspot_regions)

    # Now we are going to be working directly on the merged regions, opening output filename
    with open(fastafile+".seq_cutoff_("+threshold+").fa", "w") as out_handle:

        #print(type(concatenated_regions))
        for x in concatenated_regions:
            print("Region coordinates:", x)
            #print(type(x))
            i = x.start
            j = x.end
            #i = Interval.start()
            #j = Interval.end()
            frag = cur_record.seq[i:j]
            window_start = int(int(genomic_start)+(int(length)-int(j)))
            print(window_start)
            window_end = int(int(genomic_start)+(int(length)-int(i)))-1
            print(window_end)
            record = SeqRecord(frag, '%s dna:chromosome chromosome:GRCh38:%s:%s:%s:%s1' % (chromosome, chromosome, window_start, window_end, strand), '', '')
            frags.append(record)

    SeqIO.write(frags, fastafile+".seq_cutoff_("+threshold+").fa", "fasta")
    #print(MFE_list)
    print("Mean MFE of concatenated windows = ", np.mean(MFE_list))
    print("Mean z-score of concatenated windows = ", np.mean(z_score_list))
    print("Mean p-value of concatenated windows = ", np.mean(p_score_list))
    print("Mean ED of concatenated windows = ", np.mean(ED_list))
