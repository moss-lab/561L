#!/usr/local/bin/python3.6
"""
Use to prep BT motif files and create the slurm script for SimRNA.
Will calculate for ALL ".bt" files in current directory.

To run:
python SimRNA_inputs-561L.py

-Snake
"""

import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', type = str, default = 'None', help = 'input filename') #if left blank, program will scan all motif files in current directory meeting ".bt" format
parser.add_argument('-o', type = str, default = 'rna', help = 'output filename') #no headers
parser.add_argument('-s', type = str, default = 'None', help = 'location of SimRNA directory (NOT the program)')
parser.add_argument('-c', type = str, default = 'None', help = 'location of config file') #iterations defined in config.dat file

#SimRNA settings
parser.add_argument('-M', type = int, default = 8, help = 'number of model instances to run')
parser.add_argument('-R', type = int, default = 10, help = 'number of replica exchanges')

#SLURM Settings
parser.add_argument('-P', type = str, default = 'whatever', help = 'partition')
parser.add_argument('-n', type = str, default = 'biocrunch9', help = 'server node')
parser.add_argument('-t', type = int, default = 3, help = 'number of days to run') #Increase if running multiple motifs
parser.add_argument('-N', type = int, default = 1, help = 'number of nodes')
parser.add_argument('-T', type = int, default = 50, help = 'number of tasks per node')
parser.add_argument('-J', type = str, default = 'SimRNA', help = 'job name')
parser.add_argument('-e', type = str, default = 'None', help = 'user email')

args = parser.parse_args()
dbn = args.i
out = args.o
simrna = args.s
config = args.c
replicas = args.R
instances = args.M
replicas = args.R
nodes = args.N
tasks = args.T
partition = args.P
server = args.n
runtime = args.t
jobname = args.J
email = args.e

current_directory = os.getcwd()
if simrna == 'None':
        default_SimRNA_directory = '/work/LAS/wmoss-lab/scripts/SimRNA'
        SimRNA_program = default_SimRNA_directory+'/SimRNA'
        SimRNA_data = default_SimRNA_directory+'/data'
        SimRNA_clustering = default_SimRNA_directory+'/clustering'
        SimRNA_trafl2pdb = default_SimRNA_directory+'/SimRNA_trafl2pdbs'
        if config == 'None':
                SimRNA_config = default_SimRNA_directory+'/config.dat'
        else:
                SimRNA_config = config
else:
        SimRNA_program = simrna+'/SimRNA'
        SimRNA_data = simrna+'/data'
        SimRNA_clustering = simrna+'/clustering'
        SimRNA_trafl2pdb = simrna+'/SimRNA_trafl2pdbs'
        if config == 'None':
                SimRNA_config = simrna+'/config.dat'
        else:
                SimRNA_config = config
path = current_directory+'/data'
if os.path.exists(path) != True:
      os.symlink(SimRNA_data, path)

with open('runsecond.sh','w') as slurm:
        slurm.writelines('#!/bin/bash -l'+'\n')
        slurm.writelines('#SBATCH --partition='+partition+'\n')
        slurm.writelines('#SBATCH --nodelist='+server+'\n')
        slurm.writelines('#SBATCH --time='+str(runtime)+'-00:00:00'+'\n')
        slurm.writelines('#SBATCH --nodes='+str(nodes)+'\n')
        slurm.writelines('#SBATCH --ntasks-per-node='+str(tasks)+'\n')
        slurm.writelines('#SBATCH --job-name='+jobname+'\n')
        slurm.writelines('#SBATCH --mail-user='+email+'\n')
        slurm.writelines('#SBATCH --mail-type=ALL'+'\n'+'\n')
        if dbn == 'None':
                #print('No DBN provided, scanning for all .bt files in current directory')
                all_files = os.listdir(current_directory)
                for filename in all_files:
                        if filename.endswith('.bt'):
                                headers_extensions = os.path.splitext(filename)
                                for dbn_headers in headers_extensions:
                                        if dbn_headers != '.bt':
                                                #dbn_header_part = dbn_headers.split('_')
                                                #if dbn_header_part[1] == 'motif':
                                                sequence = dbn_headers+'.sequence'
                                                structure = dbn_headers+'.structure'
                                                with open(filename,'r') as dbn_in:
                                                        with open(sequence,'w') as sequence_out:
                                                                line = dbn_in.readlines()
                                                                sequence_out.write(line[1])
                                                                sequence_length = len(line[1])
                                                                rmsd = 0.0
                                                                rmsd = float(sequence_length)/10.0
                                                with open(filename,'r') as dbn_in:
                                                        with open(structure,'w') as structure_out:
                                                                line = dbn_in.readlines()
                                                                structure_out.write(line[2])
                                                if instances == 1:
                                                        slurm.writelines(SimRNA_program+' -c '+SimRNA_config+' -E '+str(replicas)+' -s '+sequence+' -S '+structure+' -o '+dbn_headers+' >& '+dbn_headers+'.log &'+'\n')
                                                elif instances < 1:
                                                        print('error in number of instances: less than one')
                                                elif instances > 100:
                                                        print('error in number of instances: too big')
                                                else:
                                                        for instance_count in range(instances):
                                                                if instance_count < 10:
                                                                        slurm.writelines(SimRNA_program+' -c '+SimRNA_config+' -E '+str(replicas)+' -s '+sequence+' -S '+structure+' -o '+dbn_headers+'_0'+str(instance_count)+' >& '+dbn_headers+'_0'+str(instance_count)+'.log &'+'\n')
                                                                else:
                                                                        slurm.writelines(SimRNA_program+' -c '+SimRNA_config+' -E '+str(replicas)+' -s '+sequence+' -S '+structure+' -o '+dbn_headers+'_'+str(instance_count)+' >& '+dbn_headers+'_'+str(instance_count)+'.log &'+'\n')
                                                slurm.writelines('wait;'+'\n')
                                                if instances != 1:
                                                        slurm.writelines('cat '+dbn_headers+'_??')
                                                        if replicas != 1:
                                                                slurm.writelines('_??')
                                                        else:
                                                                slurm.writelines('_119')
                                                        slurm.writelines('.trafl > '+dbn_headers+'_all.trafl &'+'\n')
                                                        slurm.writelines('wait;'+'\n')
                                                        slurm.writelines(SimRNA_clustering+' '+dbn_headers+'_all.trafl 0.01 5.0 7.0 10.0 '+str(float(rmsd))+' >& '+dbn_headers+'_clust.log &'+'\n')
                                                        slurm.writelines('wait;'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_01_01-000001.pdb '+dbn_headers+'_all_thrs5.00A_clust01.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_01_01-000001.pdb '+dbn_headers+'_all_thrs5.00A_clust02.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_01_01-000001.pdb '+dbn_headers+'_all_thrs5.00A_clust03.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_01_01-000001.pdb '+dbn_headers+'_all_thrs7.00A_clust01.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_01_01-000001.pdb '+dbn_headers+'_all_thrs7.00A_clust02.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_01_01-000001.pdb '+dbn_headers+'_all_thrs7.00A_clust03.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_01_01-000001.pdb '+dbn_headers+'_all_thrs10.00A_clust01.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_01_01-000001.pdb '+dbn_headers+'_all_thrs10.00A_clust02.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_01_01-000001.pdb '+dbn_headers+'_all_thrs10.00A_clust03.trafl 1 AA &'+'\n') 
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_01_01-000001.pdb '+dbn_headers+'_all_thrs'+str(float(rmsd))+'0A_clust01.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_01_01-000001.pdb '+dbn_headers+'_all_thrs'+str(float(rmsd))+'0A_clust02.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_01_01-000001.pdb '+dbn_headers+'_all_thrs'+str(float(rmsd))+'0A_clust03.trafl 1 AA &'+'\n')
                                                        slurm.writelines('wait;'+'\n')
                                                else:
                                                        if replicas != 1:
                                                                slurm.writelines('cat '+dbn_headers+'_??.trafl > '+dbn_headers+'_all.trafl &'+'\n')
                                                                slurm.writelines('wait;'+'\n')
                                                        slurm.writelines(SimRNA_clustering+' '+dbn_headers)
                                                        if replicas != 1:
                                                                slurm.writelines('_all')
                                                        else:
                                                                slurm.writelines('_119')
                                                        slurm.writelines('.trafl 0.01 5.0 7.0 10.0 '+str(float(rmsd))+' >& '+dbn_headers+'_clust.log &'+'\n')
                                                        slurm.writelines('wait;'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_119-000001.pdb '+dbn_headers+'_119_thrs5.00A_clust01.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_119-000001.pdb '+dbn_headers+'_119_thrs5.00A_clust02.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_119-000001.pdb '+dbn_headers+'_119_thrs5.00A_clust03.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_119-000001.pdb '+dbn_headers+'_119_thrs7.00A_clust01.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_119-000001.pdb '+dbn_headers+'_119_thrs7.00A_clust02.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_119-000001.pdb '+dbn_headers+'_119_thrs7.00A_clust03.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_119-000001.pdb '+dbn_headers+'_119_thrs10.00A_clust01.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_119-000001.pdb '+dbn_headers+'_119_thrs10.00A_clust02.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_119-000001.pdb '+dbn_headers+'_119_thrs10.00A_clust03.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_119-000001.pdb '+dbn_headers+'_119_thrs'+str(float(rmsd))+'0A_clust01.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_119-000001.pdb '+dbn_headers+'_119_thrs'+str(float(rmsd))+'0A_clust02.trafl 1 AA &'+'\n')
                                                        slurm.writelines(SimRNA_trafl2pdb+' '+dbn_headers+'_119-000001.pdb '+dbn_headers+'_119_thrs'+str(float(rmsd))+'0A_clust03.trafl 1 AA &'+'\n')
                                                        slurm.writelines('wait;'+'\n')
        else:
                single_dbn = dbn
                sequence = out+'.sequence'
                structure = out+'.structure'
                with open(single_dbn,'r') as dbn_in:
                        with open(sequence,'w') as sequence_out:
                                line = dbn_in.readlines()
                                sequence_out.write(line[1])
                                sequence_length = len(line[1])
                                rmsd = 0.0
                                rmsd = float(sequence_length)/10.0
                with open(single_dbn,'r') as dbn_in:
                        with open(structure,'w') as structure_out:
                                line = dbn_in.readlines()
                                structure_out.write(line[2])
                if instances == 1:
                        slurm.writelines(SimRNA_program+' -c '+SimRNA_config+' -E '+str(replicas)+' -s '+sequence+' -S '+structure+' -o '+out+' >& '+out+'.log &'+'\n')
                elif instances < 1:
                        print('error in number of instances: less than one')
                elif instances > 100:
                        print('error in number of instances: too big')
                else:
                        for instance_count in range(instances):
                                if instance_count < 10:
                                        slurm.writelines(SimRNA_program+' -c '+SimRNA_config+' -E '+str(replicas)+' -s '+sequence+' -S '+structure+' -o '+out+'_0'+str(instance_count)+' >& '+out+'_0'+str(instance_count)+'.log &'+'\n')
                                else:
                                        slurm.writelines(SimRNA_program+' -c '+SimRNA_config+' -E '+str(replicas)+' -s '+sequence+' -S '+structure+' -o '+out+'_'+str(instance_count)+' >& '+out+'_'+str(instance_count)+'.log &'+'\n')
                slurm.writelines('wait;'+'\n')
                if instances != 1:
                        slurm.writelines('cat '+out+'_??')
                        if replicas != 1:
                                slurm.writelines('_??')
                        else:
                                slurm.writelines('_119')
                        slurm.writelines('.trafl > '+out+'_all.trafl &'+'\n')
                        slurm.writelines('wait;'+'\n')
                        slurm.writelines(SimRNA_clustering+' '+out+'_all.trafl 0.01 5.0 7.0 10.0 '+str(float(rmsd))+' >& '+out+'_clust.log &'+'\n')
                        slurm.writelines('wait;'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_01_01-000001.pdb '+out+'_all_thrs5.00A_clust01.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_01_01-000001.pdb '+out+'_all_thrs5.00A_clust02.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_01_01-000001.pdb '+out+'_all_thrs5.00A_clust03.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_01_01-000001.pdb '+out+'_all_thrs7.00A_clust01.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_01_01-000001.pdb '+out+'_all_thrs7.00A_clust02.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_01_01-000001.pdb '+out+'_all_thrs7.00A_clust03.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_01_01-000001.pdb '+out+'_all_thrs10.00A_clust01.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_01_01-000001.pdb '+out+'_all_thrs10.00A_clust02.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_01_01-000001.pdb '+out+'_all_thrs10.00A_clust03.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_01_01-000001.pdb '+out+'_all_thrs'+str(float(rmsd))+'0A_clust01.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_01_01-000001.pdb '+out+'_all_thrs'+str(float(rmsd))+'0A_clust02.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_01_01-000001.pdb '+out+'_all_thrs'+str(float(rmsd))+'0A_clust03.trafl 1 AA &'+'\n')
                        slurm.writelines('wait;'+'\n')
                else:
                        if replicas != 1:
                                slurm.writelines('cat '+out+'_??.trafl > '+out+'_all.trafl &'+'\n')
                                slurm.writelines('wait;'+'\n')
                        slurm.writelines(SimRNA_clustering+' '+out)
                        if replicas != 1:
                                slurm.writelines('_all')
                        else:
                                slurm.writelines('_119')
                        slurm.writelines('.trafl 0.01 5.0 7.0 10.0 '+str(float(rmsd))+' >& '+out+'_clust.log &'+'\n')
                        slurm.writelines('wait;'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_119-000001.pdb '+out+'_119_thrs5.00A_clust01.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_119-000001.pdb '+out+'_119_thrs5.00A_clust02.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_119-000001.pdb '+out+'_119_thrs5.00A_clust03.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_119-000001.pdb '+out+'_119_thrs7.00A_clust01.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_119-000001.pdb '+out+'_119_thrs7.00A_clust02.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_119-000001.pdb '+out+'_119_thrs7.00A_clust03.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_119-000001.pdb '+out+'_119_thrs10.00A_clust01.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_119-000001.pdb '+out+'_119_thrs10.00A_clust02.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_119-000001.pdb '+out+'_119_thrs10.00A_clust03.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_119-000001.pdb '+out+'_119_thrs'+str(float(rmsd))+'0A_clust01.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_119-000001.pdb '+out+'_119_thrs'+str(float(rmsd))+'0A_clust02.trafl 1 AA &'+'\n')
                        slurm.writelines(SimRNA_trafl2pdb+' '+out+'_119-000001.pdb '+out+'_119_thrs'+str(float(rmsd))+'0A_clust03.trafl 1 AA &'+'\n')
                        slurm.writelines('wait;'+'\n')
