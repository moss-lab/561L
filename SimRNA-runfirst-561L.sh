#!/bin/bash -l
#SBATCH --partition=whatever            #Partition to submit to
#SBATCH --nodelist=biocrunch9
#SBATCH --time=0-01:00:00               #Time limit for this job
#SBATCH --nodes=1                       #Nodes to be used for this job during runtime
#SBATCH --ntasks-per-node=1             #Number of tasks to fire
#SBATCH --job-name=SimRNA-prep          #Name of this job in work queue
#SBATCH --mail-user=@iastate.edu  #Email to send notifications to
#SBATCH --mail-type=ALL                 #Email notification type (BEGIN, END, FAIL, ALL)

python /work/LAS/wmoss-lab/scripts/SimRNA_inputs-561L.py -t 7 -e @iastate.edu &
wait;

#Materials:
#SimRNA-runfirst.sh (this file)
#".bt" files (individual motif files, will run all ".bt" files in current directory)

#Methods:
#Add email twice to this script (line 8, after "-e" on line 11)
#command line: sbatch SimRNA-runfirst-561L.sh
#validate files and script (increase time for multiple motifs, or use "-t" option above)
#command line: sbatch runsecond.sh
#profit

#-Snake