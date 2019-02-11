#!/bin/bash

#SBATCH --mail-user=ykent@iastate.edu 
#SBATCH --mail-type=end
#SBATCH --partition=compute 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=32
#SBATCH --job-name=MnO

# record start date and host
date
hostname

/home/ykent/GitLab/comrisb/bin/comrisb.py 

# record end date
date

