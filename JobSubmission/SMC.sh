#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=0:6:0
#$ -l mem=1G
#$ -l tmpfs=10G
#$ -t 1-5000
#$ -N smc_mod

#$ -wd abcsmc-batch
module load gsl

#execute code
./ABCSMC Seed=$SGE_TASK_ID Particles=100 Final=0.1 Job=$SGE_TASK_ID Alpha=0.95 Nparameters=9 Time=0
