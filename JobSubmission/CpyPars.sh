#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=0:6:0
#$ -l mem=1G
#$ -l tmpfs=10G
#$ -t 1-40
#$ -N cpy_mod
#$ -hold_jid smc_mod

#$ -wd abcsmc-batch
module load gsl

#execute code
./copypars Job=$SGE_TASK_ID
