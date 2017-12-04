#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=0:6:0
#$ -l mem=1G
#$ -l tmpfs=10G
#$ -t 1-5
#$ -N cpy_mod2
#$ -hold_jid smc_mod2

#$ -wd abcsmc-batch
module load gsl

#execute code
./copypars Job=$SGE_TASK_ID
