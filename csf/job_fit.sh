#!/bin/bash --login
#SBATCH --job-name=srl_fit
#SBATCH -p multicore
#SBATCH -n 4
#SBATCH -t 1-0
#SBATCH --mem=16G
#SBATCH --output=srl_fit_%j.log

module purge
module load apps/gcc/R/4.5.0
module load libs/gcc/flexiblas/3.4.5
module load libs/gcc/glpk/5.0

cd ~/scratch/srl_project

Rscript run.R
