#!/bin/bash --login
#SBATCH --job-name=lfo_hyp
#SBATCH -p multicore
#SBATCH -n 4
#SBATCH -t 1-12
#SBATCH --mem=16G
#SBATCH --output=lfo_hyp_%j.log

module purge
module load apps/gcc/R/4.5.0
module load libs/gcc/flexiblas/3.4.5
module load libs/gcc/glpk/5.0

cd ~/scratch/srl_project
Rscript run_lfo.R hyp
