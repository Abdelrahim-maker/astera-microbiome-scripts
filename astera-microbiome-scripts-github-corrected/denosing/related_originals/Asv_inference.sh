#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=2:00:00
#SBATCH --account=PAS1212

module load gnu/9.1.0
module load mkl/2019.0.5
module load R/4.0.2

Rscript /fs/ess/PAS1212/Afaf_rich/000_scripts/final_dada2_post_asv_script.R