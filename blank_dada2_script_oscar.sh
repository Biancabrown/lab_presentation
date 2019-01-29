#!/bin/bash

#SBATCH -J dada2
#SBATCH -o dada2_out
#SBATCH -e dada2_error
#SBATCH -n 8
#SBATCH --mem=100G
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=ALL


module load gcc/5.2.0
module load  R/3.4.3 
module load dada2


#run dada2 script on trimmed sequences

Rscript --vanilla --max-ppsize=5000000 dada2_script.R