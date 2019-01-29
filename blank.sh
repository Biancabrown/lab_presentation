#!/bin/bash

#SBATCH -J dada2_yellow_stone
#SBATCH -o dada2_yellow_stone_out
#SBATCH -e dada2_yellow_stone_error
#SBATCH -n 8
#SBATCH --mem=100G
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=ALL




module load program

code