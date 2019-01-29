#!/bin/bash

#SBATCH -J dada2_cutadapt
#SBATCH -o dada2_cutadapt_out
#SBATCH -e dada2_cutadapt_error
#SBATCH -n 8
#SBATCH --mem=100G
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=ALL


module load gcc/5.2.0
module load  R/3.4.3 
module load dada2

module load cutadapt

#unzip 

gunzip *.fastq.gz

#cut adapt to remove forward, reverse and their reverse complements. 

for i in *_R1_001.fastq;

do
  SAMPLE=$(echo ${i} | sed "s/_R1_\001\.fastq//") 
  echo ${SAMPLE}_R1_001.fastq ${SAMPLE}_R2_001.fastq
  cutadapt -a regular 3 end primer -g regular 5 end primer -A reverse comp -a -G reverse comp of -g -o  ${SAMPLE}_trimmed_R1_001.fastq -p ${SAMPLE}_trimmed_R2_001.fastq  ${SAMPLE}_R1_001.fastq  ${SAMPLE}_R2_001.fastq 

done

#run dada2 script on trimmed sequences

Rscript --vanilla --max-ppsize=5000000 dada2_script.R