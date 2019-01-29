Bianca Palmer Brown Lab Meeting Presentation 

##DADA2

Traditionally, sequence reads are clustered into operational taxonomic units (OTUs) at a defined identity threshold to avoid sequencing errors generating spurious taxonomic units. 

There have been numerous bioinformatic methods recently released that attempt to correct sequencing errors to determine real biological sequences at single nucleotide resolution by generating amplicon sequence variants (ASVs).

*remember to cite*


1. DADA2 is for modeling and correcting Illumina-sequenced amplicon errors 

2. Infers sample sequences exactly and resolves differences of as little as 1 nucleotide.

*remember to cite*


### Target amplicon with primers



![alt-text-1](https://raw.githubusercontent.com/Biancabrown/lab_presentation/master/primer.jpeg)

primers are in green and red
amplicon is in black 
Total length is 250 bp


### Seuqencing 

![alt-text-1](https://raw.githubusercontent.com/Biancabrown/lab_presentation/master/image_2.jpeg "title-1")

### Multiple sequences of different amplicons  
![alt-text-1](https://raw.githubusercontent.com/Biancabrown/lab_presentation/master/image_3.jpeg "title-1")




We get those sequences in a fastq file. We get a fastq file for each of samples 

```
Sample 1

BBBAA?@@DAAFGGFGCAAECGHDCEEFHHHGDCDBFEHHHFCCAAE@BE@DHCC>FFHHHHHHFGHHHGG?EFHHHHHHHGFCEEGHHFHHHHHHHGDDGHGEGFHGHHGFHHHGGHGFFFGDHHGCGC@GGDC@CCGEFA/:G=:CCD:CC0FGCG:AEFFFFGEFFGGGEFGGGFGGGFBFFFD@@@AEFA:9@?FFFFFFFFFDD?DFFFF;DFBBFADDDDB9FFFEA9DDFD:9:=CFFAE/0A
@HWI-M02808:154:ANHE0:1:1101:17870:2090 1:N:0:AGGCATCTTACG
TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGTTGCTAAGTCAGTGGTAAAAGCGTGTGGCTCAACCATACCAAGCCATTGAAACTGGCGATCTTGAGTGTAAACGAGGTAGGCGGAATGTGACGTGTAGCGGTGAAATGCTTAGATATGTCACAGAACCCCGATTGCGAAGGCAGCTTACCAGCATACAACTGACGCTGAGGCACGAAAGCGTGGGTATCAAACA
+
BBBBBBB@AFBFFGGGFAEFAFGEFEEGHHHHFFGBFFHHHFFCBEEEBEEDHGFE>FFGHGHHFGHHHHGFFHHGHHE?EEFFEFGFHHHFFGGHCGGFFHHHHGGHHF?GGGADGHHHFC=>>GDFHDCDC.<DEHGGCDEG;CGHDCE:/9;:B9A.9;B;;0CB90;0;C;CBFFB/CCCDD-;;DE/9---A9ADAFBFFFFF.FBFFFFFFB/F.;AD.?EAEF####################
@HWI-M02808:154:ANHE0:1:1101:17877:2107 1:N:0:AGGCATCTTACG
TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGTTGCTAAGTCAGTGGTAAAAGCGTGTGGCTCAACCATACCAAGCCATTGAAACTGGCGATCTTGAGTGTAAACGAGGTAGGCGGAATGTGACGTGTAGCGGTGAAATGCTTAGATATGTCACAGAACCCCGATTGCGAAGGCAGCTTACCAGCATACAACTGACGCTGAGGCACGAAAGCGTGGGGATCAAACA

```

```
Sample 2

BBBAA?@@DAAFGGFGCAAECGHDCEEFHHHGDCDBFEHHHFCCAAE@BE@DHCC>FFHHHHHHFGHHHGG?EFHHHHHHHGFCEEGHHFHHHHHHHGDDGHGEGFHGHHGFHHHGGHGFFFGDHHGCGC@GGDC@CCGEFA/:G=:CCD:CC0FGCG:AEFFFFGEFFGGGEFGGGFGGGFBFFFD@@@AEFA:9@?FFFFFFFFFDD?DFFFF;DFBBFADDDDB9FFFEA9DDFD:9:=CFFAE/0A
@HWI-M02808:154:ANHE0:1:1101:17870:2090 1:N:0:AGGCATCTTACG
TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGTTGCTAAGTCAGTGGTAAAAGCGTGTGGCTCAACCATACCAAGCCATTGAAACTGGCGATCTTGAGTGTAAACGAGGTAGGCGGAATGTGACGTGTAGCGGTGAAATGCTTAGATATGTCACAGAACCCCGATTGCGAAGGCAGCTTACCAGCATACAACTGACGCTGAGGCACGAAAGCGTGGGTATCAAACA
+
BBBBBBB@AFBFFGGGFAEFAFGEFEEGHHHHFFGBFFHHHFFCBEEEBEEDHGFE>FFGHGHHFGHHHHGFFHHGHHE?EEFFEFGFHHHFFGGHCGGFFHHHHGGHHF?GGGADGHHHFC=>>GDFHDCDC.<DEHGGCDEG;CGHDCE:/9;:B9A.9;B;;0CB90;0;C;CBFFB/CCCDD-;;DE/9---A9ADAFBFFFFF.FBFFFFFFB/F.;AD.?EAEF####################
@HWI-M02808:154:ANHE0:1:1101:17877:2107 1:N:0:AGGCATCTTACG
TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGTTGCTAAGTCAGTGGTAAAAGCGTGTGGCTCAACCATACCAAGCCATTGAAACTGGCGATCTTGAGTGTAAACGAGGTAGGCGGAATGTGACGTGTAGCGGTGAAATGCTTAGATATGTCACAGAACCCCGATTGCGAAGGCAGCTTACCAGCATACAACTGACGCTGAGGCACGAAAGCGTGGGGATCAAACA

```

WE use dada2 to do the following:

1. Check sequence quality
2. Trim primers (can also be used with cut adapt
3. Remove sequences that are not a particular length
4. Check for errors
5. Dereplicate
6. Remove chimeras 
7. Assign taxonomy

A similar class of methods developed for 454-scale data was typically used to ‘denoise’ sequencing data prior to constructing OTUs (Quince et al., 2011), while new ASV methods are explicitly intended to replace OTUs as the atomic unit of analysis.


###General code 
R based

```
library(dada2); packageVersion("dada2")
library(Rcpp); packageVersion("Rcpp")
library(ShortRead); packageVersion("ShortRead")

path <- "fastq files"


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
#filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
#filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,trimLeft = c(19, 20), truncLen=c(231,230),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
#saveRDS(out, "out.rds")

#errF <- learnErrors(filtFs, multithread=TRUE)
#saveRDS(errF, "errF.rds")

#errR <- learnErrors(filtRs, multithread=TRUE)
#saveRDS(errR, "errR.rds")

#derepFs <- derepFastq(filtFs, verbose=TRUE)
#saveRDS(derepFs, "derepFs.rds")

#derepRs <- derepFastq(filtRs, verbose=TRUE)
#saveRDS(derepRs, "derepRs.rds")

#Name the derep-class objects by the sample names

#names(derepFs) <- sample.names

#names(derepRs) <- sample.names

#dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
#saveRDS(dadaFs, "dadaFs.rds")

#dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#saveRDS(dadaRs, "dadaRs.rds")

#mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
#saveRDS(mergers, "mergers.rds")

#Inspect the merger data.frame from the first sample

#seqtab <- makeSequenceTable(mergers)
#saveRDS(seqtab, "seqtab.rds")

#seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#saveRDS(seqtab.nochim, "seqtab.nochim.rds")

#taxa <- assignTaxonomy(seqtab.nochim, "/users/bbrown3/data/bbrown3/2018_22_03_SMH_microbiome/gg_13_8_train_set_97.fa.gz", multithread=TRUE)
#saveRDS(taxa, "taxa.rds")


```


Where our code differs and why?  

```

Bianca microbiome

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,trimLeft = c(19, 20), truncLen=c(231,230),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
saveRDS(out, "out.rds")

Courtney Plant

filterAndTrim(file.path(path,fns), file.path(filtpath,fns),maxN=0, maxEE=1, truncQ=2, rm.phix=TRUE,compress=TRUE, verbose=TRUE, multithread=TRUE)

Patrick 

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(2,2), truncQ=2, minLen=8, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
saveRDS(out, "out.rds")

Brian CO1

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
    truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)


```



Main reason:

1. Variation in size of our amplicon regions
2. How we choose to remove our primers


###Above example



![alt-text-1](https://raw.githubusercontent.com/Biancabrown/lab_presentation/master/primer.jpeg)




```
trim left =c(20, 20)

Remove primers from the beginning of the left forward (20) and reverse (20) sequences. 

TruncLen=c(210,210) 
remove sequences that are shorter than those length 

```


 
![alt-text-1](https://raw.githubusercontent.com/Biancabrown/lab_presentation/master/image_4.jpeg)


###Example Number 2

Image if we had a 300 bp sequences

![alt-text-1](https://raw.githubusercontent.com/Biancabrown/lab_presentation/master/image_5.jpeg)


We are still doing Miseq 250 paired end

```
trim left =c(25,25)

#Trim forward primers and reverse primer. However, it will assume that each sequence has a forward and reverse primer


```



![alt-text-1](https://raw.githubusercontent.com/Biancabrown/lab_presentation/master/image_6.jpeg)



So we will get odd overlaps etc. 

Instead using dada2 to remove primers, some people are usng cut adapt. 


##Cut Adapt in R

```{r}

path <- "files"
list.files(path)

fnFs <- sort(list.files(path, pattern="R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2_001.fastq", full.names = TRUE))

FWD <- "GGWACWGGWTGAACWGTWTAYCCYCC"
REV <- "TANACYTCNGGRTGNCCRAARAAYCA"

allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
        RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))
```



###Cut adapt Oscar 

```

#!/bin/bash

#SBATCH -J cutadapt
#SBATCH -o cutadapt_out
#SBATCH -e cutadapt_out
#SBATCH -n 8
#SBATCH --mem=100G
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=ALL


module unload anaconda/3-4.4.0
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
  cutadapt -a GGGCAATCCTGAGCCAA -g TTGGCTCAGGATTGCCC -A CCATTGAGTCTCTGCACCTATC -G GATAGGTGCAGAGACTCAATGG -o  ${SAMPLE}_trimmed_R1_001.fastq -p ${SAMPLE}_trimmed_R2_001.fastq  ${SAMPLE}_R1_001.fastq  ${SAMPLE}_R2_001.fastq 

done

```


#After running cut adapt, you can run your dada2 code:

If running it on your computer, you can use Brian's script with some modification: 


```{r}

path <- "file"
list.files(path)

fnFs <- sort(list.files(path, pattern="R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2_001.fastq", full.names = TRUE))

FWD <- "GGWACWGGWTGAACWGTWTAYCCYCC"
REV <- "TANACYTCNGGRTGNCCRAARAAYCA"

allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
        RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))
```

```{r setup, include=FALSE}

rm(list = ls())

library(dada2)
library("plyr")
library("ggplot2")
library("phyloseq")
devtools::install_github("GuillemSalazar/FastaUtils")
require(FastaUtils)
library(ShortRead)
library(Biostrings)
```

###DADA2 Pipeline

```{r}
path.cut <- "/Users/MacBookPro/Dropbox (Brown)/Active/SIDE) UF Frog Diet/Demultiplexed/cutadapt"
list.files(path.cut)
cutFs <- sort(list.files(path.cut, pattern="R1_001.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern="R2_001.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(cutFs), "_"), `[`, 1)

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
    truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

head(track)

barplot(track[,6], cex.names=0.5, las=2)

setwd("~/Desktop")
save.image("DADA2_Output_2.R")
```



##If running on oscar then we you need to do the following steps:

1. Write a R script dada2
2. Write a bash script and include the cutadapt code and dada2 code within that R script. 
3. Call call bash script


###Structure of R script

R script are written in R

```
#!/usr/bin/env Rscript


```


###Other things I will like to point out using our dada2 rscript


```
#!/usr/bin/env Rscript

library(dada2); packageVersion("dada2")
library(Rcpp); packageVersion("Rcpp")
library(ShortRead); packageVersion("ShortRead")

path <- "fast q"


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,trimLeft = c(19, 20), truncLen=c(231,230),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
saveRDS(out, "/gpfs_home/bbrown3/data/bbrown3/MSU_READS/dada_20161227_16S_V4_PE_2/out.rds")

errF <- learnErrors(filtFs, multithread=TRUE)
saveRDS(errF, "errF.rds")

errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errR, "/gpfs_home/bbrown3/data/bbrown3/MSU_READS/dada_20161227_16S_V4_PE_2/errR.rds")

derepFs <- derepFastq(filtFs, verbose=TRUE)
saveRDS(derepFs, "/gpfs_home/bbrown3/data/bbrown3/MSU_READS/dada_20161227_16S_V4_PE_2/derepFs.rds")

derepRs <- derepFastq(filtRs, verbose=TRUE)
saveRDS(derepRs, "/gpfs_home/bbrown3/data/bbrown3/MSU_READS/dada_20161227_16S_V4_PE_2/derepRs.rds")

#Name the derep-class objects by the sample names

names(derepFs) <- sample.names

names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
saveRDS(dadaFs, "/gpfs_home/bbrown3/data/bbrown3/MSU_READS/dada_20161227_16S_V4_PE_2/dadaFs.rds")

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
saveRDS(dadaRs, "/gpfs_home/bbrown3/data/bbrown3/MSU_READS/dada_20161227_16S_V4_PE_2/dadaRs.rds")

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
saveRDS(mergers, "/gpfs_home/bbrown3/data/bbrown3/MSU_READS/dada_20161227_16S_V4_PE_2/mergers.rds")

#Inspect the merger data.frame from the first sample

seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "/gpfs_home/bbrown3/data/bbrown3/MSU_READS/dada_20161227_16S_V4_PE_2/seqtab.rds")

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim, "/gpfs_home/bbrown3/data/bbrown3/MSU_READS/dada_20161227_16S_V4_PE_2/seqtab.nochim.rds")

taxa <- assignTaxonomy(seqtab.nochim, "/users/bbrown3/data/bbrown3/2018_22_03_SMH_microbiome/gg_13_8_train_set_97.fa.gz", multithread=TRUE)
saveRDS(taxa, "/gpfs_home/bbrown3/data/bbrown3/MSU_READS/dada_20161227_16S_V4_PE_2/taxa.rds")



```

1. save all your output

saveRDS(errF, "errF.rds")

In case your program crash you can start again from that point. 

TO retrieve your output

readRDS("err.rds")


2. Parts of the script that can be run in parallel: 

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
saveRDS(dadaRs, "/gpfs_home/bbrown3/data/bbrown3/MSU_READS/dada_20161227_16S_V4_PE_2/dadaRs.rds")




###After writing R script, you write a bash script. If we were only going to call the bash script 

General Structure of a bash script


```
#!/bin/bash

#SBATCH -J job_name
#SBATCH -o out
#SBATCH -e error
#SBATCH -n 8
#SBATCH --mem=100G
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=ALL

module load cutadapt

blah blah code

```

Structure of a bash script using cut adapt and dada2 

##Patrick's Example

```

#!/bin/bash

#SBATCH -J dada2_yellow_stone
#SBATCH -o dada2_yellow_stone_out
#SBATCH -e dada2_yellow_stone_error
#SBATCH -n 8
#SBATCH --mem=100G
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=ALL


module load gcc/5.2.0
module load  R/3.4.3 
module load dada2

module load cutadapt

#unzip 

gunzip /gpfs_home/pfreema1/data/pfreema1/YellowstoneTK101_196/*.fastq.gz

#cut adapt to remove forward, reverse and their reverse complements. 

for i in /gpfs_home/pfreema1/data/pfreema1/YellowstoneTK101_196/*_R1_001.fastq;

do
  SAMPLE=$(echo ${i} | sed "s/_R1_\001\.fastq//") 
  echo ${SAMPLE}_R1_001.fastq ${SAMPLE}_R2_001.fastq
  cutadapt -a GGGCAATCCTGAGCCAA -g TTGGCTCAGGATTGCCC -A CCATTGAGTCTCTGCACCTATC -G GATAGGTGCAGAGACTCAATGG -o  ${SAMPLE}_trimmed_R1_001.fastq -p ${SAMPLE}_trimmed_R2_001.fastq  ${SAMPLE}_R1_001.fastq  ${SAMPLE}_R2_001.fastq 

done

#run dada2 script on trimmed sequences

Rscript --vanilla --max-ppsize=5000000 dada2_yellow_stone.R


```







