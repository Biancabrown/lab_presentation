#!/usr/bin/env Rscript

library(dada2); packageVersion("dada2")
library(Rcpp); packageVersion("Rcpp")
library(ShortRead); packageVersion("ShortRead")

path <- "path_to_your_file"


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,trimLeft = c(), truncLen=c(),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
saveRDS(out, "out.rds")

errF <- learnErrors(filtFs, multithread=TRUE)
saveRDS(errF, "errF.rds")

errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errR, "errR.rds")

derepFs <- derepFastq(filtFs, verbose=TRUE)
saveRDS(derepFs, "derepFs.rds")

derepRs <- derepFastq(filtRs, verbose=TRUE)
saveRDS(derepRs, "derepRs.rds")

#Name the derep-class objects by the sample names

names(derepFs) <- sample.names

names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
saveRDS(dadaFs, "dadaFs.rds")

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
saveRDS(dadaRs, "dadaRs.rds")

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
saveRDS(mergers, "mergers.rds")

#Inspect the merger data.frame from the first sample

seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "seqtab.rds")

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim, "seqtab.nochim.rds")

