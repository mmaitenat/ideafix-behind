#!/usr/bin/env Rscript
prepareDataMatrices <- function(samplename,
                                outfolder,
                                deam.filename = NULL,
                                somatic.mut.filename = NULL,
                                somatic.mut.source = NULL,
                                somatic.mut.HRun.filename = NULL,
                                somatic.mut.Polyx.filename = NULL,
                                deam.HRun.filename = NULL,
                                deam.Polyx.filename = NULL,
                                read.length = 100) {
  ### This function creates X and Y matrices from a real-data vcf file, be it deaminations or somatic mutations.
  
  library(tidyverse)
  library(logging)
  source("ExtractFromVcf_bcftools.R")
  source("ExtractFromFa.R")
  source("getHomopolymers.R")
  
  #### Log
  basicConfig()
  log.filename <- file.path(outfolder, "prepareDataMatrices.log")
  addHandler(writeToFile, file = log.filename)
  invokement <- paste(R.utils::commandArgs(), collapse=" ")
  loginfo("How R was invoked:\n", logger = "")
  loginfo(invokement, logger = "")
  
  #### Helper functions
  
  getDeamPosFromFile <- function(deam_filename) {
    data <- as_tibble(read.table(deam_filename, sep = ":")[, c(2:5, 8)])
  }
  
  fixFragmentLength0 <- function(X.mat) {
    X.mat %>%
      mutate(norm.pos.from.end = replace(norm.pos.from.end, both.reads.aligned == 0, NA),
             frag.length = replace(frag.length, both.reads.aligned == 0, NA),
             frag.length.frac = replace(frag.length.frac, both.reads.aligned == 0, NA)) -> X.mat
  }
  
  #### Common data
  # Change here accordingly
  fasta.filename <- "~/genomes/hg19/ucsc.hg19.fasta"
  
  #### Obtain X and Y for the provided arguments
  if (!is.null(deam.filename)) {
    
    # Log
    log.msg <- sprintf("Deamination X matrix will be created from file %s and sent to %s. Y will also be created.", deam.filename, file.path(outfolder, sprintf("%s_deaminations_X.rds", samplename)))
    loginfo(log.msg, logger = "")
    
    # Get sample.id
    sample.id <- sub('\\_.*', '', sub('\\.*', '', basename(deam.filename)))
    
    # Load X
    X.deam <- ExtractFromVcf(vcf.filename = deam.filename, 
                             sample.id = sample.id,
                             read.length = read.length,
                             is.paired = TRUE,
                             allele.freq = TRUE, 
                             alt.bases = TRUE, 
                             ref.bases = TRUE, 
                             ref.allele = TRUE, 
                             alt.allele = TRUE, 
                             base.qual = TRUE, 
                             base.qual.frac = TRUE, 
                             frag.length = TRUE, 
                             pos.from.end = TRUE, 
                             map.qual = TRUE,
                             FdeamC = TRUE,
                             SB = TRUE)
    
    out.fasta.deam <- ExtractFromFa(vcf.filename = deam.filename,
                                    fa.filename = fasta.filename,
                                    k = 2)
    homopolymer.deam <- getHomopolymers(HRun.vcf = deam.HRun.filename,
                                        VCFPolyx.vcf = deam.Polyx.filename)
    
    X.deam <- full_join(X.deam, out.fasta.deam, by = "id")
    X.deam <- left_join(X.deam, homopolymer.deam, by = "id")
    
    # Filter out those that are not C:G>T:A
    X.deam %>%
      mutate(mut = paste(ref.allele, alt.allele, sep = ":")) %>%
      filter(mut %in% c("C:T", "G:A")) %>%
      select(-mut) -> X.deam
    
    # Fix according to unknown fragment length
    X.deam <- fixFragmentLength0(X.deam)
    
    # Homopolymers: replace NA (not homopolymer) with 0 and discard those called by Polyx that do not meet our defined requirements
    X.deam %>%
      mutate(hp.length = replace_na(hp.length, 0)) %>%
      mutate(hp.length = ifelse(hp.source %in% "Polyx" & alt.allele != toupper(before.1) & alt.allele != toupper(after.1), 0, hp.length)) %>%
      select(-hp.source) -> X.deam
    
    # Create Y
    Y.deam <- tibble(id = X.deam$id, isDeam = "1", isSomatic = "0", isSNP = "0") 
    
    X.deam %>%
      filter(both.reads.aligned == 0) %>%
      pull(id) -> frag.zero.deam.ids
    
    Y.deam %>%
      mutate(isNoise = ifelse(id %in% frag.zero.deam.ids, "1", "0")) %>%
      mutate(isDeam = replace(isDeam, isNoise == "1", "0")) -> Y.deam
    
    col.idx <- 2:5
    Y.deam <- mutate_at(Y.deam, col.idx, funs(as.factor))
    
    # Save data
    saveRDS(X.deam, file = file.path(outfolder, sprintf("%s_deaminations_X.rds", samplename)))
    saveRDS(Y.deam, file = file.path(outfolder, sprintf("%s_deaminations_Y.rds", samplename)))
  }
    
  if (!is.null(somatic.mut.filename)) {
    if (is.null(somatic.mut.source)) {
      stop("Somatic mutation vcf filename source (FF or FFPE) should be provided")
    }
    
    # Log
    log.msg <- sprintf("Somatic mutation X matrix will be created from file %s and sent to %s. Y will also be created.", somatic.mut.filename, file.path(outfolder, sprintf("%s_real-mutations_%s_X.rds", samplename, somatic.mut.source)))
    loginfo(log.msg, logger = "")
    
    # Get sample.id
    sample.id <- sub('\\_.*', '', basename(somatic.mut.HRun.filename)) # could also use VCFPolyx filename
    
    # Load X
    X.mut <- ExtractFromVcf(vcf.filename = somatic.mut.filename, 
                            sample.id = sample.id,
                            read.length = read.length,
                            is.paired = FALSE,
                            allele.freq = TRUE, 
                            alt.bases = TRUE, 
                            ref.bases = TRUE, 
                            ref.allele = TRUE, 
                            alt.allele = TRUE, 
                            base.qual = TRUE, 
                            base.qual.frac = TRUE, 
                            frag.length = TRUE, 
                            pos.from.end = TRUE, 
                            map.qual = TRUE,
                            FdeamC = TRUE,
                            SB = TRUE)
    
    out.fasta.mut <- ExtractFromFa(vcf.filename = somatic.mut.filename,
                                   fa.filename = fasta.filename,
                                   k = 2)
    
    homopolymer.mut <- getHomopolymers(HRun.vcf = somatic.mut.HRun.filename,
                                        VCFPolyx.vcf = somatic.mut.Polyx.filename)
    
    X.mut <- full_join(X.mut, out.fasta.mut, by = "id")
    X.mut <- left_join(X.mut, homopolymer.mut, by = "id")
    
    # Fix according to unknown fragment length
    X.mut <- fixFragmentLength0(X.mut)
    
    # Homopolymers: replace NA (not homopolymer) with 0 and discard those called by Polyx that do not meet our defined requirements 
    X.mut %>%
      mutate(hp.length = replace_na(hp.length, 0)) %>%
      mutate(hp.length = ifelse(hp.source %in% "Polyx" & alt.allele != toupper(before.1) & alt.allele != toupper(after.1), 0, hp.length)) %>%
      select(-hp.source) -> X.mut
    
    # Create Y
    Y.mut <- tibble(id = X.mut$id, isDeam = "0", isSomatic = "1")
    
    X.mut %>%
      filter(both.reads.aligned == 0) %>%
      pull(id) -> frag.zero.mut.ids
    
    Y.mut %>%
      mutate(isNoise = ifelse(id %in% frag.zero.mut.ids, "1", "0")) %>%
      mutate(isSNP = X.mut$isSNP) %>%
      mutate(isNoise = replace(isNoise, isSNP == "1", "0")) %>%
      mutate(isSomatic = replace(isSomatic, isSNP == "1" | isNoise == "1", "0")) -> Y.mut
    
    col.idx <- 2:5
    Y.mut <- mutate_at(Y.mut, col.idx, funs(as.factor))
    
    # Save data
    saveRDS(X.mut, file = file.path(outfolder, sprintf("%s_real-mutations_%s_X.rds", samplename, somatic.mut.source)))
    saveRDS(Y.mut, file = file.path(outfolder, sprintf("%s_real-mutations_%s_Y.rds", samplename, somatic.mut.source)))
  }
}

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
prepareDataMatrices(samplename = args$samplename, 
                    deam.filename = args$deam.filename, 
                    somatic.mut.filename = args$somatic.mut.filename, 
                    somatic.mut.source = args$somatic.mut.source, 
                    somatic.mut.HRun.filename = args$somatic.mut.HRun.filename,
                    somatic.mut.Polyx.filename = args$somatic.mut.Polyx.filename,
                    deam.HRun.filename = args$deam.HRun.filename,
                    deam.Polyx.filename = args$deam.Polyx.filename,
                    outfolder = args$outfolder, 
                    read.length = ifelse(is.null(args$read.length), 100, args$read.length))