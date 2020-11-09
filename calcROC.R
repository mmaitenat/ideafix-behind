#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(purrr)

### Helper functions

calcFeaturePower <- function(feature.value, labels) {
  library(caret)
  ### Helper functions
  ROCarea <- function(tpr, fpr) { # area del trapezoide
    height <- (tpr[-1] + tpr[-length(tpr)])/2
    width <- -diff(fpr) 
    sum(height*width, na.rm = TRUE)
  }
  
  ### Run code
  data <- bind_cols(x = feature.value, y = labels)
  
  data %>%
    arrange(x) %>%
    pull(x) %>%
    unique() -> thresholds
  
  thresholds %>%
    set_names(as.character(1:length(thresholds))) %>% # need to do this for using map_dfr
    map_dfr(function(thr) {
      pred.labels <- ifelse(data$x > thr, 0, 1) %>%
        factor(levels = c("0", "1"))
      out <- tibble(tpr = sensitivity(pred.labels, data$y),
                    fpr = 1 - specificity(pred.labels, data$y))
    }) -> roc.data
  AUC_value <- ROCarea(tpr = roc.data$tpr, fpr = roc.data$fpr)
  return(AUC_value)
}

### Run code

# Parse arguments

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1) {
  stop("Only one argument should be provided: feature name to calculate power of")
}

feature <- args[1]

features_to_colnames <- list(
  allele_freq = "allele.freq",
  alt_bases = "alt.bases",
  norm_alt_bases = "norm.alt.bases",
  ref_bases = "ref.bases",
  norm_ref_bases = "norm.ref.bases",
  base_qual = "base.qual",
  base_qual_frac = "base.qual.frac",
  frag_length = "frag.length",
  frag_length_frac = "frag.length.frac",
  pos_from_end = "pos.from.end",
  map_qual = "map.qual",
  FdeamC = "FdeamC",
  SOB = "SOB",
  SBGuo = "SBGuo",
  SBGATK = "SBGATK",
  norm_pos_from_end = "norm.pos.from.end",
  is_repeat_region = "is.repeat.region",
  hp_length = "hp.length"
)

`%nin%` = Negate(`%in%`)
if (feature %nin% names(features_to_colnames)) {
  stop(sprintf("Feature should be one of %s", paste(names(features_to_colnames), collapse = " ")))
}


# Calculate power

Y.deam.filenames <- list.files(path = ".", pattern = "_deaminations_Y.rds", full.names = TRUE)
Y.deam <- Y.deam.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.)))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(isDeam == "1")

Y.mut.filenames <- list.files(path = ".", pattern = "_real-mutations_FFPE_Y.rds", full.names = TRUE)
Y.mut <- Y.mut.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.))))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(isSomatic == "1" | isSNP == "1")

X.deam.filenames <- list.files(path = ".", pattern = "_deaminations_X.rds", full.names = TRUE)
X.deam <- X.deam.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.)))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(complete_id %in% Y.deam$complete_id)

X.mut.filenames <- list.files(path = ".", pattern = "_real-mutations_FFPE_X.rds", full.names = TRUE)
X.mut <- X.mut.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.))))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(complete_id %in% Y.mut$complete_id)

X <- bind_rows(list(deam = X.deam, mut = X.mut), .id = "source") %>% mutate_if(is.character, as.factor)
Y <- bind_rows(list(deam = Y.deam, mut = Y.mut), .id = "source")

# Tidy X and Y matrices
# Keep only C>T/G>A.
nonCT.idx <- which(is.na(X), arr.ind = TRUE)[,1]
X <- X[-nonCT.idx,] %>%
  droplevels()
Y <- Y[-nonCT.idx,]

# Check if a mutation is repeated i.e. considered both as somatic mutation and deamination. If it is, remove it from deamination list
X %>%
  group_by(sample) %>%
  filter(duplicated(id)) %>%
  pull(complete_id) -> dup.ids

X %>%
  filter(!(complete_id %in% dup.ids) | source != "deam") -> X

Y %>%
  filter(!(complete_id %in% dup.ids) | source != "deam") -> Y

# Frequency filter
AF.cutoff <- 0.3
X %>%
  filter(allele.freq <= AF.cutoff) -> X
Y %>%
  filter(complete_id %in% X$complete_id) -> Y

# Pull necessary features
X %>%
  pull(features_to_colnames[[feature]]) -> feature_values

miniY <- Y %>%
  pull(isDeam) %>%
  factor(levels = c("0", "1"))

AUC_power_feature_allsamples <- calcFeaturePower(feature.value = feature_values, labels = miniY)
