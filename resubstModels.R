#!/usr/bin/env Rscript

####################################################################
# READ BEFORE USING
# Script needs three arguments: learning model to be used (NB, xgbTree, NN, logReg or RF), id (arbitrary number to be attached to the output filename) and samplename
###################################################################

library(dplyr)
library(tidyr)
library(purrr)
library(lattice)
library(caret)
library(pROC)
require(xgboost)
source("RF_resubst.R")
source("logReg_resubst.R")
source("NB_resubst.R")
source("NN_resubst.R")
source("xgbTree_resubst.R")


args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 3) {
  stop("Only three arguments should be provided: model(NB, logReg, RF, kNN, NN or svm), task_id and samplename")
}

model <- args[1]
id.seed <- args[2]
samplename <- args[3]

# Load data and keep only deaminations and mutations
Y.deam.filenames <- list.files(path = ".", pattern = "_deaminations_Y.rds")
Y.deam.filenames <- Y.deam.filenames[grepl(sprintf("^(?!%s_d)", samplename), Y.deam.filenames, perl = TRUE)]
Y.deam <- Y.deam.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.)))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(isDeam == "1")

Y.mut.filenames <- list.files(path = ".", pattern = "_real-mutations_FFPE_Y.rds")
Y.mut.filenames <- Y.mut.filenames[grepl(sprintf("^(?!%s_re)", samplename), Y.mut.filenames, perl = TRUE)] 
Y.mut <- Y.mut.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.))))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(isSomatic == "1" | isSNP == "1")

X.deam.filenames <- list.files(path = ".", pattern = "_deaminations_X.rds")
X.deam.filenames <- X.deam.filenames[grepl(sprintf("^(?!%s_d)", samplename), X.deam.filenames, perl = TRUE)]
X.deam <- X.deam.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.)))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(complete_id %in% Y.deam$complete_id)

X.mut.filenames <- list.files(path = ".", pattern = "_real-mutations_FFPE_X.rds")
X.mut.filenames <- X.mut.filenames[grepl(sprintf("^(?!%s_re)", samplename), X.mut.filenames, perl = TRUE)] 
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

# We'll limit the analysis to instances with allele-frequency < 0.3. After that we'll subsample deamination or somatic mutation data, so that a ratio of 1:50 is maintained between deaminations and somatic mutations. 
AF.cutoff <- 0.3
X %>%
  filter(allele.freq <= AF.cutoff) -> X
Y %>%
  filter(complete_id %in% X$complete_id) -> Y

# subsampMatrices <- function(X, Y, n.keep, mut.source) {
#   if (nrow(X) != nrow(Y)) {
#     stop("X and Y matrices should have the same number of instances")
#   }
#   idx.to.subsample <- which(X$source == mut.source)
#   idx.to.keep <- setdiff(1:nrow(Y), idx.to.subsample)
#   subsample.idx <- base::sample(idx.to.subsample, n.keep)
#   all.idx.to.keep <- sort(c(idx.to.keep, subsample.idx))
#   X <- X[all.idx.to.keep,]
#   Y <- Y[all.idx.to.keep,]
#   return(list(X = X, Y = Y))
# }


# Pull necessary features
X.cols <- setdiff(colnames(X), c("source", "current", "bases", "read.length", "both.reads.aligned", "sample", "frag.length.frac", "isSNP"))
X <- X[, X.cols]
Y <- Y %>%
  pull(isDeam) %>%
  as.factor()
levels(Y) <- make.names(levels(Y))

# Call models

models.to.funs <- list(
  NB = "trainNB",
  logReg = "trainLogReg",
  RF = "trainRF",
  xgbTree = "trainxgbTree",
  NN = "trainNN"
)

# Some filtering and adaption depending on learning model

oneHotEncoding <- function(X) {
  # This function takes categorical variables/columns from X matrix and creates multiple binary variables/columns. This is done for learning models that cannot deal with categorical variables
  categorical.vars <- X %>%
    select_if(function(x) is.factor(x) | is.character(x)) %>%
    colnames()
  categorical.vars <- setdiff(categorical.vars, c("id", "complete_id"))
  lapply(categorical.vars, function(var) {
    X %>%
      mutate(yesno = 1) %>%
      spread(var, yesno, fill = 0, sep = "_")
  }) -> dedup.X.by.var
  dedup.X.by.var %>%
    purrr::reduce(.f = full_join) -> dedup.X
  dedup.X %>%
    select_if(is.numeric) -> dedup.X # in this final step we discard id/complete_id column
  dedup.X
}

# Some filtering and adaption depending on learning model
X <- if(model %in% c("NN", "xgbTree")) {
  oneHotEncoding(X)
} else {
  dplyr::select(X, -ends_with("id"))
}

# Call learning model
fun <- get(models.to.funs[[model]])
fun(X, Y, id.seed)
