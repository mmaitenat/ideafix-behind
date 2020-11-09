#!/usr/bin/env Rscript

####################################################################
# READ BEFORE USING
# Script needs four arguments: learning model to be used (NB, xgbTree, NN, logReg or RF), id (arbitrary number to be attached to the output filename), train samplename and test samplename
###################################################################

library(dplyr)
library(tidyr)
library(purrr)
library(lattice)
library(caret)
library(pROC)
require(xgboost)
source("RF.R")
source("logReg.R")
source("NB.R")
source("NN.R")
source("xgbTree.R")


args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 4) {
  stop("Only four arguments should be provided: model(NB, logReg, RF, NN or xgbTree), task_id, train samplename and test samplename")
}

model <- args[1]
id.seed <- args[2]
train_samplename <- args[3]
test_samplename <- args[4]

# Load data and keep only deaminations and mutations
Y.deam.filenames <- list.files(path = ".", pattern = "_deaminations_Y.rds")
Y.deam.filenames <- Y.deam.filenames[grepl(sprintf("^(?!%s_d)", test_samplename), Y.deam.filenames, perl = TRUE)]
Y.deam.train <- Y.deam.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.)))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(isDeam == "1")

Y.mut.filenames <- list.files(path = ".", pattern = "_real-mutations_FFPE_Y.rds")
Y.mut.filenames <- Y.mut.filenames[grepl(sprintf("^(?!%s_re)", test_samplename), Y.mut.filenames, perl = TRUE)] 
Y.mut.train <- Y.mut.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.))))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(isSomatic == "1" | isSNP == "1")

X.deam.filenames <- list.files(path = ".", pattern = "_deaminations_X.rds")
X.deam.filenames <- X.deam.filenames[grepl(sprintf("^(?!%s_d)", test_samplename), X.deam.filenames, perl = TRUE)]
X.deam.train <- X.deam.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.)))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(complete_id %in% Y.deam.train$complete_id)

X.mut.filenames <- list.files(path = ".", pattern = "_real-mutations_FFPE_X.rds")
X.mut.filenames <- X.mut.filenames[grepl(sprintf("^(?!%s_re)", test_samplename), X.mut.filenames, perl = TRUE)] 
X.mut.train <- X.mut.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.))))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(complete_id %in% Y.mut.train$complete_id)

X.train <- bind_rows(list(deam = X.deam.train, mut = X.mut.train), .id = "source") %>% mutate_if(is.character, as.factor)
Y.train <- bind_rows(list(deam = Y.deam.train, mut = Y.mut.train), .id = "source")

# Tidy X and Y matrices
# Keep only C>T/G>A.
nonCT.idx.train <- which(is.na(X.train), arr.ind = TRUE)[,1]
X.train <- X.train[-nonCT.idx.train,] %>%
  droplevels()
Y.train <- Y.train[-nonCT.idx.train,]

# Check if a mutation is repeated i.e. considered both as somatic mutation and deamination. If it is, remove it from deamination list
X.train %>%
  group_by(sample) %>%
  filter(duplicated(id)) %>%
  pull(complete_id) -> train.dup.ids

X.train %>%
  filter(!(complete_id %in% train.dup.ids) | source != "deam") -> X.train

Y.train %>%
  filter(!(complete_id %in% train.dup.ids) | source != "deam") -> Y.train

# We'll limit the analysis to instances with allele-frequency < 0.3. After that we'll subsample deamination or somatic mutation data, so that a ratio of 1:50 is maintained between deaminations and somatic mutations. 
AF.cutoff <- 0.3
X.train %>%
  filter(allele.freq <= AF.cutoff) -> X.train
Y.train %>%
  filter(complete_id %in% X.train$complete_id) -> Y.train

subsampMatrices <- function(X, Y, n.keep, mut.source) {
  if (nrow(X) != nrow(Y)) {
    stop("X and Y matrices should have the same number of instances")
  }
  idx.to.subsample <- which(X$source == mut.source)
  idx.to.keep <- setdiff(1:nrow(Y), idx.to.subsample)
  subsample.idx <- base::sample(idx.to.subsample, n.keep)
  all.idx.to.keep <- sort(c(idx.to.keep, subsample.idx))
  X <- X[all.idx.to.keep,]
  Y <- Y[all.idx.to.keep,]
  return(list(X = X, Y = Y))
}

# subsample in case it's needed
max.ratio <- 50
deam.mut.ratio.train <- table(X.train$source)
if (deam.mut.ratio.train["deam"] > max.ratio*deam.mut.ratio.train["mut"]) { # if TRUE subsample deaminations
  cat("Train data: Subsampling deaminations\n")
  deam.toKeep.train <- deam.mut.ratio.train["mut"]*max.ratio
  subs.train.XY <- subsampMatrices(X = X.train, Y = Y.train, n.keep = deam.toKeep.train, mut.source = "deam")
  X.train <- subs.train.XY[["X"]]
  Y.train <- subs.train.XY[["Y"]]
} else if (deam.mut.ratio.train["deam"] < max.ratio*deam.mut.ratio.train["mut"]) { # if true subsample mutation
  cat("Train data: Subsampling mutations\n")
  mut.toKeep.train <- ceiling(deam.mut.ratio.train["deam"]/max.ratio)
  subs.train.XY <- subsampMatrices(X = X.train, Y = Y.train, n.keep = mut.toKeep.train, mut.source = "mut")
  X.train <- subs.train.XY[["X"]]
  Y.train <- subs.train.XY[["Y"]]
}

# Pull necessary features
X.cols <- setdiff(colnames(X.train), c("source", "current", "bases", "read.length", "both.reads.aligned", "sample", "frag.length.frac", "isSNP"))
X.train <- X.train[, X.cols]
Y.train <- Y.train %>%
  pull(isDeam) %>%
  as.factor()
levels(Y.train) <- make.names(levels(Y.train))

# Test data
# Load data and keep only deaminations and mutations
Y.deam.test <- readRDS(sprintf("%s_deaminations_Y.rds", test_samplename)) %>%
  filter(isDeam == "1")
Y.mut.test <- readRDS(sprintf("%s_real-mutations_FFPE_Y.rds", test_samplename)) %>%
  filter(isSomatic == "1" | isSNP == "1")
X.deam.test <- readRDS(sprintf("%s_deaminations_X.rds", test_samplename)) %>%
  filter(id %in% Y.deam.test$id)
X.mut.test <- readRDS(sprintf("%s_real-mutations_FFPE_X.rds", test_samplename))  %>%
  filter(id %in% Y.mut.test$id)

X.test <- bind_rows(list(deam = X.deam.test, mut = X.mut.test), .id = "source") %>% mutate_if(is.character, as.factor)
Y.test <- bind_rows(list(deam = Y.deam.test, mut = Y.mut.test), .id = "source")

# Tidy X and Y matrices
# Keep only C>T/G>A.
nonCT.idx.test <- which(is.na(X.test), arr.ind = TRUE)[,1]
X.test <- X.test[-nonCT.idx.test,] %>%
  droplevels()
Y.test <- Y.test[-nonCT.idx.test,]

# Check if a mutation is repeated i.e. considered both as somatic mutation and deamination. If it is, remove it from deamination list
test.dup.ids <- X.test$id[duplicated(X.test$id)]
X.test %>%
  filter(!(id %in% test.dup.ids) | source != "deam") -> X.test
Y.test %>%
  filter(!(id %in% test.dup.ids) | source != "deam") -> Y.test

# Frequency filtering
# We'll limit the analysis to instances with allele-frequency < 0.3.

AF.cutoff <- 0.3
X.test %>%
  filter(allele.freq <= AF.cutoff) -> X.test
Y.test %>%
  filter(id %in% X.test$id) -> Y.test

# We wont subsample test data

# Pull necessary features
X.cols.test <- setdiff(colnames(X.test), c("source", "current", "bases", "read.length", "both.reads.aligned", "sample", "frag.length.frac", "isSNP"))
X.test <- X.test[, X.cols.test]
Y.test <- Y.test %>%
  pull(isDeam) %>%
  as.factor()
levels(Y.test) <- make.names(levels(Y.test))

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

# Train data
X.train <- if(model %in% c("NN", "xgbTree")) {
  oneHotEncoding(X.train)
} else {
  dplyr::select(X.train, -ends_with("id"))
}

# Test data
X.test <- if(model %in% c("NN", "xgbTree")) {
  oneHotEncoding(X.test)
} else {
  dplyr::select(X.test, -ends_with("id"))
}


# Call learning model
fun <- get(models.to.funs[[model]])
fun(X.train, Y.train, X.test, Y.test, id.seed)
