#!/usr/bin/env Rscript

library(tidyverse)
library(caret)
DATADIR <- "~/data/ENA_SRP044740/tidydata/"

#### Define functions ---------------

trainRF <- function(X, Y) {
  library(h2o)
  Sys.unsetenv("http_proxy")
  seed.id <- as.numeric(1)
  
  ## Create an H2O cloud 
  h2o.init(
    port = 54321,
    nthreads = 4,            ## -1: use all available threads
    max_mem_size = "40G")    ## specify the memory size for the H2O cloud
  # h2o.removeAll() # Clean slate - just in case the cluster was already running
  
  ## For h2o, we need X and Y matrices to be together in the same matrix
  train <- bind_cols(X, Y = Y)
  
  # Type transform
  train.h2o <- as.h2o(train)
  
  rate.per.class.list <- c(0.6320000291, 0.6320000291)
  
  ## run our first predictive model
  model <- h2o.randomForest(
    training_frame = train.h2o,    ## the H2O frame for training
    x = 1:(ncol(train.h2o)-1),     ## the predictor columns, by column index
    y = ncol(train.h2o),           ## the target index (what we are predicting)
    categorical_encoding = "AUTO",
    model_id = "RF_final_model",    ## name the model in H2O
    sample_rate_per_class = rate.per.class.list,  ## setting this to 1 wont do any sampling 
    ntrees = 500,
    mtries = -1, # If set to -1, defaults to sqrtp for classification and p/3 for regression (where p is the # of predictors Defaults to -1.
    max_depth = 20, # same as nodesize
    stopping_rounds = 0,  
    seed = 1)   ## Set the random seed so that this can be reproduced
  
  h2o.saveModel(model, path = "results/models/RF_final_model", force = TRUE)
}

trainxgbTree <- function(X, Y) {
  set.seed(1)
  
  # these are xgbtree's default parameters
  tunegrid <- expand.grid(nrounds = 100,
                          max_depth = 6,
                          eta = 0.3,
                          gamma = 0,
                          colsample_bytree = 1,
                          min_child_weight = 1,
                          subsample = 1)
  
  fitControl.noCV <- trainControl(
    method = "none",
    summaryFunction = twoClassSummary,
    classProbs = TRUE,
    verboseIter = TRUE)
  
  model <- caret::train(x = as.matrix(X), 
                                y = Y, 
                                method = "xgbTree", 
                                tuneGrid = tunegrid,
                                trControl = fitControl.noCV, 
                                preProcess = c("center", "scale"),
                                metric = "ROC",
                        importance = TRUE)
  
  # save
  saveRDS(model, file = "results/models/xgbTree_final_model.RDS")
}

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

oneHotEncoding <- function(X) {
  # This function takes categorical variables/columns from X matrix and creates multiple binary variables/columns. This is done for learning models that cannot deal with categorical variables
  # Note that the function dummyVars from caret does it for you
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


#### Load data  ---------------

Y.deam.filenames <- list.files(path = DATADIR, pattern = "_deaminations_Y.rds", full.names = TRUE)
Y.deam <- Y.deam.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.)))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(isDeam == "1")

Y.mut.filenames <- list.files(path = DATADIR, pattern = "_real-mutations_FFPE_Y.rds", full.names = TRUE)
Y.mut <- Y.mut.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.))))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(isSomatic == "1" | isSNP == "1")

X.deam.filenames <- list.files(path = DATADIR, pattern = "_deaminations_X.rds", full.names = TRUE)
X.deam <- X.deam.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.)))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(complete_id %in% Y.deam$complete_id)

X.mut.filenames <- list.files(path = DATADIR, pattern = "_real-mutations_FFPE_X.rds", full.names = TRUE)
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

# Check if a mutation is repeated i.e. considered both as somatic mutation and deamination. If it is, remove it from deaminationlist
X %>%
  group_by(sample) %>%
  filter(duplicated(id)) %>%
  pull(complete_id) -> dup.ids

X %>%
  filter(!(complete_id %in% dup.ids) | source != "deam") -> X

Y %>%
  filter(!(complete_id %in% dup.ids) | source != "deam") -> Y

AF.cutoff <- 0.3
X %>%
  filter(allele.freq <= AF.cutoff) -> X
Y %>%
  filter(complete_id %in% X$complete_id) -> Y


# subsample in case it's needed
max.ratio <- 50
deam.mut.ratio <- table(X$source)
if (deam.mut.ratio["deam"] > max.ratio*deam.mut.ratio["mut"]) { # if TRUE subsample deaminations
  cat("Subsampling deaminations\n")
  deam.toKeep <- deam.mut.ratio["mut"]*max.ratio
  subs.XY <- subsampMatrices(X = X, Y = Y, n.keep = deam.toKeep, mut.source = "deam")
  X <- subs.XY[["X"]]
  Y <- subs.XY[["Y"]]
} else if (deam.mut.ratio["deam"] < max.ratio*deam.mut.ratio["mut"]) { # if true subsample mutation
  cat("Subsampling mutations\n")
  mut.toKeep <- ceiling(deam.mut.ratio["deam"]/max.ratio)
  subs.XY <- subsampMatrices(X = X, Y = Y, n.keep = mut.toKeep, mut.source = "mut")
  X <- subs.XY[["X"]]
  Y <- subs.XY[["Y"]]
}

# Pull necessary features
X.cols <- setdiff(colnames(X), c("source", "current", "bases", "read.length", "both.reads.aligned", "sample", "frag.length.frac", "isSNP"))
X <- X[, X.cols]
Y <- Y %>%
  pull(isDeam) %>%
  as.factor()
levels(Y) <- make.names(levels(Y))

# Some filtering and adaption depending on learning model
X.xgbtree <- oneHotEncoding(X)
X.RF <- dplyr::select(X, -ends_with("id"))

#### Train and save final models  ---------------

trainxgbTree(X = X.xgbtree, Y = Y)
trainRF(X = X.RF, Y = Y)
