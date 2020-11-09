#!/usr/bin/env Rscript

##################################################################################################
# In this script I am comparing my model with 3 tools or practices: 
# 1. Using coverage thresholds
# 2. Using VAF thresholds
# 3. SOBDetector
##################################################################################################

library(tidyverse)
library(h2o)
library(caret)
library(xgboost)
library(ggsci)

# Load test data --------------

test_datadir <- "~/data/colon_liver/tidydata"

Y.deam.filenames <- list.files(path = test_datadir, pattern = "_deaminations_Y.rds", full.names = TRUE)
Y.deam <- Y.deam.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.)))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(isDeam == "1")

Y.mut.filenames <- list.files(path = test_datadir, pattern = "_real-mutations_FFPE_Y.rds", full.names = TRUE)
Y.mut <- Y.mut.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.))))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(isSomatic == "1" | isSNP == "1")

X.deam.filenames <- list.files(path = test_datadir, pattern = "_deaminations_X.rds", full.names = TRUE)
X.deam <- X.deam.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.)))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(complete_id %in% Y.deam$complete_id)

X.mut.filenames <- list.files(path = test_datadir, pattern = "_real-mutations_FFPE_X.rds", full.names = TRUE)
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

# Pull necessary features
X.cols <- setdiff(colnames(X), c("source", "current", "bases", "read.length", "both.reads.aligned", "sample", "frag.length.frac", "isSNP"))
X <- X[, X.cols]
Y <- Y %>%
  pull(isDeam) %>%
  as.factor()
levels(Y) <- make.names(levels(Y))


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

X.xgbTree <- oneHotEncoding(X)
X.RF <- dplyr::select(X, -ends_with("id"))

## Ideafix: XGBoost and RF ----

ROCarea <- function(tpr, fpr) { # area del trapezoide
  height <- (tpr[-1] + tpr[-length(tpr)])/2
  width <- diff(fpr)
  sum(height*width, na.rm = TRUE)
}

model_datadir <- "results/models/"
xgbTree_model <- readRDS(file.path(model_datadir, "xgbTree_final_model.RDS"))
Sys.unsetenv("http_proxy")
h2o.init(
  port = 54325,
  nthreads = 2,  
  max_mem_size = "20G") 
RF_model <- h2o.loadModel(file.path(model_datadir, "RF_final_model"))

## XGBoost 
xgbTree_pred <- predict(xgbTree_model, as.matrix(X.xgbTree), type = "prob")
xgbTree_out <- pROC::roc(Y, xgbTree_pred[, 1], plot = FALSE)
ROC_XGBoost <- ROCarea(tpr = rev(xgbTree_out$sensitivities), fpr = 1-rev(xgbTree_out$specificities))

## RF
test <- bind_cols(X.RF, Y = Y)
test.h2o <- as.h2o(test)
RF_pred <- h2o.predict(object = RF_model, 
                       newdata = test.h2o)
RF_out <- pROC::roc(Y, as.data.frame(RF_pred)$X1, plot = FALSE)
ROC_RF <- ROCarea(tpr = rev(RF_out$sensitivities), fpr = 1-rev(RF_out$specificities))

## Coverage and VAF filters

X %>%
  mutate(depth = alt.bases + ref.bases) -> X_with_depth

X_with_depth %>%
  arrange(depth) %>%
  pull(depth) %>%
  unique() -> depth_thresholds

depth_thresholds <- c(0, depth_thresholds)

X %>%
  arrange(allele.freq) %>%
  pull(allele.freq) %>%
  unique() -> VAF_thresholds

VAF_thresholds <- c(0, VAF_thresholds)

depth_thresholds %>%
  set_names(as.character(1:length(depth_thresholds))) %>% # need to do this for using map_dfr
  map_dfr(function(thr) {
    pred.labels <- ifelse(X_with_depth$depth > thr, "X1", "X0") %>%
      factor(levels = c("X0", "X1"))
    out <- tibble(tpr = sensitivity(data = pred.labels, 
                                    reference = Y),
                  fpr = 1 - specificity(data = pred.labels, 
                                        reference = Y))
  }) -> depth_roc_data

ROC_depth <- ROCarea(tpr = depth_roc_data$tpr,
                     fpr = depth_roc_data$fpr)

VAF_thresholds %>%
  set_names(as.character(1:length(VAF_thresholds))) %>% # need to do this for using map_dfr
  map_dfr(function(thr) {
    pred.labels <- ifelse(X$allele.freq > thr, "X0", "X1") %>%
      factor(levels = c("X0", "X1"))
    out <- tibble(tpr = sensitivity(data = pred.labels, 
                                    reference = Y),
                  fpr = 1 - specificity(data = pred.labels, 
                                        reference = Y))
  }) -> vaf_roc_data

ROC_vaf <- ROCarea(tpr = vaf_roc_data$tpr,
                   fpr = vaf_roc_data$fpr)

## SOBDetector ----
# Prepare input data
# Need to create a data table for each sample
SOBDetector_dir <- "results/SOBDetector/"
bamdir <- "~/data/colon_liver/"

# To be run only once
X %>%
  separate(complete_id, into = c("sample", "chr", "pos"), sep = ":") %>%
  select(sample, chr, pos, ref.allele, alt.allele) %>%
  group_split(sample) %>%
  walk(., function(x) {
    samplename = x$sample[1]
    write.table(select(x, -sample),
                file = file.path(SOBDetector_dir, paste0(samplename, ".table")),
                quote = FALSE,
                sep = '\t',
                row.names = FALSE)
  })

# Also to be run only once
# Change SOBDetector path accordingly
call_SOBDetector <- function(samplename) {
  infilename <- file.path(SOBDetector_dir, paste0(samplename, ".table"))
  bamfilename <- file.path(bamdir, samplename, "aligned", paste0(samplename, ".bam"))
  outfilename <- file.path(SOBDetector_dir, "output", paste0(samplename, "_SOB.out"))
  cmd <- sprintf("java -jar ~/SOBDetector/SOBDetector_v1.0.2.jar --input-type Table --input-variants %s --input-bam %s --output-variants %s", infilename, bamfilename, outfilename)
  print(cmd)
  system(cmd, intern = FALSE)
}

# Read SOBDetector data
list.files(file.path(SOBDetector_dir, "output"), full.names = TRUE) %>%
  set_names(sub("_[^_]+$", "", basename(.))) %>%
  map_dfr(read_table2, col_names = TRUE, col_types = cols(.default = "c"), .id = "sample") %>%
  mutate(complete_id = paste(sample, chr, pos, sep = ":")) -> SOBDetector_out

# I discarded ids from Y vector so I can't obtain real labels. Need to rebuild Y hence
Y <- bind_rows(list(deam = Y.deam, mut = Y.mut), .id = "source")
Y <- Y[-nonCT.idx,]
Y %>%
  filter(!(complete_id %in% dup.ids) | source != "deam") -> Y
Y %>%
  filter(complete_id %in% X$complete_id) -> Y
Y %>%
  left_join(SOBDetector_out, by = "complete_id") %>%
  select(isDeam, pArtifact, artiStatus) -> SOBDetector_predictions
SOBDetector_predictions %>%
  mutate(isDeam = factor(isDeam),
         pArtifact = as.numeric(pArtifact)) -> SOBDetector_predictions

# ROC curve
pArtifacts_thresholds <- SOBDetector_predictions %>%
  select(pArtifact) %>%
  arrange(pArtifact) %>%
  pull(pArtifact) %>%
  unique()

pArtifacts_thresholds <- c(0, pArtifacts_thresholds)

pArtifacts_thresholds <- pArtifacts_thresholds[!is.na(pArtifacts_thresholds)]
pArtifacts_thresholds %>%
  set_names(as.character(1:length(pArtifacts_thresholds))) %>% # need to do this for using map_dfr
  map_dfr(function(thr) {
    pred.labels <- ifelse(SOBDetector_predictions$pArtifact > thr, "1", "0") %>%
      factor(levels = c("0", "1"))
    out <- tibble(tpr = sensitivity(data = pred.labels, 
                                    reference = SOBDetector_predictions$isDeam),
                  fpr = 1 - specificity(data = pred.labels, 
                                        reference = SOBDetector_predictions$isDeam))
  }) -> SOBDetector_roc_data

ROC_SOBDetector <- ROCarea(tpr = SOBDetector_roc_data$tpr,
                           fpr = SOBDetector_roc_data$fpr)



# ROC curve plots
xgbTree_roc_data <- data.frame(tpr = rev(xgbTree_out$sensitivities), fpr = 1-rev(xgbTree_out$specificities))
RF_roc_data <- data.frame(tpr = rev(RF_out$sensitivities), fpr = 1-rev(RF_out$specificities))
vaf_roc_data <- data.frame(tpr = rev(vaf_roc_data$tpr), fpr = rev(vaf_roc_data$fpr))

all_ROC_values <- bind_rows(VAF = vaf_roc_data, 
                            depth = depth_roc_data, 
                            SOBDetector = SOBDetector_roc_data,
                            XGBoost = xgbTree_roc_data, 
                            RF = RF_roc_data,
                            .id = "model") 
ggplot(all_ROC_values, aes(x = fpr, y = tpr, group = model, colour = model)) +
  theme_bw() +
  scale_color_npg() +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  geom_step() +
  labs(x = "1 - Specificity",
       y = "Sensitivity") -> roc_curves

plot_outdir <- "results/figs/"
pdf(file.path(plot_outdir, "fig7.pdf"))
plot(roc_curves)
dev.off()