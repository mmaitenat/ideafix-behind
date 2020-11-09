#!/usr/bin/env Rscript

library(tidyverse)
library(ggsci)

datadir <- "results/"
sample_list <- list.dirs(datadir, full.names = TRUE)
names(sample_list) <- basename(sample_list)
sample_list <- sample_list[grepl("sample.*[0-9]$", sample_list)]
plot_outdir <- "figs/"

# Load data ----

# Test error

test_list <- sample_list[grepl("train.*test", sample_list)]
test_data <- list.files(path = test_list, pattern = "*.test-error.*.RDS", full.names = TRUE) %>%
  set_names(sprintf("%s-%s", basename(dirname(.)), basename(.))) %>%
  map(readRDS) %>%
  map_chr("auc") %>%
  enframe() %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = gsub("_test-error.*", "", name)) %>%
  separate(name, into = c("sample", "model"), sep = "-") %>%
  mutate(sample = gsub("_[^_]+$", "", sample)) %>%
  rename(ROC = value)

# Resubstitution error

resubst_data <- list.files(path = sample_list, pattern = "*.resubst-error.*.RDS", full.names = TRUE) %>%
  set_names(sprintf("%s-%s", basename(dirname(.)), basename(.))) %>%
  map(readRDS) %>%
  map_chr("auc") %>%
  enframe() %>%
  mutate(value = as.numeric(value)) %>%
  mutate(name = gsub("_resubst-error.*", "", name)) %>%
  separate(name, into = c("sample", "model"), sep = "-") %>%
  mutate(sample = gsub("_[^_]+$", "", sample)) %>%
  rename(ROC = value)

# kcv error
# We need to treat RF and the rest of models separately
## Non-RF
kcv_data_nonRF <- list.files(path = sample_list, pattern = "[xgbTree|logReg|NB|NN]_model_[0-9].*.RDS", full.names = TRUE) %>%
  set_names(sprintf("%s-%s", basename(dirname(.)), basename(.))) %>%
  map_dfr(readRDS, .id = "name") %>%
  mutate(name = gsub("_model.*", "", name)) %>%
  separate(name, into = c("sample", "model"), sep = "-") %>%
  mutate(sample = gsub("_[^_]+$", "", sample))

## RF
kcv_data_RF <- list.files(path = sample_list, pattern = "RF_model_[0-9].*.RDS", full.names = TRUE) %>%
  set_names(sprintf("%s-%s", basename(dirname(.)), basename(.))) %>%
  map_dfr(function(x) {
    data <- readRDS(x)
    data@model$cross_validation_metrics_summary["auc",][-(1:2)]
  }, .id = "name") %>%
  gather(Resample, ROC, cv_1_valid:cv_10_valid) %>%
  mutate(ROC = as.numeric(ROC)) %>%
  mutate(Resample = case_when(
    Resample == "cv_1_valid" ~ "Fold01",
    Resample == "cv_2_valid" ~ "Fold02",
    Resample == "cv_3_valid" ~ "Fold03",
    Resample == "cv_4_valid" ~ "Fold04",
    Resample == "cv_5_valid" ~ "Fold05",
    Resample == "cv_6_valid" ~ "Fold06",
    Resample == "cv_7_valid" ~ "Fold07",
    Resample == "cv_8_valid" ~ "Fold08",
    Resample == "cv_9_valid" ~ "Fold09",
    Resample == "cv_10_valid" ~ "Fold10")) %>%
  mutate(name = gsub("_model.*", "", name)) %>%
  separate(name, into = c("sample", "model"), sep = "-") %>%
  mutate(sample = gsub("_[^_]+$", "", sample))

## Join data

kcv_data <- bind_rows(kcv_data_nonRF, kcv_data_RF)

## all data

all_data <- bind_rows(resubst_data, select(kcv_data, ROC, sample, model), .id = "origin")
all_data %>%
  mutate(origin = ifelse(origin == "1", "resubstitution", "kcv")) -> all_data
test_data %>%
  mutate(origin = "sample-out") -> test_data
all_data <- bind_rows(all_data, test_data)

all_data %>%
  mutate(idx = gsub("sample", "", sample) %>%
           as.numeric())  %>%
  mutate(model = ifelse(model == "xgbTree", "XGBoost", model)) %>%
  arrange(idx) %>%
  mutate(sample = factor(sample, levels = unique(sample))) %>%
  select(-idx) -> all_data

# Numeric values ----

all_data %>%
  group_by(model, origin) %>%
  summarize(median(ROC))

# Plots ----

## Figure 5
all_data %>%
  filter(origin == "kcv") %>%
  ggplot(aes(x = model, y = ROC, fill = model, colour = model)) + 
  geom_violin() +
  theme_bw() +
  scale_fill_npg() +
  scale_colour_npg() +
  ylab("kcv AUC (CVE)") +
  xlab("") +
  theme(legend.position = "none",
        legend.title = element_blank()) -> fig5
pdf(file = file.path(plot_outdir, "fig5.pdf"), width = 4, height = 2)
print(fig5)
dev.off()

# Figure 4
all_data %>%
  group_by(model, origin, sample) %>%
  summarize(AUC = median(ROC)) %>%
  ungroup() %>%
  ggplot(aes(x = model, y = AUC, fill = origin, colour = origin)) + 
  geom_violin(position = "dodge", scale = 'width', lwd = 0.3) +
  theme_bw() +
  scale_fill_npg() +
  scale_colour_npg() +
  ylab("AUC") +
  xlab("") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.size = unit(1,"line"),
        legend.text = element_text(size = 6)) -> fig4
pdf(file = file.path(plot_outdir, "fig4.pdf"), width = 4, height = 2)
print(fig4)
dev.off()

## Supplementary Figure 6
all_data %>%
  filter(origin == "sample-out",
         model %in% c("RF", "XGBoost")) %>%
  mutate(sample = gsub("train_ALL_test_sample", "", sample)) %>%
  mutate(idx = as.numeric(gsub("_.*", "", sample))) %>%
  arrange(idx) %>%
  mutate(sample = factor(sample, levels = unique(sample))) %>%
  ggplot(aes(x = sample, y = ROC, fill = model)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = c("#3C5488FF", "#F39B7FFF")) +
  ylab("sample-out AUC (SOE)") +
  xlab("sample") +
  theme(legend.position = "bottom",
        legend.title = element_blank()) -> supp_fig6
pdf(file = file.path(plot_outdir, "supp_fig6.pdf"), width = 9, height = 3)
print(supp_fig6)
dev.off()


## Figure 6

all_data %>%
  group_by(model, sample, origin) %>%
  summarize(med_ROC = median(ROC)) %>%
  ungroup() %>%
  mutate(sample = gsub(".*sample", "", sample)) %>%
  mutate(idx = as.numeric(gsub("_.*", "", sample))) %>%
  arrange(idx) %>%
  mutate(sample = factor(sample, levels = unique(sample))) %>%
  filter(model == "RF") %>%
  ggplot(aes(x = sample, y = med_ROC, fill = origin)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Average RF AUC") +
  xlab("sample") +
  coord_cartesian(ylim = c(0.7, 1)) + 
  scale_fill_manual(values = c("#89CFF0", "#4682B4", "#008081")) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank())-> fig6
pdf(file = file.path(plot_outdir, "fig6.pdf"), width = 9, height = 3)
print(fig6)
dev.off()

## Figure 1: Deamination vs non-deamination VAF

# Load data --

datadir <- "~/data/ENA_SRP044740/tidydata/"
Y.deam.filenames <- list.files(path = datadir, pattern = "_deaminations_Y.rds", full.names = TRUE)
Y.deam <- Y.deam.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.)))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(isDeam == "1")

Y.mut.filenames <- list.files(path = datadir, pattern = "_real-mutations_FFPE_Y.rds", full.names = TRUE)
Y.mut <- Y.mut.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.))))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(isSomatic == "1" | isSNP == "1")

X.deam.filenames <- list.files(path = datadir, pattern = "_deaminations_X.rds", full.names = TRUE)
X.deam <- X.deam.filenames %>%
  set_names(sub("_[^_]+$", "", sub("_[^_]+$", "", basename(.)))) %>%
  map_dfr(readRDS, .id = "sample") %>%
  mutate(complete_id = paste(sample, id, sep = ":")) %>%
  filter(complete_id %in% Y.deam$complete_id)

X.mut.filenames <- list.files(path = datadir, pattern = "_real-mutations_FFPE_X.rds", full.names = TRUE)
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

# Plot
X %>%
  ggplot(aes(x = allele.freq, fill = source)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  xlab("VAF") +
  ylab("") +
  scale_fill_npg(name = "Dose", labels = c("deamination", "non-deamination")) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 8))-> fig1
pdf(file = file.path(plot_outdir, "fig1.pdf"), width = 3, height = 3)
print(fig1)
dev.off()

# Figure 3

AF.cutoff <- 0.3
X %>%
  filter(allele.freq <= AF.cutoff) -> X
Y %>%
  filter(complete_id %in% X$complete_id) -> Y

X %>%
  group_by(sample, source) %>%
  summarize(var_count = n()) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(all_var_count = sum(var_count),
         deam_fraction = var_count/all_var_count) %>%
  ungroup() %>%
  filter(source == "deam") %>%
  arrange(all_var_count) %>%
  mutate(sample = gsub("sample", "", sample)) %>%
  mutate(sample = factor(sample, levels = unique(sample))) %>%
  ggplot(aes(x = sample, y = all_var_count, fill = deam_fraction)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  xlab("sample") +
  ylab("Variant count") +
  scale_fill_gradient(low = "#4DBBD5FF", high = "#E64B35FF") +
  labs(fill = "deamination fraction") +
  theme(legend.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90,  hjust = 1)) -> fig3

pdf(file = file.path(plot_outdir, "fig3.pdf"))
print(fig3)
dev.off()