#!/usr/bin/env Rscript

library(tidyverse)
library(reshape2)
library(ggsci)
library(ggpubr)

## Some common data ------------------------------
# Change all the filepaths accordingly
samplenames_filename <- "sample_pairs_with_names.csv"
samplenames <- read_csv(samplenames_filename, col_names = c("FF", "FFPE", "sample"))
plot_outdir <- "results/figs/"

## Load the data and calculate values  ------------------------------

# Median fragment length ----

calc_median_fraglength <- function(filename) {
  cmd <- sprintf("awk '{ if ($1 > 0) { print } }' %s | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'", filename)
  system(cmd, intern = TRUE) %>% 
    as.numeric()
}

fragment_datadir <- "~/data/ENA_SRP044740/"
fraglength_filenames <- list.files(fragment_datadir, pattern = "*fragment_length.txt", recursive = TRUE, full.names = TRUE)
names(fraglength_filenames) <- gsub("/", "", gsub("/aligned", "", gsub(fragment_datadir, "", gsub("/[^/]+$", "", fraglength_filenames))))

median_fraglength <- map_dbl(fraglength_filenames, calc_median_fraglength) %>%
  tibble(sample = names(.), value = .) %>%
  left_join(dplyr::select(samplenames, -FF), by = c("sample" = "FFPE")) %>%
  rowwise() %>%
  mutate(samplename = ifelse(is.na(sample.y), paste(samplenames[[match(sample, samplenames$FF)[1], "sample"]], "FF", sep = "_"), sample.y))

# Median depth ----

depth_filename <- "~/data/ENA_SRP044740/ENA.depth"
sample_suffix <- 30:69
samples <- paste("SRR15232", sample_suffix, sep = "")
depth_data <- read_tsv(depth_filename, col_names = c("chr", "pos", samples))

depth_data %>%
  select(-chr, -pos) %>%
  map(function(x) summary(x)) -> median_depth_per_sample

plyr::ldply(median_depth_per_sample, `[`, "Median") %>%
  left_join(samplenames, by = c(".id" = "FFPE")) %>%
  rowwise() %>%
  mutate(sample = ifelse(is.na(FF), paste(samplenames[[match(.id, samplenames$FF)[1], "sample"]], "FF", sep = "_"), sample)) -> depth_table

depth_data %>%
  select(-chr, -pos) %>%
  map(function(x) sd(x)) -> sd_per_sample

unlist(sd_per_sample)/unlist(map(median_depth_per_sample, "Mean")) -> coeff_var_per_sample

# Duplication levels -----

FFPE_dedup_level_filename <- "FFPE_duplication_levels.csv"
FFPE_dedup_level_data <- read_delim(FFPE_dedup_level_filename, col_names = c("sample", "R1", "R2"), delim = ";") %>%
  mutate(R1 = 100 - R1,
         R2 = 100 - R2) %>%
  rowwise() %>% 
  mutate(mean_dedup = mean(c(R1, R2))) %>%
  arrange(desc(mean_dedup)) %>%
  mutate(sample = factor(sample, levels = unique(sample))) 

FF_dedup_level_filename <- "FF_duplication_levels.csv"
FF_dedup_level_data <- read_csv(FF_dedup_level_filename, col_names = c("sample", "R1", "R2")) %>%
  mutate(R1 = 100 - R1,
         R2 = 100 - R2) %>%
  rowwise() %>% 
  mutate(mean_dedup = mean(c(R1, R2))) %>%
  arrange(desc(mean_dedup)) %>%
  mutate(sample = factor(sample, levels = unique(sample))) 

FF_dedup_level_data %>%
  rename(sample0 = sample) %>%
  left_join(samplenames, by = c("sample0" = "FFPE"))  %>%
  rowwise() %>%
  mutate(sample = paste(samplenames[[match(sample0, samplenames$FF)[1], "sample"]], "FF", sep = "_")) %>%
  select(sample, R1, R2, mean_dedup) -> FF_dedup_level_data

all_dedup_level_data <- bind_rows(FF_dedup_level_data,
                                  FFPE_dedup_level_data) %>%
  ungroup()

# GIV score ---------

GIV_datadir <- "~/data/ENA_SRP044740"
filenames <- list.files(GIV_datadir, pattern = "damage.out", recursive = TRUE, full.names = TRUE)
names(filenames) <- gsub(GIV_datadir, "", gsub("/[^/]+$", "", filenames))

map_dfr(filenames, function(x) {read_tsv(x, col_names = c("abs", "type", "experiment", "count", "family", "damage")) %>%
    filter(type == "C_T")}
) -> data

# No damage: GIV score = 1
# Some damage: GIV score > 1.5
# Extensive damage: GIV score > 2

thr <- 1.5
data %>%
  select(experiment, damage) %>%
  rename(FFPE = experiment) -> giv_data

giv_data %>%
  left_join(samplenames) %>%
  rowwise() %>%
  mutate(sample = ifelse(is.na(FF), paste(samplenames[[match(FFPE, samplenames$FF)[1], "sample"]], "FF", sep = "_"), sample)) %>%
  ungroup() %>%
  select(sample, damage) %>%
  arrange(damage) -> giv_data

## Calculate numbers for table  ------------------------------

# Median fragment length ----

# FF
median(median_fraglength$value[1:13])
max(median_fraglength$value[1:13])
min(median_fraglength$value[1:13])

#FFPE
median(median_fraglength$value[14:40])
max(median_fraglength$value[14:40])
min(median_fraglength$value[14:40])

# Median depth ----

# FF
median(coeff_var_per_sample[1:13])
max(coeff_var_per_sample[1:13])
min(coeff_var_per_sample[1:13])
# FFPE
median(coeff_var_per_sample[14:40])
max(coeff_var_per_sample[14:40])
min(coeff_var_per_sample[14:40])

# Duplication levels -----

# FF
median(all_dedup_level_data$mean_dedup[1:13])
max(all_dedup_level_data$mean_dedup[1:13])
min(all_dedup_level_data$mean_dedup[1:13])

# FFPE
median(all_dedup_level_data$mean_dedup[14:40])
max(all_dedup_level_data$mean_dedup[14:40])
min(all_dedup_level_data$mean_dedup[14:40])


# GIV score ---------

# FF
median(giv_data$damage[1:13])
max(giv_data$damage[1:13])
min(giv_data$damage[1:13])

# FFPE
median(giv_data$damage[14:40])
max(giv_data$damage[14:40])
min(giv_data$damage[14:40])


## Plots  ------------------------------

# Prepare data

depth_table %>%
  arrange(desc(Median)) %>%
  mutate(sample = factor(sample, levels = unique(sample))) %>%
  select(sample, Median) %>%
  ungroup() -> median_depth

median_fraglength %>%
  select(samplename, value) %>%
  rename(median_fraglength = value) %>%
  left_join(median_depth, by = c("samplename" = "sample")) %>%
  rename(median_depth = Median) %>%
  ungroup() -> all_damage_data0

all_damage_data0 %>%
  left_join(all_dedup_level_data, by = c("samplename" = "sample")) %>%
  select(-R1, -R2) -> all_damage_data1

all_damage_data1 %>%
  left_join(giv_data, by = c("samplename" = "sample")) %>%
  rowwise() %>%
  mutate(source = ifelse(grepl("FF", samplename), "FF", "FFPE")) %>%
  ungroup() %>%
  mutate(samplename = factor(samplename, levels = unique(samplename))) -> all_damage_data

# Damage composite plot

# First, define a color palette. We'll do it by giving different transparency levels to an existing palette
mypal <- c(pal_npg("nrc", alpha = 0.2)(10),
           pal_npg("nrc", alpha = 0.5)(10),
           pal_npg("nrc", alpha = 0.7)(10),
           pal_npg("nrc", alpha = 0.9)(10))

# Down: depth of coverage
depth_boxplot_table <- bind_rows(median_depth_per_sample) %>%
  as.matrix() %>%
  t()
depth_boxplot_table <- as_tibble(depth_boxplot_table, rownames = "sample") 
colnames(depth_boxplot_table) <- c("sample", "min", "low", "mid", "mean", "top", "max")

# Add names
depth_boxplot_table %>%
  left_join(samplenames, by = c("sample" = "FFPE")) %>%
  rowwise() %>%
  mutate(sample.y = ifelse(is.na(sample.y), paste(samplenames[[match(sample, samplenames$FF)[1], "sample"]], "FF", sep = "_"), sample.y)) %>%
  ungroup() %>%
  select(-sample, -FF) %>%
  rename(sample = sample.y) %>%
  mutate(sample = factor(sample, levels = unique(sample))) -> depth_boxplot_table

depth_boxplot_table %>%
  ggplot(aes(x = sample, ymin = min, lower = low, middle = mid, upper = top, ymax = max, fill = sample)) +
  geom_boxplot(stat = "identity") +
  theme_bw() + 
  ylab("Depth of coverage") +
  xlab("") +
  scale_fill_manual(values = mypal) +
  coord_cartesian(ylim = c(0, 400)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")  -> boxplot_depth

# Above: fragment length
all_damage_data %>%
  ggplot(aes(x = samplename, y = median_fraglength, fill = samplename)) +
  geom_bar(stat = "identity") +
  theme_bw() + 
  xlab("") +
  ylab("Median fragment length") +
  scale_fill_manual(values = mypal) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none")  -> barplot_fraglength

# Above: duplication level
all_damage_data %>%
  ggplot(aes(x = samplename, y = mean_dedup, fill = samplename)) +
  geom_bar(stat = "identity") +
  theme_bw() + 
  xlab("") +
  ylab("Mean library duplication level") +
  scale_fill_manual(values = mypal) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none")  -> barplot_dedup

# Top: GIV score
all_damage_data %>%
  ggplot(aes(x = samplename, y = damage, colour = samplename)) +
  geom_point(stat = "identity") +
  theme_bw() + 
  xlab("") +
  ylab("GIV score") +
  geom_hline(yintercept = 2, color = "azure4") +
  scale_colour_manual(values = mypal) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none")  -> giv_dotplot


pdf(file = file.path(plot_outdir, "fig2.pdf"), width = 9, height = 3)
print(giv_dotplot)
print(barplot_fraglength)
print(barplot_dedup)
print(boxplot_depth)
dev.off()
