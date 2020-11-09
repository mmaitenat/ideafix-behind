#!/usr/bin/env Rscript

library(tidyverse)
library(reshape2)
library(ggsci)
library(ggpubr)

## Some common data ------------------------------
# Change all the filepaths accordingly
plot_outdir <- "results/figs/"
colonliver_samplenames_filename <- "maxwell_FFPE_FF_tumor_pairs.csv"
colonliver_samplenames <- read_csv(colonliver_samplenames_filename, col_names = c("FFPE", "FF", "tissue", "remove1", "remove2")) %>%
  select(FFPE, FF, tissue) %>%
  mutate(sample = paste0(tissue, rep(1:2, times = 2))) %>%
  select(-tissue)

# Median fragment length ---------

calc_median_fraglength <- function(filename) {
  cmd <- sprintf("awk '{ if ($1 > 0) { print } }' %s | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'", filename)
  system(cmd, intern = TRUE) %>% 
    as.numeric()
}

colonliver_fragment_datadir <- "~/data/colon_liver/"
colonliver_fraglength_filenames <- list.files(colonliver_fragment_datadir, pattern = "*fragment_length.txt", recursive = TRUE, full.names = TRUE)
names(colonliver_fraglength_filenames) <- gsub("/", "", gsub("/aligned", "", gsub(colonliver_fragment_datadir, "", gsub("/[^/]+$", "", colonliver_fraglength_filenames))))


colonliver_median_fraglength <- map_dbl(colonliver_fraglength_filenames, calc_median_fraglength) %>%
  tibble(sample = names(.), value = .) %>%
  left_join(dplyr::select(colonliver_samplenames, -FF), by = c("sample" = "FFPE")) %>%
  rowwise() %>%
  mutate(samplename = ifelse(is.na(sample.y), paste(colonliver_samplenames[[match(sample, colonliver_samplenames$FF)[1], "sample"]], "FF", sep = "_"), sample.y)) %>%
  ungroup() %>%
  select(-sample.y)

# Duplication levels -----

colonliver_dedup_level_filename <- "duplication_levels.csv"
colonliver_dedup_level_data <- read_csv(colonliver_dedup_level_filename, col_names = c("sample", "R1", "R2")) %>%
  mutate(R1 = 100 - R1,
         R2 = 100 - R2) %>%
  rowwise() %>% 
  mutate(mean_dedup = mean(c(R1, R2))) %>%
  arrange(desc(mean_dedup)) %>%
  ungroup() %>%
  mutate(sample = factor(sample, levels = unique(sample)))

# GIV score ---------

colonliver_GIV_datadir <- "~/data/colon_liver/"
colonliver_GIV_filenames <- list.files(colonliver_GIV_datadir, pattern = "damage.out", recursive = TRUE, full.names = TRUE)
names(colonliver_GIV_filenames) <- gsub(colonliver_GIV_datadir, "", gsub("/[^/]+$", "", colonliver_GIV_filenames))

map_dfr(colonliver_GIV_filenames, function(x) {read_tsv(x, col_names = c("abs", "type", "experiment", "count", "family", "damage")) %>%
    filter(type == "C_T")}) -> colonliver_GIV_data

# No damage: GIV score = 1
# Some damage: GIV score > 1.5
# Extensive damage: GIV score > 2

thr <- 1.5
colonliver_GIV_data %>%
  select(experiment, damage) %>%
  rename(FFPE = experiment) -> colonliver_GIV_data

colonliver_GIV_data %>%
  left_join(colonliver_samplenames) %>%
  rowwise() %>%
  mutate(sample = ifelse(is.na(FF), paste(colonliver_samplenames[[match(FFPE, colonliver_samplenames$FF)[1], "sample"]], "FF", sep = "_"), sample)) %>%
  ungroup() %>%
  select(sample, damage) %>%
  arrange(damage) -> colonliver_GIV_data

# Median depth ----

depth_filename <- "~/data/colon_liver/colon_liver.depth"
samplenames <- readLines("run_accessions.txt")
depth_data <- read_tsv(depth_filename, col_names = c("chr", "pos", samplenames))

depth_data %>%
  select(-chr, -pos) %>%
  map(function(x) summary(x)) -> median_depth_per_sample

depth_data %>%
  select(-chr, -pos) %>%
  map(function(x) sd(x)) -> sd_per_sample

unlist(sd_per_sample)/unlist(map(median_depth_per_sample, "Mean")) -> coeff_var_per_sample

plyr::ldply(median_depth_per_sample, `[`, "Median") %>%
  left_join(colonliver_samplenames, by = c(".id" = "FFPE")) %>%
  rowwise() %>%
  mutate(sample = ifelse(is.na(FF), paste(colonliver_samplenames[[match(.id, colonliver_samplenames$FF)[1], "sample"]], "FF", sep = "_"), sample)) %>%
  ungroup() %>%
  arrange(desc(Median)) %>%
  mutate(sample = factor(sample, levels = unique(sample)))  -> depth_table

## Calculate numbers for table  ------------------------------

# Median fragment length ----

# FF
median(colonliver_median_fraglength$value[1:4])
max(colonliver_median_fraglength$value[1:4])
min(colonliver_median_fraglength$value[1:4])

#FFPE
median(colonliver_median_fraglength$value[5:8])
max(colonliver_median_fraglength$value[5:8])
min(colonliver_median_fraglength$value[5:8])

# Median depth ----

# FF
median(coeff_var_per_sample[5:8])
max(coeff_var_per_sample[5:8])
min(coeff_var_per_sample[5:8])
# FFPE
median(coeff_var_per_sample[1:4])
max(coeff_var_per_sample[1:4])
min(coeff_var_per_sample[1:4])

# Duplication levels -----

# FF
median(colonliver_dedup_level_data$mean_dedup[c(2, 4, 5, 6)])
max(colonliver_dedup_level_data$mean_dedup[c(2, 4, 5, 6)])
min(colonliver_dedup_level_data$mean_dedup[c(2, 4, 5, 6)])

# FFPE
median(colonliver_dedup_level_data$mean_dedup[c(1, 3, 7, 8)])
max(colonliver_dedup_level_data$mean_dedup[c(1, 3, 7, 8)])
min(colonliver_dedup_level_data$mean_dedup[c(1, 3, 7, 8)])


# GIV score ---------

# FF
median(colonliver_GIV_data$damage[5:8])
max(colonliver_GIV_data$damage[5:8])
min(colonliver_GIV_data$damage[5:8])

# FFPE
median(colonliver_GIV_data$damage[1:4])
max(colonliver_GIV_data$damage[1:4])
min(colonliver_GIV_data$damage[1:4])


## Plots  ------------------------------

# Prepare data

# Bind all colonliver matrices
depth_table %>%
  left_join(colonliver_median_fraglength, by = c(".id" = "sample")) %>%
  rename(median_fraglength = value,
         median_depth = Median) %>%
  left_join(colonliver_dedup_level_data, by = c(".id" = "sample")) %>%
  rename(mean_dup = mean_dedup) %>%
  left_join(colonliver_GIV_data, by = c("samplename" = "sample")) %>%
  select(-sample, -R1, -R2) -> all_colonliver_damage_data

all_colonliver_damage_data %>%
  mutate(samplename = factor(samplename, levels = unique(samplename))) -> all_colonliver_damage_data

all_colonliver_damage_data %>%
  mutate(source = c(rep("FF", 4),
                    rep("FFPE", 4))) -> all_colonliver_damage_data

# Damage composite plot

# First, define a color palette. We'll do it by giving different transparency levels to an existing palette
mypal_colonliver <- pal_npg("nrc", alpha = 1)(8)

# Down: depth of coverage

depth_boxplot_table <- bind_rows(median_depth_per_sample) %>%
  as.matrix() %>%
  t() 
depth_boxplot_table <- as_tibble(depth_boxplot_table, rownames = "sample") 
colnames(depth_boxplot_table) <- c("sample", "min", "low", "mid", "mean", "top", "max")

depth_boxplot_table <- depth_boxplot_table %>%
  left_join(colonliver_median_fraglength) %>%
  select(-value) %>%
  mutate(samplename = factor(samplename, levels = all_colonliver_damage_data$samplename))

depth_boxplot_table %>%
  ggplot(aes(x = samplename, ymin = min, lower = low, middle = mid, upper = top, ymax = max, fill = samplename)) +
  geom_boxplot(stat = "identity") +
  theme_bw() + 
  ylab("Depth of coverage") +
  xlab("") +
  scale_fill_manual(values = mypal_colonliver) +
  coord_cartesian(ylim = c(0, 400)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")  -> colonliver_boxplot_depth

# Above: fragment length
all_colonliver_damage_data %>%
  ggplot(aes(x = samplename, y = median_fraglength, fill = samplename)) +
  geom_bar(stat = "identity") +
  theme_bw() + 
  xlab("") +
  ylab("Median fragment length") +
  scale_fill_manual(values = mypal_colonliver) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none")  -> colonliver_barplot_fraglength

# Above: duplication level
all_colonliver_damage_data %>%
  ggplot(aes(x = samplename, y = mean_dup, fill = samplename)) +
  geom_bar(stat = "identity") +
  theme_bw() + 
  xlab("") +
  ylab("Mean library duplication level") +
  scale_fill_manual(values = mypal_colonliver) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none")  -> colonliver_barplot_dup

# Top: GIV score
all_colonliver_damage_data %>%
  ggplot(aes(x = samplename, y = damage, colour = samplename)) +
  geom_point(stat = "identity") +
  theme_bw() + 
  xlab("") +
  ylab("GIV score") +
  geom_hline(yintercept = 2, color = "azure4") +
  scale_colour_manual(values = mypal_colonliver) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none")  -> colonliver_giv_dotplot


pdf(file = file.path(plot_outdir, "supp_fig4.pdf"), width = 9, height = 3)
print(colonliver_giv_dotplot)
print(colonliver_barplot_fraglength)
print(colonliver_barplot_dup)
print(colonliver_boxplot_depth)
dev.off()


# Make colonliver deamination vs non deamination vaf plot ----

test_datadir <-  "~/data/colon_liver/tidydata"
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

X %>%
  ggplot(aes(x = allele.freq, fill = source)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  xlab("VAF") +
  ylab("") +
  scale_fill_npg(name = "Dose", labels = c("deamination", "non-deamination")) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 8))-> fig_vaf_colonliver_a

X %>%
  filter(allele.freq <= 0.3) %>%
  ggplot(aes(x = allele.freq, fill = source)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  xlab("VAF") +
  ylab("") +
  scale_fill_npg(name = "Dose", labels = c("deamination", "non-deamination")) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 8))-> fig_vaf_colonliver_b

vaf_figure_final <- ggarrange(fig_vaf_colonliver_a, fig_vaf_colonliver_b,
                              labels = c("A", "B"),
                              ncol = 2, nrow = 1, 
                              common.legend = TRUE, 
                              legend = "bottom")

pdf(file = file.path(plot_outdir, "supp_fig2.pdf"), width = 6, height = 3)
print(vaf_figure_final)
dev.off()