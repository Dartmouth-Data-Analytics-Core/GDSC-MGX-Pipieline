#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# Title: plot_alpha.R
# Description: Generate alpha diversity plots
#
# Author: Mike Martinez
# Lab: GDSC
# Project: Clover-Seq
# Date created: 05/16/25
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD LIBRARIES AND SET PATHS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(dplyr)
library(vegan)
library(ggplot2)
library(data.table)
library(tidyr)
library(ggrepel)

#----- Specify outDir path
resultsDir <- "results/"

#----- Set command line args
args <- commandArgs(trailingOnly = TRUE)

#----- Check that all arguments are supplied
if (length(args) < 1 | length(args) > 1) {
  stop("Usage: RScript plot_alpha.R <sample_sheet.csv>")
}

#----- Set variables based on command line args
sampleSheet <- args[1]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN THE DATA, TIDY, AND ARRANGE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Read in the sample sheet
sampleList <- fread(sampleSheet) |>
    as.data.frame()

#----- Get vector of samples and form column names
samples <- sampleList[,1]
specCols <- c(samples)
genusCols <- c(samples)

#----- Read in the species only data
spec <- fread(paste0(resultsDir, "merged_species_abundance_table.txt")) |>
    as.data.frame()
rownames(spec) <- spec[,1]
spec[,1] <- NULL
colnames(spec) <- specCols
fwrite(spec, file = paste0(resultsDir, "named_merged_species_abundance_table.txt"))

#----- Transpose the data (taxa as columns)
specT <- as.data.frame(t(spec))

#----- Read in the genus only data
genus <- fread(paste0(resultsDir, "merged_genus_abundance_table.txt")) |>
    as.data.frame()
rownames(genus) <- genus[,1]
genus[,1] <- NULL
colnames(genus) <- genusCols
fwrite(genus, file = paste0(resultsDir, "named_merged_genus_abundance_table.txt"))

#----- Transpose the data (taxa as columns)
genusT <- as.data.frame(t(genus))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# CALCULATE ALPHA DIVERSITY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- On species
specShan <- vegan::diversity(specT, index = "shannon")
specSimp <- vegan::diversity(specT, index = "simpson")
specObs <- rowSums(specT > 0)

#----- Combine into dataframe
specAlpha <- data.frame(
    Sample = rownames(specT),
    Shannon = specShan,
    Simpson = specSimp,
    Observed = specObs
)

#----- On Genus
genShan <- vegan::diversity(genusT, index = "shannon")
genSimp <- vegan::diversity(genusT, index = "simpson")
genObs <- rowSums(genusT > 0)

#----- Combine into dataframe
genusAlpha <- data.frame(
    Sample = rownames(genusT),
    Shannon = genShan,
    Simpson = genSimp,
    Observed = genObs
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PLOT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Pivot longer for plotting (species)
specAlphaLong <- specAlpha |>
  pivot_longer(cols = c(Shannon, Simpson, Observed),
               names_to = "Metric",
               values_to = "Value") |>
    as.data.frame()
fwrite(specAlphaLong, file = paste0(resultsDir, "Species_alpha_diversity_for_plotting.csv"))

#----- Pivot longer for plotting (genus)
genusAlphaLong <- genusAlpha |>
  pivot_longer(cols = c(Shannon, Simpson, Observed),
               names_to = "Metric",
               values_to = "Value") |>
    as.data.frame()
fwrite(genusAlphaLong, file = paste0(resultsDir, "Genus_alpha_diversity_for_plotting.csv"))

#----- Plot species shannon
x <- ggplot(specAlpha, aes(x = Sample, y = Shannon)) +
  geom_col(color = "black", width = 0.7, alpha = 0.8) +
  geom_text(aes(label = round(Shannon, 2)), vjust = -0.4, size = 3.5) +
  theme_classic(base_size = 15) +
  labs(
    title = "Shannon Alpha Diversity per Sample",
    x = "",
    y = "Shannon Index"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  )
ggsave(paste0(resultsDir, "Species_Shannon_Diversity_Per_Sample.png"), x)

#----- Plot Genus Shannon
y <- ggplot(genusAlpha, aes(x = Sample, y = Shannon)) +
    geom_col(color = "black", width = 0.7, alpha = 0.8) +
    geom_text(aes(label = round(Shannon, 2)), vjust = -0.4, size = 3.5) +
    theme_classic(base_size = 15) +
    labs(
        title = "Shannon Alpha Diversity per Sample",
        x = "",
        y = "Shannon Index"
    ) +
    scale_fill_brewer(palette = "Set2") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.title = element_blank()
    )
ggsave(paste0(resultsDir, "Genus_Shannon_Diversity_Per_Sample.png"), x)