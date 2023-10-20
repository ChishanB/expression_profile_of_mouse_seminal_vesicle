
# Written: 13/10/23
# Most recent update: 20/10/23

# Purpose: Produce PCA plots from DEseq-analysed counts data. DEseq filters for the top 500 most 
# differentially expressed. This was done because an initial PCA (not included in these scripts)
# containing ALL sncRNAs possessing >=5 counts across >=3 samples showed no distinct separation.

# Author: Chishan Burch
# Contact: chishanburch@gmail.com


###############################################################################


# Load packages
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(DESeq2)
library(stats)
library(ggrepel)
library(ggfortify)

# Read in filtered data. In 1_Extracting_sncRNA_from_excel.R, the raw counts were filtered
# for sncRNAs possessing >=5 counts across >=3 samples. Then these sncRNAs were selected 
# within the TPM dataset. This was an initial quality check to reduce noise.


# Raw
sncRNA_counts <- fread("./1_data/sncRNA/sncRNA_raw_counts_filtered.csv") # Created in 1_Extracting_sncRNA_from_excel.R

sncRNA_counts <- sncRNA_counts %>%
  column_to_rownames(var = "V1")

# TPM
sncRNA_TPM <- fread("./1_data/sncRNA/sncRNA_TPM_counts_filtered.csv") # Created in 1_Extracting_sncRNA_from_excel.R

sncRNA_TPM <- 
  sncRNA_TPM %>%
  data.frame() %>%
  column_to_rownames(var = "V1")


# meta data (samples)
sample_metadata <- fread("./1_data/sncRNA/sncRNA_sample_metadata.csv") %>% 
  as.data.frame() %>%
  mutate(Treatment_group = factor(Treatment_group, levels = c("Control", "Acr")))

# TPM PCA plots

log2_TPM <- log2(sncRNA_TPM + 1)

Gene_Matrix = log2_TPM

# Row normalised by Mean centering (x-RowMean)
Mean_centred <- Gene_Matrix - rowMeans(Gene_Matrix) 

# Row normalised by Z score (x-Rowmean/ SD)
Z_scored <- t(scale(t(Gene_Matrix)))

# Plotting PCAs for TPM

level_colors <- c("Acr" = "#CC79A7", "Control" = "#009E73")
level_colors_darker <- c("Acr" = "#663C53", "Control" = "#004f39")

meancentred_PCA <- prcomp(t(Mean_centred), center = F)

meancentred_PCA_plot <-
  autoplot(meancentred_PCA, 
           data = sample_metadata, 
           colour = "Treatment_group", 
           frame = TRUE, 
           frame.type = 'norm') + scale_fill_manual(values = level_colors) +
  scale_color_manual(values = level_colors_darker) + labs("Treatment group") +
  guides(fill = guide_legend(title = "Treatment group")) + theme_minimal()

#ggsave(plot = meancentred_PCA_plot, filename = "./2_figures/meancentred_TPM_PCA.png", width = 8, height = 6, dpi = 300)

zscored_PCA <- prcomp(t(Z_scored), center = F)

zscored_PCA_plot <-
  autoplot(zscored_PCA, 
           data = sample_metadata, 
           colour = "Treatment_group", 
           frame = TRUE, 
           frame.type = 'norm') + scale_fill_manual(values = level_colors) +
  scale_color_manual(values = level_colors_darker) + labs("Treatment group") +
  guides(fill = guide_legend(title = "Treatment group")) + theme_minimal()

#ggsave(plot = zscored_PCA_plot, filename = "./2_figures/zscored_TPM_PCA.png", width = 8, height = 6, dpi = 300)


# Clear separation was not observed between the two treatment groups for the zscored 
# or mean centered TPM counts, so it was decided that a PCA should be made utilising the 
# top 500 most differentially expressed sncRNAs. These are identified with DE analysis.

# DE analysis needs to be performed on raw counts, and not TPM because DESeq has its own forms of normalisation.


# Perform DESeq normalisation (# DO NOT log2 raw counts before DEseq)

dds <- DESeqDataSetFromMatrix(countData = sncRNA_counts,
                              colData = sample_metadata,
                              desig = ~ Treatment_group)


# Raw counts PCA plot

# Perform variance standardisation with blind = T so the computer is not aware 
# of which treatment group the samples belong to (i.e., unsupervised analysis).

vsd <- vst(dds, blind = T)
vst_matrix <- assay(vsd)

colours <- c("Control" = "#D55E00", "Acr" = "#332288")


PCA_plot <- 
  plotPCA(vsd, intgroup = "Treatment_group", ntop = 500) +
  scale_color_manual(values = colours) +
  stat_ellipse(type = "norm", level = 0.95, alpha = 0.8) + theme_minimal()+
    geom_label_repel(aes(label=colnames(vsd)))

# Save the raw counts PCA plot
#ggsave(plot = PCA_plot, filename = "./2_figures/Principal_component_analysis.png", width = 9, height = 6, dpi = 300)



