
# Written: 13/10/23
# Most recent update: 20/10/23

# Purpose: Read in filtered (>=5 counts across >=3 samples) raw and TPM counts.
# Produce PCA plots from log2 transformed filtered TPM counts that have been 
# mean-centered (first PCA) and z-scored (second PCA). 
# Next, produce PCA from the top 500 most differentially expressed sncRNAs as determined
# by performing DE analysis on the raw counts. Perform vst normalisation on the
# dds object before plotting. Lastly, creating 2 csvs containing DEseq results for later
# use. "sncRNA_DeSeq_padj_under_0.05.csv" and "sncRNA_DeSeq_complete.csv"

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
# top 500 most differentially expressed sncRNAs as assessed by DE analysis. ntop = 500 is specified
# when plotting the PCA.

# DE analysis needs to be performed on raw counts, and not TPM because DESeq has its own forms of normalisation.


# Perform DESeq normalisation (# DO NOT log2 raw counts before DEseq)

dds <- DESeqDataSetFromMatrix(countData = sncRNA_counts,
                              colData = sample_metadata,
                              desig = ~ Treatment_group)


# Perform variance standardisation with blind = T so the computer is not aware 
# of which treatment group the samples belong to (i.e., unsupervised analysis).

vsd <- vst(dds, blind = T)
vst_matrix <- assay(vsd)


# Raw counts PCA plot

colours <- c("Control" = "#D55E00", "Acr" = "#332288")

# ntop 500 specifies that we want the top 500 most differentially expressed sncRNAs to form the raw counts PCA
PCA_plot <- 
  plotPCA(vsd, intgroup = "Treatment_group", ntop = 500) +
  scale_color_manual(values = colours) +
  stat_ellipse(type = "norm", level = 0.95, alpha = 0.8) + theme_minimal()+
    geom_label_repel(aes(label=colnames(vsd)))

# Save the raw counts PCA plot
#ggsave(plot = PCA_plot, filename = "./2_figures/Principal_component_analysis.png", width = 9, height = 6, dpi = 300)

# Re-create dds without the vst normalisation (pie chart won't require the vst, and the vst will simply be
# re-done for the pheatmap. Save remade DEseq results for later use in pheatmap, piecharts etc

dds <- DESeqDataSetFromMatrix(countData = sncRNA_counts, # Needs to be raw counts with rownames as RNAs
                              colData = sample_metadata, #data frame of meta data (rownames are samples)
                              desig = ~ Treatment_group) # Column in metadata for treatment


dds <- DESeq(dds)

res <- results(dds) 

res_df <- as.data.frame(res) # Coerce into a more easily accessible format.

# write.csv(res_df, file = "./1_data/sncRNA/sncRNA_DeSeq_complete.csv", row.names = TRUE)

# Filter dds results based on a significance threshold of p-adjusted < 0.05.
res_df_significant <- res_df %>%
  filter(padj <= 0.05)

# Save dds results meeting significance threshold. 
# write.csv(res_df_significant, file = "./1_data/sncRNA/sncRNA_DeSeq_padj_under_0.05.csv", row.names = FALSE)

# How many RNAs upregulated in response to acrylamide exposure (treated)?
# nrow(res_df_significant [res_df_significant $ log2FoldChange > 0, ])
#    12 upregulated RNAs
# How many RNAs downregulated in response to acrylamide exposure (treated)?
# nrow(res_df_significant [res_df_significant $ log2FoldChange < 0, ])
#    46 downregulated RNAs


