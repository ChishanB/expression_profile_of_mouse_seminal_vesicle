
# Written: 5/10/23
# Purpose: pheatmap
# Author: Chishan Burch
# Contact: chishanburch@gmail.com

##################

#install.packages("pheatmap")

# Load packages
library(pheatmap)
library(data.table)
library(tibble)
library(DESeq2)
library(ggplot2)

# Re-create the DEseq object like in 2_PCA_and_Deseq.R

sncRNA_counts <- fread("./1_data/sncRNA/sncRNA_raw_counts_filtered.csv") # Created at the end of 1_Extracting_sncRNA_from_excel.R
#sncRNA_counts <- fread("./1_data/sncRNA/sncRNA_TPM_counts_filtered.csv")
Deseq_results <- fread("./1_data/sncRNA/sncRNA_DeSeq_padj_under_0.05.csv")
gene_metadata <- fread("./1_data/sncRNA/sncRNA_gene_metadata.csv")

sncRNA_counts <- 
  sncRNA_counts %>%
  data.frame() %>%
  column_to_rownames(var = "V1")

sample_metadata <- fread("./1_data/sncRNA/sncRNA_sample_metadata.csv")
sample_metadata <- as.data.frame(sample_metadata) %>% column_to_rownames(., var = "sample.IDs")


dds <- DESeqDataSetFromMatrix(countData = sncRNA_counts, # Needs to be raw counts with rownames set as RNA ids
                              colData = sample_metadata, # Data frame of meta data
                              desig = ~ Treatment_group) 
dds <- DESeq(dds)


# vst normalise the DEseq object
vsd <- vst(dds, blind = T)

# Convert it to a data frame
str(assay(vsd))
Deseq_vst <- assay(vsd)

# Subset for DE genes
Deseq_vector <- Deseq_results$V1

vst_subset <- Deseq_vst[rownames(Deseq_vst) %in% Deseq_vector, ]


# Create Mean centred and z-score (to create two different plots)
Mean_centred <- vst_subset- rowMeans(vst_subset)
Z_scored <- t(scale(t(vst_subset)))

# Use pheatmap to plot, using default hierarchical clustering and using the following variables:

 ## Declare custom colours
custom_colours = list(
  Treatment_group = c(Acr = "#332288", Control = "#D55E00"),
  Type = c(hairpin = "#0072B2", mature = "#CC79A7", piRNA = "#F0E442", precursor = "#882255", 
           rRNA = "#009E73", snoRNA = "#E69F00") 
        )


pheatmap_1 <- pheatmap(vst_subset, scale = "none", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", annotation_col = sample_annotation,
         annotation_row = miRNA_annotation, annotation_colors = custom_colours)

# Create function to save heatmap
  save_pheatmap_png <- function(x, filename, width=980, height=1200, res = 150) {
    png(filename, width = width, height = height, res = res)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
# Save heatmap
  save_pheatmap_png(pheatmap_1, "./2_figures/pheatmap_1.png")
  
# Select the types from gene_metadata, and assign ids as rownames
gene_metadata_df <- 
  gene_metadata %>%
  filter(sRNA.ID_type %in% rownames(vst_subset)) %>%
  data.frame() %>%
  column_to_rownames(var = "sRNA.ID_type") %>%
  dplyr::select(Type)

sample_annotation <- sample_metadata 
miRNA_annotation <- gene_metadata_df


# Mean-centred pheatmap
pheatmap_mean_centred <- pheatmap(Mean_centred, scale = "none", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", annotation_col = sample_annotation,
         annotation_row = miRNA_annotation, annotation_colors = custom_colours)
# Save
save_pheatmap_png(pheatmap_mean_centred, "./2_figures/pheatmap_mean_centred.png")

# Z-scored pheatmap
pheatmap_zscored <- pheatmap(Z_scored, scale = "none", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", annotation_col = sample_annotation,
         annotation_row = miRNA_annotation, annotation_colors = custom_colours)

# Save
save_pheatmap_png(pheatmap_zscored, "./2_figures/pheatmap_zscored.png")


# TPM values pheatmap (ensure the TPM filtered counts csv was used in Deseq etc. 
#Other 3 plots use raw filtered counts csv)

vst_subset <- log2(vst_subset)

pheatmap_TPMcounts <- pheatmap(vst_subset, scale = "none", clustering_distance_rows = "euclidean", 
                                  clustering_distance_cols = "euclidean", annotation_col = sample_annotation,
                                  annotation_row = miRNA_annotation, annotation_colors = custom_colours)
# Save
save_pheatmap_png(pheatmap_TPMcounts, "./2_figures/pheatmap_TPMcounts_log2.png")

           
           

