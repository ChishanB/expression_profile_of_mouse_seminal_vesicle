# Written: 5/10/23
# Most recent update: 
# Purpose: Part 2 of miR_target_prediction. Now creating venn for downregulated miRNAs 
# Author: Chishan Burch
# Contact: chishanburch@gmail.com

#######################################################################################

# Load libraries
library(data.table)
library(dplyr)

library(ggplot2)
library(ggvenn)

library(multiMiR)
#library(edgeR)
#library(limma)

miRNA_df <- fread("./1_data/sncRNA/sncRNA_DeSeq_padj_under_0.05.csv")

# The threshold this time is DOWNREGULATED miRNAs (negative log2foldchange) and padj < 0.05 (already applied)
miRNA_DOWN <- miRNA_df %>%
  filter(log2FoldChange < 0)

# Set the column name to miRNA ids

miRNA_DOWN <- setnames(miRNA_DOWN, "V1", "miRNA_ids")

miRNA_DOWN <-
  miRNA_DOWN %>%
  mutate(
    Type = case_when(
      grepl(pattern = "hairpin", x = miRNA_ids) ~ "hairpin",
      grepl(pattern = "mature", x = miRNA_ids) ~ "mature"
    )
  )

# Select only rows containing mature miRNAs
mature_miRNAs <- miRNA_DOWN[Type == "mature"]

#Now remove "_type" from each of the miRNA_ids so the ids can be used for searching
mature_miRNA_ids <- mature_miRNAs %>%
  dplyr::select(miRNA_ids)%>%
  .[["miRNA_ids"]] %>%
  gsub("_mature", "", .)

# multimiR analysis
multimiR_results <- get_multimir(org = "mmu",
                                 mirna = mature_miRNA_ids,
                                 table = c("predicted"),
                                 summary = TRUE)

# Validated mRNA targets
multimiR_results_validated <- get_multimir(org = "mmu",
                                           mirna = mature_miRNA_ids,
                                           table = c("validated"),
                                           summary = TRUE)

# Coerce predicted and validated mRNA target multimiR results into data frames
multimiR_results <- as.data.frame(multimiR_results@data)
multimiR_results_validated <- as.data.frame(multimiR_results_validated@data) #2

# Save
write.csv(multimiR_results, file = "./1_data/multimiR_miRNA_DOWN.csv")

# Read in the dataset containing differentially regulated mRNAs (experimental data). 
# We want to find any overlap between DOWNREGULATED miRNAs and UPREGULATED mRNAs
mRNA_df <- fread("./1_data/mRNA/GSE183632_SkerrettByrne_DESeq_list_Geo_dataset.csv")

# Filter based on thresholds which represent the overlap we are seeking. 
#I.e., padj < 0 , log2foldchange < 0
mRNA_UP <- mRNA_df %>%
  filter(log2FoldChange > 0 & Padj < 0.05)
# Save
write.csv(mRNA_UP, file = "./1_data/mRNA/mRNA_up_SkerrettByrne.csv", row.names = FALSE)

# Create venn diagram showing overlap of validated mature miRNA targets (DOWNREGULATED in dataset)
# and mRNAs UPREGULATED in the dataset due to acrylamide exposure 

# Remove any with blanks in target_entrez, as this is the id we are choosing to use for
# consistency
multimiR_results %>% filter(target_entrez != "")

miRNA_targets <- dplyr::select(multimiR_results, target_entrez) %>%
  .[["target_entrez"]]
mRNA_UP <- dplyr::select(mRNA_UP, target_entrez) %>%
  .[["target_entrez"]] 

miRNA_targets_validated <- dplyr::select(multimiR_results_validated, target_entrez) %>%
  .[["target_entrez"]] #2

sum(mRNA_UP %in% miRNA_targets)
sum(mRNA_UP %in% miRNA_targets_validated) #2

ggven_list <- list(predicted_miRNA_targets = as.character(miRNA_targets),
                   mRNA = as.character(mRNA_UP))

ggven_list_validated <- list(validated_miRNA_targets = as.character(miRNA_targets_validated),
                             mRNA = as.character(mRNA_UP))

ggvenn(ggven_list) + scale_fill_manual(values = c("#D55E00", "#44BB99"))
ggsave(filename = "./2_figures/miRNA_down_mRNA_up_pre.png", width = 6, height = 4)

ggvenn(ggven_list_validated) + scale_fill_manual(values = c("#F0E442", "#44BB99"))
ggsave(filename = "./2_figures/miRNA_down_mRNA_up_val.png", width = 6, height = 4)

# Get the ids of mRNAs in the intersection of predicted miRNA targets and mRNA upregulated
# predicted
intersecting_mRNAs_v1 <- miRNA_targets %in% mRNA_UP
sum(intersecting_mRNAs_v1)
intersecting_mRNAs_v1 <- miRNA_targets[intersecting_mRNAs_v1] # Result
# validated
intersecting_mRNAs_v2 <- miRNA_targets_validated %in% mRNA_UP
sum(intersecting_mRNAs_v2)
intersecting_mRNAs_v2 <- miRNA_targets_validated[intersecting_mRNAs_v2] # Result

# Get their names using the target_entrez ids
mRNA_df_2 <- mRNA_df %>%
  filter(log2FoldChange > 0)

Symbol <- mRNA_df_2$Symbol 
target_entrez <- mRNA_df_2$target_entrez

mRNA_ids <- data.frame(Symbol, target_entrez)

mRNA_ids <- mRNA_ids %>%
  filter(target_entrez %in% intersecting_mRNAs_v1)

# Now the validated targets
mRNA_ids_validated <- data.frame(Symbol, target_entrez)

mRNA_ids_validated <- mRNA_ids_2 %>%
  filter(target_entrez %in% intersecting_mRNAs_v2)

# Save ids
write.csv(mRNA_ids, file = "./1_data/mRNA/predictednames_miDOWN.csv", row.names = FALSE)
write.csv(mRNA_ids_validated, file = "./1_data/mRNA/validatednames_miDOWN.csv", row.names = FALSE)

predicted_names <- fread("./1_data/mRNA/predictednames_miDOWN.csv")
validated_names <- fread("./1_data/mRNA/validatednames_miDOWN.csv")
# # # # # # # # # #



combined_targets <- c(miRNA_targets_validated, miRNA_targets)
combined_targets <- unique(combined_targets)

mRNA_UP

ggven_list_combined <- list(all_miRNA_targets = as.character(combined_targets),
                   UPREG_mRNA = as.character(mRNA_UP))


ggvenn(ggven_list_combined) + scale_fill_manual(values = c("#AA4499", "#44BB99"))
ggsave(filename = "./2_figures/miRNA_down_mRNA_up_all.png", width = 6, height = 4)






#From Paul Tol: https://personal.sron.nl/~pault/
Tol_bright <- c('#EE6677', '#228833', '#4477AA', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB')
Tol_muted <- c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', '#882255', '#AA4499', '#DDDDDD')
Tol_light <- c('#BBCC33', '#AAAA00', '#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF', '#44BB99', '#DDDDDD')
#From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
