# Written: 28/9/23
# Updated: 5/10/23
# Most recent update: 
# Purpose: Target mRNA prediction for mature miRNAs meeting the padj < 0.05 threshold.  
# Author: Chishan Burch
# Contact: chishanburch@gmail.com

#######################################################################################

# install multiMiR package for target micro RNA prediction

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("multiMiR")
#BiocManager::install("edgeR")

# Also installed dependency "limma"

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("biomaRt")

# View documentation for the version of this package installed 
#browseVignettes("multiMiR")
#browseVignettes("biomaRt")
#install.packages("dplyr")
# Load libraries
#library(openxlsx)
library(data.table)
library(dplyr)
#library(tibble)

library(ggplot2)
library(ggvenn)

library(multiMiR)
#library(edgeR)
#library(limma)


# Take the DeSeq data generated in the previous script, and filter by thresholds.
#miRNA_df <- fread("./1_data/sncRNA/sncRNA_DeSeq_object.csv")
miRNA_df <- fread("./1_data/sncRNA/sncRNA_DeSeq_padj_under_0.05.csv")
miRNA_df_all <- fread("./1_data/sncRNA/sncRNA_DeSeq_complete.csv") #this one will be used for background_list in 7_GO_target_prediction.R

# The threshold this time is UPREGULATED miRNAs (positive log2foldchange) and padj < 0.05 (already applied)
miRNA_UP <- miRNA_df %>%
filter(log2FoldChange > 0)

# Set the column name to miRNA_ids
miRNA_UP <- setnames(miRNA_UP, "V1", "miRNA_ids")


miRNA_all <- setnames(miRNA_df_all, "miRNA_ids_underscores", "miRNA_ids", skip_absent = TRUE) #For GO prediction

miRNA_UP <-
  miRNA_UP %>%
  mutate(
    Type = case_when(
      grepl(pattern = "hairpin", x = miRNA_ids) ~ "hairpin",
      grepl(pattern = "mature", x = miRNA_ids) ~ "mature"
    )
  )

miRNA_all <-
  miRNA_all %>%
  mutate(
    Type = case_when(
      grepl(pattern = "hairpin", x = miRNA_ids) ~ "hairpin",
      grepl(pattern = "mature", x = miRNA_ids) ~ "mature"
    )
  )

# Select only rows containing mature miRNAs
mature_miRNAs <- miRNA_UP[Type == "mature"]
mature_miRNAs_all <- miRNA_all[Type == "mature"]

#Select the miRNA_ids and remove the Type ("_mature", "_hairpin" etc so the ID can be used for searching)
mature_miRNA_ids <- mature_miRNAs %>%
  dplyr::select(miRNA_ids)%>%
  .[["miRNA_ids"]] %>%
  gsub("_mature", "", .)
str()


mature_miRNA_ids_2 <- mature_miRNAs_all %>%
  dplyr::select(miRNA_ids)%>%
  .[["miRNA_ids"]] %>%
  gsub("_mature", "", .)
str()

# need to specify org = "mmu" every time the function is used
multimiR_results <- get_multimir(org = "mmu",
                                 mirna = mature_miRNA_ids,
                                 table = c("predicted"),
                                 summary = TRUE)
# Validated mRNA targets
multimiR_results_validated <- get_multimir(org = "mmu",
                                 mirna = mature_miRNA_ids,
                                 table = c("validated"),
                                 summary = TRUE)

multimiR_results_all <- get_multimir(org = "mmu",
                                 mirna = mature_miRNA_ids_2,
                                 table = c("predicted"),
                                 summary = TRUE)
# Validated mRNA targets
multimiR_results_validated_all <- get_multimir(org = "mmu",
                                           mirna = mature_miRNA_ids_2,
                                           table = c("validated"),
                                           summary = TRUE)


# Pass multimiR results as a data frame.
# multimiR_resuslts contains the predicted mRNA targets for upregulated miRNAs
multimiR_results <- as.data.frame(multimiR_results@data) 
multimiR_results_validated <- as.data.frame(multimiR_results_validated@data) #2

multimiR_results_all <- as.data.frame(multimiR_results_all@data) 
multimiR_results_validated_all <- as.data.frame(multimiR_results_validated_all@data) 

# Save
# Are the content of the top 2 even right anymore?
#write.csv(multimiR_results, file = "./1_data/mRNA/multimiR_miRNA_UP.csv")
#write.csv(multimiR_results_validated, file = "./1_data/mRNA/multimiR_miRNA_UP_validated.csv")

#write.csv(multimiR_results_all, file = "./1_data/mRNA/multimiR_miRNA_UP_all.csv")
#write.csv(multimiR_results_validated_all, file = "./1_data/mRNA/multimiR_miRNA_UP_validated_all.csv")

#####

# Read in the dataset containing differentially regulated mRNAs (experimental data). 
# We want to find any overlap between UPREGULATED miRNAs and DOWNREGULATED mRNAs
mRNA_df <- fread("./1_data/mRNA/GSE183632_SkerrettByrne_DESeq_list_Geo_dataset.csv")

# Filter based on thresholds which represent the overlap we are seeking. 
#I.e., padj < 0 , log2foldchange < 0
mRNA_DOWN <- mRNA_df %>%
  filter(log2FoldChange < 0 & Padj < 0.05)
# Save
#write.csv(mRNA_DOWN, file = "./1_data/mRNA/mRNA_down_SkerrettByrne.csv", row.names = FALSE)


# Create venn diagram showing overlap of validated mature miRNA targets (UPREGULATED in dataset) 
# and mRNAs DOWNREGULATED in the dataset due to acrylamide exposure

# Remove those with blanks in target_entrez, as this is the id we are choosing to use for
# consistency
multimiR_results %>% filter(target_entrez != "")
multimiR_results_validated %>% filter(target_entrez != "")

miRNA_targets <- dplyr::select(multimiR_results, target_entrez) %>%
  .[["target_entrez"]]
mRNA_DOWN <- dplyr::select(mRNA_DOWN, target_entrez) %>%
  .[["target_entrez"]]

miRNA_targets_validated <- dplyr::select(multimiR_results_validated, target_entrez) %>%
  .[["target_entrez"]] #2



# sum() the number of variables in mRNA_DOWN matching to variables in miRNA_targets. 
# Because TRUE = 1, sum(dataframe1 %in% dataframe2) will sum the number of times values in 
# dataframe1 correspond to dataframe2.

sum(mRNA_DOWN %in% miRNA_targets)
sum(mRNA_DOWN %in% miRNA_targets_validated) #2

ggven_list <- list(predicted_miRNA_targets = as.character(miRNA_targets),
                   mRNA = as.character(mRNA_DOWN))

ggven_list_validated <- list(validated_miRNA_targets = as.character(miRNA_targets_validated),
                   mRNA = as.character(mRNA_DOWN))

ggvenn(ggven_list) + scale_fill_manual(values = c("#AA4499", "#66CCEE"))
#ggsave(filename = "./2_figures/miRNA_up_mRNA_down_pre.png", width = 6, height = 4)

ggvenn(ggven_list_validated) + scale_fill_manual(values = c("#332288", "#66CCEE"))
#ggsave(filename = "./2_figures/miRNA_up_mRNA_down_val.png", width = 6, height = 4)

# Get the ids of mRNAs in the intersection of predicted miRNA targets and mRNA downregulated 
intersecting_mRNAs_v1 <- miRNA_targets %in% mRNA_DOWN
sum(intersecting_mRNAs_v1)
intersecting_mRNAs_v1 <- miRNA_targets[intersecting_mRNAs_v1] # Result

# Get their names using the target_entrez ids
mRNA_df_2 <- mRNA_df %>%
  filter(log2FoldChange < 0)

Symbol <- mRNA_df_2$Symbol # These 2 lines will be reused for the validated targets
target_entrez <- mRNA_df_2$target_entrez #

mRNA_ids <- data.frame(Symbol, target_entrez)

mRNA_ids <- mRNA_ids %>%
  filter(target_entrez %in% intersecting_mRNAs_v1)
# Save ids
write.csv(mRNA_ids, file = "./1_data/mRNA/predictednames_miUP.csv", row.names = FALSE)

# Do the same for the validated targets
intersecting_mRNAs_v2 <- miRNA_targets_validated %in% mRNA_DOWN
intersecting_mRNAs_v2 <- miRNA_targets_validated[intersecting_mRNAs_v2]

mRNA_ids_validated <- data.frame(Symbol, target_entrez)
mRNA_ids_validated <- mRNA_ids_validated %>%
  filter(target_entrez %in% intersecting_mRNAs_v2)
# Save ids
write.csv(mRNA_ids_validated, file = "./1_data/mRNA/validatednames_miUP.csv", row.names = FALSE)




# # # # # #



combined_targets <- c(miRNA_targets_validated, miRNA_targets)
combined_targets <- unique(combined_targets)

mRNA_DOWN

ggven_list_combined <- list(all_miRNA_targets = as.character(combined_targets),
                          DOWNREG_mRNA = as.character(mRNA_DOWN))


ggvenn(ggven_list_combined) + scale_fill_manual(values = c("#EE6677", "#44BB99"))
ggsave(filename = "./2_figures/miRNA_up_mRNA_down_all.png", width = 6, height = 4)









                                      
