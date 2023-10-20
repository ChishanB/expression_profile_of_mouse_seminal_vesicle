
# Written: 31/8/23
# Most recent update: 20/10/23

# Purpose: Extracting the sheet from the dataset containing the 
# sncRNA RNA-Seq dataset, and splitting the columns into separate csv files
# for use in later scripts. Then, performing some initial quality checks using bar graphs.
# I.e., getting a 'first look' at the data before performing PCA.

# ****NOTE: Mixomics has not been utilised on this data or within this project. Mixomics 
# analysis is one of the future directions of this project.

# Author: Chishan Burch
# Contact: chishanburch@gmail.com

#######################################################################################

# Load Packages: 
library(openxlsx)
library(dplyr)
library(data.table)

# Read in the Excel file and specify the sheet desired. Store it in the variable "sncRNA"
sncRNA <- read.xlsx("./1_data/mixomics_data.xlsx", sheet = "sncRNA")

# In sncRNA, select() the columns containing gene metadata (i.e., non-numerical information) 
#  and store them inside "gene_metadata".
gene_metadata <-
  sncRNA %>% 
  dplyr::select(sRNA.id, Type, sRNA.ID_type, Description)

# In gene_metadata:
# Replace spaces with underscores so that it does not become a problem later
gene_metadata$sRNA.ID_type <- sub(" ", "_", gene_metadata$sRNA.ID_type)
gene_metadata$Type <- sub(" ", "_", gene_metadata$Type)

# Save the gene metadata data frame as a csv
#write.csv(gene_metadata, file = "./1_data/sncRNA/sncRNA_gene_metadata.csv", row.names = FALSE)

# In sncRNA, select() the columns containing sncRNA.ID_type and the raw counts information,
# and store them in "raw_counts". Rename counts columns so the data frame looks cleaner
# in the bar plot to be constructed (e.g., changing A4 to A4_count). 

raw_counts <- 
  sncRNA %>% 
  dplyr::select(sRNA.ID_type, A2 = A2_count, A3 = A3_count, A4 = A4_count, A5 = A5_count, P1 = P1_count, P3 = P3_count, P4 = P4_count, P5 = P5_count)

# Replace spaces with underscores so that it does not become a problem later
raw_counts$sRNA.ID_type <- sub(" ", "_", raw_counts$sRNA.ID_type)

# Save the raw counts data frame as a csv
#write.csv(raw_counts, file = "./1_data/sncRNA/sncRNA_raw_counts.csv", row.names = FALSE)

# In sncRNA, select() the columns containing sncRNA.ID_type and the raw counts information,
#  and store them in "TPM_counts". Rename counts columns so the data frame looks cleaner
TPM_counts <- 
  sncRNA %>% 
  dplyr::select(sRNA.ID_type, A2 = A2_TPM, A3 = A3_TPM, A4 = A4_TPM, A5 = A5_TPM, P1 = P1_TPM, P3 = P3_TPM, P4 = P4_TPM, P5 = P5_TPM)

# Replace spaces with underscores so that it does not become a problem later
TPM_counts$sRNA.ID_type <- sub(" ", "_", TPM_counts$sRNA.ID_type)

# Take the column names from TPM_counts and store them in "sample_IDs"
sample_IDs <- colnames(TPM_counts)[2:ncol(TPM_counts)] 

# Save the TPM counts data frame as a csv
#write.csv(TPM_counts, file = "./1_data/sncRNA/sncRNA_TPM_counts.csv", row.names = FALSE)

# Create a data frame containing the sample ids, call it and Treatment group and 
# store it in the variable "sample_metadata"
sample_metadata <- data.frame("sample.IDs" = sample_IDs) 
sample_metadata <- 
  sample_metadata %>%
  mutate(
    Treatment_group = case_when(
      grepl(pattern = "A", x = sample_IDs) ~ "Acr",
      grepl(pattern = "P", x = sample_IDs) ~ "Control"
    )
  )

# Save the sample metadata data frame as a csv
#write.csv(sample_metadata, file = "./1_data/sncRNA/sncRNA_sample_metadata.csv", row.names = FALSE)


## Initial filtering and quality check (prior to PCA)

# Filter to reduce the amount of noise in that will be present in DeSeq analysis 
# and PCA later. DESeq2 is designed to work with raw count data and internally 
# applies normalization and statistical modeling to account for differences in 
# library size and composition. Hence, the DeSeq later will use raw counts only. 
# However, making QC bar graphs with both raw counts and TPM allows us to obtain 
# a deeper understanding of the data.

# Define thresholds 
min.counts <- 5
min.samples <- 3

# Set the columns as the row names
raw_counts <- 
  raw_counts %>%
  data.frame() %>%
  column_to_rownames(var = "sRNA.ID_type")

### Note: The filtering was not working until I set the row names as the sncRNA
### ids. I think it was counting the numeric row names as the row sums.

raw_counts_filtered <- raw_counts %>%
  filter(rowSums(. >= min.counts) >= min.samples)

# Save the raw_counts data frame for later analysis
#write.csv(raw_counts_filtered, file = "./1_data/sncRNA/sncRNA_raw_counts_filtered.csv", row.names = TRUE)

# Flitering step (above) explained:
# rowSums(. >= min.counts): This part calculates the row sums for each row where 
# the condition >= min.counts is met. In other words, it calculates the number of samples (columns) where the count is greater than or equal to min.counts.
# rowSums(...) returns a numeric vector with the row sums.
# filter(... >= min.samples): This part filters the rows based on the condition 
# that the row sum (number of samples with counts greater than or equal to 
# min.counts) should be greater than or equal to min.samples. 
# Rows that meet this condition are retained, while others are filtered out.

# Prepare data for plotting
raw_colsums <- colSums(raw_counts_filtered)
raw_colsums <- cbind(raw_colsums, sample_metadata)

# Using factor() here and limits = c() in geom_bar to specify that we want the control group on the left
raw_colsums$sample.IDs <- factor(raw_colsums$sample.IDs, levels = c("P1", "P3", "P4", "P5", "A2", "A3", "A4", "A5"))

# Plot the sequencing depth of sncRNAs meeting the raw counts threshold (by sample)
custom_colours <- c("#332288", "#D55E00")

ggplot(raw_colsums, aes(x = sample.IDs, y = raw_colsums, fill = Treatment_group)) +
  geom_bar(stat = "identity", width = 0.4) + scale_fill_manual(values = custom_colours) +
  labs(x = "Sample IDs", y = "Sequencing depth (raw counts)", fill = "Treatment Group") +
  scale_x_discrete(limits = c("P1", "P3", "P4", "P5", "A2", "A3", "A4", "A5"))

#ggsave(filename = "./2_figures/1_rawcounts_QC.png", width = 6, height = 4)

# Finding these sncRNAs in the TPM dataset

# Now that we've identified the rows of sncRNAs meeting the raw counts threshold, 
# in raw_counts, we need to seek these same sncRNAs within the TPM_counts data frame.

# To do that, we need to get the ids of the sncRNAs passing the threshold.
RNAs_passing_QC <- rownames(raw_counts_filtered)

# Prepare the TPM data frame for filtering by setting ids as row names
TPM_counts <- 
  TPM_counts %>%
  data.frame() %>%
  column_to_rownames(var = "sRNA.ID_type")

# Filter and prepare for plotting
TPM_counts_filtered <- TPM_counts[rownames(TPM_counts) %in% RNAs_passing_QC, ]

# (Save for pheatmap later)
#write.csv(TPM_counts_filtered, file = "./1_data/sncRNA/sncRNA_TPM_counts_filtered.csv", row.names = TRUE)

TPM_colsums <- colSums(TPM_counts_filtered)
TPM_colsums <- cbind(TPM_colsums, sample_metadata)
TPM_colsums$sample.IDs <- factor(TPM_colsums$sample.IDs, levels = c("P1", "P3", "P4", "P5", "A2", "A3", "A4", "A5"))

# Plot the TPM sequencing depth of sncRNAs meeting the raw counts threshold (by sample)
custom_colours2 <- c("#CC79A7", "#009E73")

ggplot(TPM_colsums, aes(x = sample.IDs, y = TPM_colsums, fill = Treatment_group)) +
  geom_bar(stat = "identity", width = 0.4) + scale_fill_manual(values = custom_colours2) +
  labs(x = "Sample IDs", y = "Sequencing depth (TPM)", fill = "Treatment Group") +
  scale_x_discrete(limits = c("P1", "P3", "P4", "P5", "A2", "A3", "A4", "A5"))

#ggsave(filename = "./2_figures/2_TPMcounts_QC.png", width = 6, height = 4)



