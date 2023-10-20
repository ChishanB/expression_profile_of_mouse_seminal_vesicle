
# Written: 5/10/23
# Most recent update: 20/10/23

# Purpose: Comparing proportions of significantly upregulated and downregulated 
# sncRNAs (adjusted p-value < 0.05) and constructing pie charts to represent this.

# Author: Chishan Burch
# Contact: chishanburch@gmail.com

##################################

# Load packages
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)

# Obtain the filtered, significant DEseq results
sncRNA_counts <- fread("./1_data/sncRNA/sncRNA_DeSeq_padj_under_0.05.csv") # Created in 2_PCA_and_DeSeq.R
gene_metadata <- fread("./1_data/sncRNA/sncRNA_gene_metadata.csv")

sample_metadata <- fread("./1_data/sncRNA/sncRNA_sample_metadata.csv") %>% 
  as.data.frame() %>%
  mutate(Treatment_group = factor(Treatment_group, levels = c("Control", "Acr"))) %>%
  column_to_rownames(., var = "sample.IDs")

sncRNA_counts <- sncRNA_counts %>%
  rename("V1" = "sRNA.ID_type")


sncRNA_counts <- sncRNA_counts %>%
  mutate(
    Expression = ifelse(log2FoldChange > 0, "UP",
                        ifelse(log2FoldChange < 0, "DOWN", "Meets_no_criteria")
    ))

# Obtain no. upregulated RNAs
nrow(sncRNA_counts[sncRNA_counts$log2FoldChange > 0, ])
#45 upregulated
nrow(sncRNA_counts[sncRNA_counts$log2FoldChange < 0, ])
#11 downregulated

UPREG <- sncRNA_counts %>% filter(sncRNA_counts$log2FoldChange > 0)
DOWNREG <- sncRNA_counts %>% filter(sncRNA_counts$log2FoldChange < 0)


# Prepare data for pie charts
UP_sncRNAs <- gene_metadata %>%
  filter(gene_metadata$sRNA.ID_type %in% UPREG$sRNA.ID_type) %>%
  data.frame() %>%
  column_to_rownames(var = "sRNA.ID_type") %>%
  dplyr::select(Type)

#Get number
types_table_1 <- table(UP_sncRNAs$Type)

DOWN_sncRNAs <- gene_metadata %>%
  filter(gene_metadata$sRNA.ID_type %in% DOWNREG$sRNA.ID_type) %>%
  data.frame() %>%
  column_to_rownames(var = "sRNA.ID_type") %>%
  dplyr::select(Type)

#Get number
types_table_2 <- table(DOWN_sncRNAs$Type)


UP_types <- data.frame(Type = names(types_table_1), Count = as.numeric(types_table_1))
DOWN_types <- data.frame(Type = names(types_table_2), Count = as.numeric(types_table_2))

# Plot upregulated sncRNAs
  ## Declare custom colours
custom_colors <- c(hairpin = "#0072B2", mature = "#CC79A7", piRNA = "#F0E442", precursor = "#882255", 
                               rRNA = "#009E73", snoRNA = "#E69F00")

UP_sncRNA_plot <- 
  ggplot(UP_types, aes(x = " ", y = Count, fill = Type), size = 13) +
  geom_bar(stat = "identity") +
  coord_polar(theta = "y") +
  labs(title = "Proportions of upregulated sncRNA subcategories", fill = "Type") + theme_void() +
  geom_text(aes(label = Count, y = Count + runif(length(Count), -0.2, 0.2)), position = position_stack(vjust = 0.5), size = 5) +
  scale_fill_manual(values = custom_colors) + theme(legend.key.size = unit(1, "cm"), legend.title = element_text(size=14),
  legend.text = element_text(size = 10)) 

DOWN_sncRNA_plot <- 
  ggplot(DOWN_types, aes(x = " ", y = Count, fill = Type), size = 20) +
  geom_bar(stat = "identity") +
  coord_polar(theta = "y") +
  labs(title = "Proportions of downregulated sncRNA subcategories", fill = "Type") + theme_void() +
  geom_text(aes(label = Count, y = Count + runif(length(Count), -0.2, 0.2)), position = position_stack(vjust = 0.5), size = 5) +
  scale_fill_manual(values = custom_colors) + theme(legend.key.size = unit(1, "cm"), legend.title = element_text(size=14),
  legend.text = element_text(size = 10)) 

# save plot1 and plot2
#ggsave(plot = UP_sncRNA_plot, filename = "./2_figures/upregulated_proportions_sncRNA.png", width = 6, height = 6, dpi = 300)
#ggsave(plot = DOWN_sncRNA_plot, filename = "./2_figures/downregulated_proportions_sncRNA.png", width = 6, height = 6, dpi = 300)


