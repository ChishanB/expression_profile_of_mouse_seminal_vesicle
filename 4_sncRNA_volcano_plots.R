
# Written: 9/10/23
# Most recent update: 20/10/23

# Purpose: Constructing volcano plots to represent differentially expressed miRNA, piRNA and rRNA 

# Author: Chishan Burch
# Contact: chishanburch@gmail.com

##################################

# Load packages
library(data.table)
library(tibble)
library(dplyr)
library(ggplot2)


# Set thresholds
p_adj_threshold <- 0.05
log2fc_threshold_low <- -1
log2fc_threshold_high <- 1

# # Read in the sncRNA DeSeq data (incl. those meeting the significance threshold and those not)
sncRNA_DeSeq <- fread("./1_data/sncRNA/sncRNA_DeSeq_complete.csv")

#sncRNA_DeSeq <- 
#  sncRNA_DeSeq %>%
#  column_to_rownames(var = "V1")

# (dplyr method) Use the mutate function and a nested if else statement to create a new column indicating whether
# or not the RNA meets the threshold padj < 0.05)
sncRNA_DeSeq <- sncRNA_DeSeq %>%
  mutate(
    "Statistical_significance" = ifelse(padj < p_adj_threshold, "Significant",
                                        ifelse(padj >= p_adj_threshold, "Not_Significant", "NS")
    ))
# using dplyr again, make a column indicating whether they are being upregulated or downregulated
sncRNA_DeSeq <- sncRNA_DeSeq %>%
  mutate(
    "Expression" = ifelse(log2FoldChange >= log2fc_threshold_high, "UP_Reg",
                          ifelse(log2FoldChange <= log2fc_threshold_low, "DOWN_Reg", "NS")
    ))

# saving the csv
#write.csv(sncRNA_DeSeq, file = "./1_data/sncRNA/sncRNA_DeSeq_complete.csv", row.names = TRUE)

#volcano_plot <-ggplot(sncRNA_DeSeq, aes(x = log2FoldChange, y = -log10(padj), col = Statistical_significance)) + 
#  geom_point(size = 1) +
#  geom_hline(yintercept = -log10(p_adj_threshold), col = "#000000", linetype = "dotted", size = 1) +
 # scale_color_manual(values = c('#D55E00', '#009E73') ) + geom_vline(xintercept = c(log2fc_threshold_low, log2fc_threshold_high), 
 #                     col = "#000000", linetype = "dotted", linewidth = 1)
??geom_hline


# The goal now is to make volcano plots with the miRNAs, snoRNAs and piRNAs

# Read the file in fresh 
sncRNA_DeSeq <- fread("./1_data/sncRNA/sncRNA_DeSeq_complete.csv")

sncRNA_DeSeq <- sncRNA_DeSeq %>%
  mutate(
    "Statistical_significance" = ifelse(padj < p_adj_threshold, "Significant",
                                        ifelse(padj >= p_adj_threshold, "Not_Significant", "NS")
    ))
# using dplyr again, make a column indicating whether they are being upregulated or downregulated
sncRNA_DeSeq <- sncRNA_DeSeq %>%
  mutate(
    "Expression" = ifelse(log2FoldChange >= log2fc_threshold_high, "UP_Reg",
                          ifelse(log2FoldChange <= log2fc_threshold_low, "DOWN_Reg", "NS")
    ))




# If needed, Remove the underscored so we can grep what we want
#y <- sncRNA_DeSeq %>% dplyr::select(V1)%>%
#  .[["V1"]] %>%
# gsub("_", " ", .) 

#z <- cbind(y, sncRNA_DeSeq)

#sncRNA_DeSeq <- setnames(z, "y", "miRNA_ids_spaces")

str()

sncRNA_DeSeq <-
 sncRNA_DeSeq %>%
  mutate(
    Type = case_when(
      grepl(pattern = "hairpin", x = miRNA_ids_spaces) ~ "hairpin",
      grepl(pattern = "mature", x = miRNA_ids_spaces) ~ "mature",
      grepl(pattern = "piRNA", x = miRNA_ids_spaces) ~ "piRNA",
      grepl(pattern = "snoRNA", x = miRNA_ids_spaces) ~ "snoRNA",
      grepl(pattern = "snRNA", x = miRNA_ids_spaces) ~ "snRNA",
      grepl(pattern = "rRNA", x = miRNA_ids_spaces) ~ "rRNA",
      grepl(pattern = "Rfam other", x = miRNA_ids_spaces) ~ "Rfam other than sncRNA",
   )
  )

tail(sncRNA_DeSeq)

# Adding a -log10(padj) column to more easily pick out specific miRNAs
sncRNA_DeSeq <- sncRNA_DeSeq %>% mutate(log_padj = -log10(padj))

# Identify the outliers meeting padj thresholds
x <- sncRNA_DeSeq %>%
  filter(padj < 0.05)
x <- x[Type == "mature"]

x %>% arrange(desc(log_padj))

# Again for rRNAs
y <- sncRNA_DeSeq %>%
  filter(padj < 0.05)
y <- y[Type == "rRNA"]


y <- sncRNA_DeSeq[Type == "rRNA"]
y <- y %>%
  filter(padj < 0.05)

y %>% arrange(desc(log_padj))


# saving the csv
#write.csv(sncRNA_DeSeq, file = "./1_data/sncRNA/sncRNA_DeSeq_complete.csv", row.names = FALSE)

# Now to construct volcano plots for some sncRNA sub categories

mature_miRNAs <- sncRNA_DeSeq[Type == "mature"]

mature_miRNAs <- mature_miRNAs %>%
  mutate(sign_DE = if_else(padj < 0.05 & log2FoldChange <= -1, "Sign_Down",
                           if_else(padj < 0.05 & log2FoldChange >= 1, "Sign_Up", "NS"))) 

#hairpin_miRNAs <- sncRNA_DeSeq[Type == "hairpin"]
#all_miRNAs <- rbind(mature_miRNAs, hairpin_miRNAs)


#piRNAs <- sncRNA_DeSeq[Type == "piRNA"]
#snoRNAs <- sncRNA_DeSeq[Type == "snoRNA"]
#rRNAs <- sncRNA_DeSeq[Type == "rRNA"]

# Declare colours and opacity

# miRNA
colours_1 <- c("Sign_Up" = "#EE6677", "Sign_Down" = "#0072B2", "NS" = "#BBBBBB")
opacity <- c("Sign_Up" = 0.5, "Sign_Down" = 0.5, "NS" = 0.5)

colours_factor <- factor(mature_miRNAs$sign_DE, 
                        levels = c("#EE6677", "#0072B2", "#BBBBBB"),
                        labels = c("Sign_Up", "Sign_Down", "NS"))

# mature miRNA plot
mature_miRNAs <- as.data.frame(mature_miRNAs)


mature_miRNAs$sign_DE <- factor(mature_miRNAs$sign_DE, levels = c("Sign_Up","Sign_Down", "NS"))

mature_miRNAs 

# IT WORKS!!

miRNA_plot <- mature_miRNAs %>%
  filter(is.na(sign_DE) != TRUE) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), col = sign_DE, alpha = sign_DE, fill = sign_DE)) +
  geom_point(size = 5, shape = 21, color = "black") +
  geom_hline(yintercept = -log10(p_adj_threshold), linetype = "dashed", size = 0.4, col = "#3A3B3C") +
  geom_vline(xintercept = c(log2fc_threshold_low, log2fc_threshold_high), linetype = "dashed", size = 0.4, col = "#3A3B3C") +
  scale_color_manual(values = colours_1) + 
  scale_fill_manual(values = colours_1) +
  scale_alpha_manual(values = opacity) + theme_minimal() #+ geom_label_repel(aes(label = miRNA_ids_underscores))

# Can we try rRNA ver


                           


ggsave(plot = miRNA_plot, filename = "./2_figures/miRNA_volcano_plot.png", width = 8, height = 6, dpi = 300)
  
# snoRNA plot
colours_2 <- c("UP_Reg" = "#D55E00", "DOWN_Reg" = "#0072B2", "NS" = "#BBBBBB")
opacity <- c("UP_Reg" = 0.5, "DOWN_Reg" = 0.5, "NS" = 0.1)


snoRNA_plot <- snoRNAs %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), col = Expression, alpha = Expression, fill = Expression)) +
  geom_point(size = 6, shape = 21, color = "black") +
  geom_hline(yintercept = -log10(p_adj_threshold), linetype = "dotted", size = 0.8) +
  geom_vline(xintercept = c(log2fc_threshold_low, log2fc_threshold_high), linetype = "dotted", size = 0.8) +
  scale_color_manual(values = colours_2) + 
  scale_fill_manual(values = colours_2) +
  scale_alpha_manual(values = opacity) + theme_minimal()

ggsave(plot = snoRNA_plot, filename = "./2_figures/snoRNA_volcano_plot.png", width = 8, height = 6, dpi = 300)


colours_3 <- c("UP_Reg" = "#F0E442", "DOWN_Reg" = "#0072B2", "NS" = "#BBBBBB")
opacity <- c("UP_Reg" = 0.5, "DOWN_Reg" = 0.5, "NS" = 0.1)


piRNA_plot <- piRNAs %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), col = Expression, alpha = Expression, fill = Expression)) +
  geom_point(size = 6, shape = 21, color = "black") +
  geom_hline(yintercept = -log10(p_adj_threshold), linetype = "dashed", size = 0.4, col = "#3A3B3C") +
  geom_vline(xintercept = c(log2fc_threshold_low, log2fc_threshold_high), linetype = "dashed", size = 0.4, col = "#3A3B3C") +
  scale_color_manual(values = colours_3) + 
  scale_fill_manual(values = colours_3) +
  scale_alpha_manual(values = opacity) + theme_minimal() 



ggsave(plot = piRNA_plot, filename = "./2_figures/piRNA_volcano_plot.png", width = 8, height = 6, dpi = 300)


# miRNA plot 2
all_miRNAs <- rbind(mature_miRNAs, hairpin_miRNAs)

# miRNA
colours_1 <- c("UP_Reg" = "#EE6677", "DOWN_Reg" = "#0072B2", "NS" = "#BBBBBB")
opacity <- c("UP_Reg" = 0.5, "DOWN_Reg" = 0.5, "NS" = 0.1)




total_miRNA_plot <- all_miRNAs %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), col = Expression, alpha = Expression, fill = Expression)) +
  geom_point(size = 6, shape = 21, color = "black") +
  geom_hline(yintercept = -log10(p_adj_threshold), linetype = "dashed", size = 0.4, col = "#3A3B3C") +
  geom_vline(xintercept = c(log2fc_threshold_low, log2fc_threshold_high), linetype = "dashed", size = 0.4, col = "#3A3B3C") +
  scale_color_manual(values = colours_1) + 
  scale_fill_manual(values = colours_1) +
  scale_alpha_manual(values = opacity) + theme_minimal()

all_miRNAs %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = Expression, alpha = Expression, col = Expression)) +
  geom_point(size = 6, shape = 21) +
  geom_hline(yintercept = -log10(p_adj_threshold), linetype = "dotted", size = 0.8) +
  geom_vline(xintercept = c(log2fc_threshold_low, log2fc_threshold_high), linetype = "dotted", size = 0.8) +
  scale_fill_manual(values = colours_1) +
  scale_color_manual(values = colours_1) +
  scale_alpha_manual(values = opacity) +
  theme_minimal() +
  geom_point(data = subset(mature_miRNAs, padj <= 0.05),  # Subset for points below threshold
             aes(fill = "Below Threshold"),  # Map to a constant value
             size = 6, shape = 21, color = "red") +  # Adjust the appearance
  guides(fill = guide_legend(title = "Expression", override.aes = list(shape = 21, color = "green")))  # Customize the legend


########################




#From Paul Tol: https://personal.sron.nl/~pault/
Tol_bright <- c('#EE6677', '#228833', '#4477AA', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB')
Tol_muted <- c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', '#882255', '#AA4499', '#DDDDDD')
Tol_light <- c('#BBCC33', '#AAAA00', '#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF', '#44BB99', '#DDDDDD')
#From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")


# snoRNA plot
snoRNA_plot <-ggplot(snoRNAs, aes(x = log2FoldChange, y = -log10(padj), col = Statistical_significance)) + 
  geom_point(size = 3) +
  geom_hline(yintercept = -log10(p_adj_threshold), col = "#000000", linetype = "dotted", size = 1) +
  scale_color_manual(values = c('#D55E00', '#009E73') ) + geom_vline(xintercept = c(log2fc_threshold_low, log2fc_threshold_high), 
                                                                     col = "#000000", linetype = "dotted", linewidth = 1)

# piRNA plot
#piRNA_plot <-ggplot(piRNAs, aes(x = log2FoldChange, y = -log10(padj), col = Statistical_significance)) + 
  geom_point(size = 10) +
  geom_hline(yintercept = -log10(p_adj_threshold), col = "#000000", linetype = "dotted", size = 1) +
  scale_color_manual(values = c('#D55E00', '#009E73') ) + geom_vline(xintercept = c(log2fc_threshold_low, log2fc_threshold_high), 
                                                                     col = "#000000", linetype = "dotted", linewidth = 1)

# rRNA plot
#rRNA_plot <-ggplot(rRNAs, aes(x = log2FoldChange, y = -log10(padj), col = Statistical_significance)) + 
  geom_point(size = 1) +
  geom_hline(yintercept = -log10(p_adj_threshold), col = "#000000", linetype = "dotted", size = 1) +
  scale_color_manual(values = c('#D55E00', '#009E73') ) + geom_vline(xintercept = c(log2fc_threshold_low, log2fc_threshold_high), 
                                                                     col = "#000000", linetype = "dotted", linewidth = 1)


