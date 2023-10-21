# Written: 6/10/23
# Most recent update: 
# Purpose: GO target prediction
# Author: Chishan Burch
# Contact: chishanburch@gmail.com

#######################################################################################

# Troubleshooting
# Package ver. clusterProfiler 4.8.3 

??clusterProfiler

# # # # # # #

# install required packages
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
#BiocManager::install("AnnotationDbi")

#BiocManager::install(c(
#  "fansi", "vctrs"
#), update = TRUE, ask = FALSE, force = TRUE)

library(clusterProfiler)
library(enrichplot)
library(data.table)
library(dplyr)
library(AnnotationDbi)
library(org.Mm.eg.db)

## This helps load the desired organism database (to install, use the line with hashtag below).

organism = "org.Mm.eg.db"
message("Loading organism databse: '", organism, "'")

#BiocManager::install(organism, character.only = TRUE)

library(organism, character.only = TRUE)

#BiocManager::install(c(
#  "fansi", "vctrs"
#), update = TRUE, ask = FALSE, force = TRUE)



# Read in the ids of mRNAs in the intersection of the venn diagram created in 5_miR_target_prediction.R
pnames <- fread(file = "./1_data/mRNA/predictednames_miUP.csv")
vnames <- fread(file = "./1_data/mRNA/validatednames_miUP.csv")

#vnames2 <- fread(file = "./1_data/mRNA/validatednames_miDOWN.csv")
#pnames2 <- fread(file = "./1_data/mRNA/predictednames_miDOWN.csv")

hit_list <- rbind(pnames, vnames)

#hit_list <- rbind(pnames2, vnames2, pnames, vnames) ##ver2

ids <- hit_list$Symbol
ids <- as.character(ids)

# Obtain ensembl ids for the intersecting miRNA ids (validated and predicted combined)
# We are using ensembl because many genes in background do not have entrez ids
hit_list <- as.vector(mapIds(org.Mm.eg.db, keys = ids, keytype = "SYMBOL", column = "ENSEMBL"))
#hit_list <- as.character(hit_list)
hit_list <- na.omit(hit_list)
str(hit_list)

#x <- c("Tpt1", "Rpl41", "Rpl12", "Selenof", "Rpl27a", "Rpl10a", "Sec62", "Rpl31", "Kcnk1", "Nme2", "Tcim", "Ppig", "Pianp", "Itln1", "9130008F23Rik")
#y <- as.data.frame(mapIds(org.Mm.eg.db, keys = x, keytype = "SYMBOL", column = "ENSEMBL"))

#write.csv(y, file = "./1_data/mRNA/intersecting_mRNA_ids.csv")

# Read in all of the mRNAs entered into deseq. Now take ensembl ids
background <- fread(file = "./1_data/mRNA/multimiR_miRNA_UP.csv")

#background2 <- fread(file = "./1_data/mRNA/multimiR_miRNA_DOWN.csv") ## ver 2

#background <- rbind(background, background2)

background <- background$target_ensembl
background_list <- na.omit(background)# if using target ensembl, na.omit removes like 600 of them lol which is 2/3

#hit_list  = mRNA intersection between upreg miRNA targets vs downreg mRNA
#background_list = mRNA targets from ALL miRNAs that were entered into the DE analysis (this is just the miRNAs after filtering lowly/non expressed).

# Step 2 - Create/run EnrichGo
go_enrich <- enrichGO(gene = hit_list,
                      universe = background_list,
                      OrgDb = organism, 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

# "If there's 'hits' that don't reach p.adj but are p<0.05 you can change the 
# adj.method ="none" and try the tree plot then. Just to have a look."
#go_enrich <- enrichGO(gene = hit_list,
                      universe = background_list,
                      OrgDb = organism, 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10,
                      pAdjustMethod = "none")

# N.B. If not significant terms come out from the initial 'hit_list' using the intersection of mRNA, try using mRNA targets using ALL (up and down reg) miRNAs from the DE analysis.

# Step 3  - Extract and save results table

go_results <- go_enrich@result
go_results$GeneRatio <- gsub("/","//",go_results$GeneRatio)
go_results$BgRatio <- gsub("/","//",go_results$BgRatio)


## Save table

# Step 4 - Plots
## Bar plot
barplot <- barplot(go_enrich,
        title = "GO Biological Pathways")

save.image(barplot, file = "./2_figures/GO_barplot.png")
#use ggsave

## Dot plot
enrichplot::dotplot(go_enrich)

save.image(file = "./2_figures/GO_barplot.png")

## To remove/account for redundant parent pathways (all those terms belonging to the same general functional group)
edox2 <- pairwise_termsim(go_enrich)

### Tree plots
p1 <- treeplot(edox2, offset = 10)

p2 <- treeplot(edox2, hclust_method = "average", offset = 10)

p3 <- treeplot(edox2, offset = 7,
               nCluster = 5, # default is 5. Changes number of groups clustered
               showCategory = 10)  # Number of 'top terms' that will be grouped and displayed in the plot



