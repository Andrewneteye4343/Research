# R version 4.3.1
# Load packages
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(GOSemSim)
library(stats)

#################################################################################################################
# Human gene database
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

# We used Differential Expression Genes (DEGs) for enrichment analysis
# DEGs can be acquired from Section8_TGFB1high_monocyte.R or Section9_GAS6_PlasmaCell.R
# Input DEGs data which only contains gene names and log2 fold change
df = read.csv("path/to/YourDirectory/DEGs.csv", header = TRUE)
original_gene_list <- df$FC
names(original_gene_list) <- df$NAME
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

# Set seed as 123
set.seed(123)

# Conduct GSEA with Gene ontology (GO) biological process (BP)
gse <- gseGO(geneList = gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 10, 
             maxGSSize = 500, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism,
             pAdjustMethod = "BH",
             seed = TRUE)

# Save the GSEA results
write.csv(x = gse, file = "path/to/YourDirectory/gsea_result.csv")

# Visualization with Dotplot
# Choose the pathway that you want to visualize
display_pathway_name = c("pathway1","pathway2","pathway3")

dotplot(gse, showCategory = display_pathway_name, font.size = 10, x = "NES")

# This is the end of the Gene Set Enrichment Analysis
#################################################################################################################