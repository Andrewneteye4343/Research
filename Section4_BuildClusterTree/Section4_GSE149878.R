# BuildClusterTree
# In order to evaluate the transcriptomics similarity between clusters, we used BuildClusterTree() function to present the similarity

# R version 4.2.2
# Load packages
library(Seurat)
library(ggplot2)

#################################################################################################################
# GSE149878 analysis

# GSE149878_norm.rds can be acquired by following the code in Section1_scRNA_snRNA_DataProcessing.R
GSE149878_norm = readRDS(file = "path/to/GSE149878/GSE149878_norm.rds")

# Subset the monocyte and macrophage cluster in GSE149878_norm
GSE149878_norm = subset(x = GSE149878_norm, subset = label == "CD14-Mo" | label == "CD14CD16-Mo" | label == "PF-M" | label == "PF-M-like")

# Change the default cluster id to our annotation name (ex: CD14-Mo, CD14CD16-Mo, ...)
GSE149878_norm@active.ident = as.factor(GSE149878_norm$label)
GSE149878_norm$seurat_clusters = GSE149878_norm$label

# Use BuildClusterTree() function to construct a phylogenetic tree based on PCA space
GSE149878_norm <- BuildClusterTree(object = GSE149878_norm, reduction = "pca")
GSE149878_tree <- Tool(object = GSE149878_norm, slot = "BuildClusterTree")

# Set the font type as "Arial"
par(family = "Arial")

# Plot the tree
ape::plot.phylo(
  x = GSE149878_tree,
  type = "phylogram",
  direction = "downwards",
  font = 1,
  cex = 2,
  label.offset = 1,
  edge.color = "black",
  edge.width = 5)
title(main = "Cluster Tree - GSE155249", xlab = "", ylab = "", cex.main = 3)

# # This is the end of the BuildClusterTree in GSE149878
#################################################################################################################