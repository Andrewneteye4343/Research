# BuildClusterTree
# In order to evaluate the transcriptomics similarity between clusters, we used BuildClusterTree() function to present the similarity

# R version 4.2.2
# Load packages
library(Seurat)
library(ggplot2)

#################################################################################################################
# GSE155249 analysis

# GSE155249_norm.rds can be acquired by following the code in Section1_scRNA_snRNA_DataProcessing.R
GSE155249_norm = readRDS(file = "path/to/GSE155249/GSE155249_norm.rds")

# Subset the monocyte and macrophage cluster in GSE155249_norm
GSE155249_norm = subset(x = GSE155249_norm, subset = label == "PF-M" | label == "TR-M" | label == "CD14CD16-Mo1" | label == "CD14CD16-Mo2" | label == "CD14CD16-Mo3")

# Change the default cluster id to our annotation name (ex: PF-M, TR-M, ...)
GSE155249_norm@active.ident = as.factor(GSE155249_norm$label)
GSE155249_norm$seurat_clusters = GSE155249_norm$label

# Use BuildClusterTree() function to construct a phylogenetic tree based on PCA space
GSE155249_norm <- BuildClusterTree(object = GSE155249_norm, reduction = "pca")
GSE155249_tree <- Tool(object = GSE155249_norm, slot = "BuildClusterTree")


# Set the font type as "Arial"
par(family = "Arial")

# Plot the tree
ape::plot.phylo(
  x = GSE155249_tree,
  type = "phylogram",
  direction = "downwards",
  font = 1,
  cex = 2,
  label.offset = 1,
  edge.color = "black",
  edge.width = 5)
title(main = "Cluster Tree - GSE155249", xlab = "", ylab = "", cex.main = 3)

# # This is the end of the BuildClusterTree in GSE155249
#################################################################################################################