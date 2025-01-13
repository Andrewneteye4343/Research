# BuildClusterTree
# In order to evaluate the transcriptomics similarity between clusters, we used BuildClusterTree() function to present the similarity

# R version 4.2.2
# Load packages
library(Seurat)
library(ggplot2)

#################################################################################################################
# GSE171524 analysis

# GSE171524_norm.rds can be acquired by following the code in Section1_scRNA_snRNA_DataProcessing.R
GSE171524_norm = readRDS(file = "path/to/GSE171524/GSE171524_norm.rds")

# Subset the monocyte and macrophage cluster in GSE171524_norm
GSE171524_norm = subset(x = data, subset = label == "CD14CD16-Mo" | label == "HYP-M" | label == "I-Mo" | label == "NA-M" | label == "PF-M" | label == "TR-M" | label == "VCAN-M")

# Change the default cluster id to our annotation name (ex: CD14CD16-Mo, HYP-M, ...)
GSE171524_norm@active.ident = as.factor(GSE171524_norm$label)
GSE171524_norm$seurat_clusters = GSE171524_norm$label

# Use BuildClusterTree() function to construct a phylogenetic tree based on PCA space
GSE171524_norm <- BuildClusterTree(object = GSE171524_norm, reduction = "pca")
GSE171524_tree <- Tool(object = GSE171524_norm, slot = "BuildClusterTree")

# Set the font type as "Arial"
par(family = "Arial")

# Plot the tree
ape::plot.phylo(
  x = GSE171524_tree,
  type = "phylogram",
  direction = "downwards",
  font = 1,
  cex = 2,
  label.offset = 1,
  edge.color = "black",
  edge.width = 5)
title(main = "Cluster Tree - GSE171524", xlab = "", ylab = "", cex.main = 3)

# # This is the end of the BuildClusterTree in GSE171524
#################################################################################################################