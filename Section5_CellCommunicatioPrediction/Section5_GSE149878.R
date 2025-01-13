# R version 4.2.2
# Load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(CellChat)
library(patchwork)
library(igraph)
library(ComplexHeatmap)

#################################################################################################################
# GSE149878 analysis

# Load the dataset
# GSE149878_norm.rds can be acquired by following the code in Section1_scRNA_snRNA_DataProcessing.R
GSE149878_norm = readRDS(file = "path/to/GSE149878/GSE149878_norm.rds")

# Create cellchat
cellchat <- createCellChat(object = GSE149878_norm, group.by = "label")
cellchat <- updateCellChat(cellchat)

# Use Secreted Signaling / ECM-Receptor / Cell-Cell Contact
# In our study, we focus on Secreted Signaling
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
cellchat = subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Check the significant cell-cell communications which occur in the dataset
cellchat@netP$pathways

# To focus on a certain communication, you can set the pathway that you are interested
# The pathway can be acquired with the above codes (ex: cellchat@netP$pathways)
pathways.show = "" 

# Check the contribution of each ligand-receptor interaction in a pathway
netAnalysis_contribution(cellchat, signaling = pathways.show)

# Visualize the computed communication probability
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 10, height = 2, font.size = 10)

# Visualize the communications with Bubble plot
# sources.use represents the ligand providing cell types in communications
# targets.use represents the receptor providing cell types in communications
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:30), remove.isolate = TRUE)
netVisual_bubble(cellchat, sources.use = c(1:30), font.size = 9, targets.use = 15, remove.isolate = TRUE)

# Draw incoming communication with heatmap
incoming_heatmap = netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 5, height = 9)

# This is the end of the Cell-cell Communication Prediction in GSE149878
#################################################################################################################