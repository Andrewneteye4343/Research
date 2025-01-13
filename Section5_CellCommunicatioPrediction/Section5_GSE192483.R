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

# GSE192483 analysis

# Load the dataset
# GSE192483_norm.rds can be acquired by following the code in Section1_scRNA_snRNA_DataProcessing.R
GSE192483_norm = readRDS(file = "path/to/GSE192483/GSE192483_norm.rds")

# Seperate GSE192483_norm into control and TB group
control = subset(GSE192483_norm, subset = group == "Control")
TB = subset(GSE192483_norm, subset = group == "TB")

# Create cellchat in each group
cellchat1 = createCellChat(object = control, group.by = "label")
cellchat2 = createCellChat(object = TB, group.by = "label")
cellchat1 <- updateCellChat(cellchat1)
cellchat2 <- updateCellChat(cellchat2)

# Use Secreted Signaling / ECM-Receptor / Cell-Cell Contact
# In our study, we focus on Secreted Signaling
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat1@DB <- CellChatDB.use
cellchat1 = subsetData(cellchat1)
cellchat2@DB <- CellChatDB.use
cellchat2 = subsetData(cellchat2)

# Predcit the cell-cell communication with CellChat
cellchat1 <- identifyOverExpressedGenes(cellchat1)
cellchat1 <- identifyOverExpressedInteractions(cellchat1)
cellchat1 <- computeCommunProb(cellchat1)
cellchat1 <- filterCommunication(cellchat1, min.cells = 10)
cellchat1 <- computeCommunProbPathway(cellchat1)
cellchat1 <- aggregateNet(cellchat1)
cellchat1 <- netAnalysis_computeCentrality(cellchat1, slot.name = "netP")

cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)
cellchat2 <- computeCommunProb(cellchat2)
cellchat2 <- filterCommunication(cellchat2, min.cells = 10)
cellchat2 <- computeCommunProbPathway(cellchat2)
cellchat2 <- aggregateNet(cellchat2)
cellchat2 <- netAnalysis_computeCentrality(cellchat2, slot.name = "netP")

# In cases, some cell type is absence in one of the group
# You can optionally use the following code to integrate clusters between groups (if requires)
# You can check the cell type in each group with table() function
table(cellchat1@idents)
table(cellchat2@idents)

# Create annotation in the more cell type group (For instance, cellchat2 has more cell types than cellchat1)
group.new = levels(cellchat2@idents)

# integrate clusters in the few cell type group
cellchat1 <- liftCellChat(cellchat1, group.new)
object.list <- list(Control = cellchat1, COVID19 = cellchat2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat@DB <- CellChatDB.use
cellchat <- updateCellChat(cellchat)

# Check the significant cell-cell communications which occur in each group
cellchat1@netP$pathways
cellchat2@netP$pathways

# To focus on a certain communication, you can set the pathway that you are interested
# The pathway can be acquired with the above codes (ex: cellchat1@netP$pathways)
pathways.show = "" 

# Check the contribution of each ligand-receptor interaction in a pathway
netAnalysis_contribution(cellchat1, signaling = pathways.show)
netAnalysis_contribution(cellchat2, signaling = pathways.show)

# Visualize the computed communication probability in each group
netAnalysis_signalingRole_network(cellchat1, signaling = pathways.show, width = 10, height = 2, font.size = 10)
netAnalysis_signalingRole_network(cellchat2, signaling = pathways.show, width = 10, height = 2, font.size = 10)

# Visualize the communications with Bubble plot
# By setting the max.dataset, you can  present the communications that only possess higher probability in that group
# sources.use represents the ligand providing cell types in communications
# targets.use represents the receptor providing cell types in communications
netVisual_bubble(cellchat, sources.use = 1, font.size = 9, comparison = c(1, 2), targets.use = c(1:30), max.dataset = 2, remove.isolate = TRUE, color.text = c("#00ced1","#fa8072"))
netVisual_bubble(cellchat, sources.use = c(1:30), font.size = 9, comparison = c(1, 2), targets.use = 1, max.dataset = 2, remove.isolate = TRUE, color.text = c("#00ced1","#fa8072"))

# This is the end of the Cell-cell Communication Prediction in GSE192483
#################################################################################################################