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
# Spatial transcriptomics from Mendeley Data (DOI:10.17632/xjtv62ncwr.1)
# Load D1 

# D1_norm.rds can be acquired by following the code in Section2_SpatialTranscriptomics_DataProcessing.R
D1_norm = readRDS(file = "path/to/V10F24-110_D1/D1_norm.rds")
data.input = GetAssayData(D1_norm, slot = "data", assay = "SCT")

# manually create a dataframe consisting of the cell labels
meta = data.frame(seurat_clusters = Idents(D1_norm), row.names = names(Idents(D1_norm))) 
meta$label = ""
meta$label[meta$seurat_clusters == "0"] = "AEP" # Alveolar epithelial cell-enriched region
meta$label[meta$seurat_clusters == "1"] = "AEP" # Alveolar epithelial cell-enriched region
meta$label[meta$seurat_clusters == "2"] = "Immune" # Immune cell-enriched region
meta$label[meta$seurat_clusters == "3"] = "CE" # Ciliated epithelial cell-enriched region

spatial.locs = GetTissueCoordinates(D1_norm, scale = NULL, cols = c("imagerow", "imagecol")) 
scale.factors = jsonlite::fromJSON(txt = file.path("path/to/V10F24-110_D1/spatial/", 'scalefactors_json.json'))
scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                     fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
)

# Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "label",
                           datatype = "spatial", coordinates = spatial.locs, scale.factors = scale.factors)
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, distance.use = TRUE, interaction.length = 200, scale.distance = 0.01)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Check the significant cell-cell communications which occur in each group
cellchat@netP$pathways

# To focus on a certain communication, you can set the pathway that you are interested
# The pathway can be acquired with the above codes (ex: cellchat1@netP$pathways)
pathways.show = ""

# Check the contribution of each ligand-receptor interaction in a pathway
netAnalysis_contribution(cellchat, signaling = pathways.show)

# Visualize the computed communication probability
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 6, height = 2.5, font.size = 10, 
                                  color.use = c("#ccebc5","#b3cde3","#fbb4ae"))

# Visualize the communications with Bubble plot
# sources.use represents the ligand providing cell types in communications
# targets.use represents the receptor providing cell types in communications
netVisual_bubble(cellchat, sources.use = c(1:3), targets.use = 1, remove.isolate = FALSE)

# This is the end of the Cell-cell Communication Prediction in Mendeley Data (DOI:10.17632/xjtv62ncwr.1) D1 slide
#################################################################################################################