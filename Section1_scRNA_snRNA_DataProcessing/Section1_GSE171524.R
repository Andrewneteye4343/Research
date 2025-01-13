# R version 4.2.2
# Load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(ggplot2)
library(ggpubr)

#################################################################################################################
# GSE171524 analysis
# Load GSE171524 data
GSE171524 <- Read10X(data.dir = "path/to/GSE171524/", gene.column = 1)

# Create Seurat object with the setting of min.cell = 3 and min.features = 200
GSE171524_object <- CreateSeuratObject(counts = GSE171524, project = "lung", min.cells = 3, min.features = 200)

# Identify the mitochondria-related genes in GSE171524_object
GSE171524_object[["percent.mt"]] <- PercentageFeatureSet(GSE171524_object, pattern = "^MT-")

# Visualization of nFeature_RNA, nCount_RNA and percent.mt for QC
VlnPlot(GSE171524_object, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
GSE171524_object <- subset(GSE171524_object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)

# Normalization following the LogNormalize method
GSE171524_norm <- NormalizeData(GSE171524_object, normalization.method = "LogNormalize", scale.factor = 10000)

# Scaling the expression with all genes in GSE171524_norm
GSE171524_norm <- FindVariableFeatures(GSE171524_norm, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GSE171524_norm)
GSE171524_norm <- ScaleData(GSE171524_norm, features = all.genes)

# PCA for dimensional reduction with the highly variable features in GSE171524_norm
GSE171524_norm <- RunPCA(GSE171524_norm, features = VariableFeatures(object = GSE171524_norm))

# Find neighbors with top20 PCs and conduct clustering with resolution 0.5
GSE171524_norm <- FindNeighbors(GSE171524_norm, dims = 1:20)
GSE171524_norm <- FindClusters(GSE171524_norm, resolution = 0.5)

# Conduct UMAP and visualization
GSE171524_norm <- RunUMAP(GSE171524_norm, dims = 1:20)
DimPlot(GSE171524_norm, reduction = "umap", shuffle = TRUE)

# Find markers in each clusters, and save the results for manually cell annotation
GSE171524_markers <- FindAllMarkers(GSE171524_norm, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(x = GSE171524_markers, file = "path/to/YourDirectory/marker.csv")

# Manually cell annotation by searching cell markers
GSE171524_norm$label = ""
GSE171524_norm$label[GSE171524_norm$seurat_clusters == "0"] = "Macrophage"
GSE171524_norm$label[GSE171524_norm$seurat_clusters == "1"] = "Tcell"
GSE171524_norm$label[GSE171524_norm$seurat_clusters == "2"] = "Myo-FB" # Myofibroblast
GSE171524_norm$label[GSE171524_norm$seurat_clusters == "3"] = "FB" # Fibroblast
GSE171524_norm$label[GSE171524_norm$seurat_clusters == "4"] = "AEP2" # Alveolar epithelial type 2 cell
GSE171524_norm$label[GSE171524_norm$seurat_clusters == "5"] = "Monocyte"
GSE171524_norm$label[GSE171524_norm$seurat_clusters == "6"] = "AEP2" # Alveolar epithelial type 2 cell
GSE171524_norm$label[GSE171524_norm$seurat_clusters == "7"] = "AEP1" # Alveolar epithelial type 1 cell
GSE171524_norm$label[GSE171524_norm$seurat_clusters == "8"] = "PC" # Plasma cell
GSE171524_norm$label[GSE171524_norm$seurat_clusters == "9"] = "ED" # Endothelial cell
GSE171524_norm$label[GSE171524_norm$seurat_clusters == "10"] = "Macrophage"
GSE171524_norm$label[GSE171524_norm$seurat_clusters == "11"] = "B/DC" # B cell and dendritic cell
GSE171524_norm$label[GSE171524_norm$seurat_clusters == "12"] = "CE" # Ciliated epithelial cell
GSE171524_norm$label[GSE171524_norm$seurat_clusters == "13"] = "Neuron" # Neuron cell
GSE171524_norm$label[GSE171524_norm$seurat_clusters == "14"] = "Mucous" # Mucous cell
GSE171524_norm$label[GSE171524_norm$seurat_clusters == "15"] = "Mast" # Mast cell
GSE171524_norm$label[GSE171524_norm$seurat_clusters == "16"] = "SMC" # Smooth muscle cell
GSE171524_norm$label[GSE171524_norm$seurat_clusters == "17"] = "Cycling" # Cycling cell

# In order to deeply investigate the population of monocyte and macrophage, we subclustered the cluster 0, 5, and 10
# Subset monocyte and macrophage clusters into one dataset
GSE171524_monocyte_lineage = subset(x = GSE171524_norm, subset = seurat_clusters == "0" | seurat_clusters == "5" | seurat_clusters == "10")
GSE171524_monocyte_lineage = FindNeighbors(GSE171524_monocyte_lineage, dims = 1:20)

# Subclustering with resolution 0.3
GSE171524_monocyte_lineage <- FindClusters(GSE171524_monocyte_lineage, resolution = 0.3)

# Conduct UMAP and visualization
GSE171524_monocyte_lineage <- RunUMAP(GSE171524_monocyte_lineage, dims = 1:20, seed.use = 1)
DimPlot(GSE171524_monocyte_lineage, reduction = "umap", label = TRUE, shuffle = TRUE)

# Find markers in each cluster
GSE171524_monocyte_lineage_markers <- FindAllMarkers(GSE171524_monocyte_lineage, only.pos = TRUE, logfc.threshold = 0)
write.csv(x = GSE171524_monocyte_lineage_markers, file = "path/to/YourDirectory/monocyte_lineage_marker.csv")

# Cell annotation in GSE171524_monocyte_lineage
GSE171524_monocyte_lineage$new_label = ""
GSE171524_monocyte_lineage$new_label[GSE171524_monocyte_lineage$RNA_snn_res.0.3 == "0"] = "PF-M" # Profibrotic macrophage
GSE171524_monocyte_lineage$new_label[GSE171524_monocyte_lineage$RNA_snn_res.0.3 == "1"] = "HYP-M" # Hypoxia-associated macrophage
GSE171524_monocyte_lineage$new_label[GSE171524_monocyte_lineage$RNA_snn_res.0.3 == "2"] = "NA-M" # Nerve- and airway-associated macrophage
GSE171524_monocyte_lineage$new_label[GSE171524_monocyte_lineage$RNA_snn_res.0.3 == "3"] = "TR-M" # Tissue-resident macrophage
GSE171524_monocyte_lineage$new_label[GSE171524_monocyte_lineage$RNA_snn_res.0.3 == "4"] = "VCAN-M" # VCAN+ macrophage
GSE171524_monocyte_lineage$new_label[GSE171524_monocyte_lineage$RNA_snn_res.0.3 == "5"] = "CD14CD16-Mo" # CD14+CD16+ monocyte

# Rename the cell annotation in GSE171524_norm based our subclustering result
PF_M = subset(GSE171524_monocyte_lineage, subset = new_label == "PF-M")
HYP_M = subset(GSE171524_monocyte_lineage, subset = new_label == "HYP-M")
NA_M = subset(GSE171524_monocyte_lineage, subset = new_label == "NA-M")
TR_M = subset(GSE171524_monocyte_lineage, subset = new_label == "TR-M")
VCAN_M = subset(GSE171524_monocyte_lineage, subset = new_label == "VCAN-M")
CD14CD16_Mo = subset(GSE171524_monocyte_lineage, subset = new_label == "CD14CD16-Mo")

allcellname = colnames(GSE171524_norm)
for (i in 1:length(allcellname)){
  if (allcellname[i] %in% colnames(PF_M)){
    GSE171524_norm$label[i] = "PF-M"
  }
  else if (allcellname[i] %in% colnames(HYP_M)){
    GSE171524_norm$label[i] = "HYP-M"
  }
  else if (allcellname[i] %in% colnames(NA_M)){
    GSE171524_norm$label[i] = "NA-M"
  }
  else if (allcellname[i] %in% colnames(TR_M)){
    GSE171524_norm$label[i] = "TR-M"
  }
  else if (allcellname[i] %in% colnames(VCAN_M)){
    GSE171524_norm$label[i] = "VCAN-M"
  }
  else if (allcellname[i] %in% colnames(CD14CD16_Mo)){
    GSE171524_norm$label[i] = "CD14CD16-Mo"
  }
}

# Label with control and disease group
id = str_extract(string = colnames(GSE171524_norm), pattern = "[_0-9]+")
GSE171524_norm$id = id
GSE171524_norm$group = ""
GSE171524_norm$group[GSE171524_norm$id == "1_1"] = "Control"
GSE171524_norm$group[GSE171524_norm$id == "1_2"] = "Control"
GSE171524_norm$group[GSE171524_norm$id == "1_3"] = "Control"
GSE171524_norm$group[GSE171524_norm$id == "1_4"] = "Control"
GSE171524_norm$group[GSE171524_norm$id == "1_5"] = "Control"
GSE171524_norm$group[GSE171524_norm$id == "1_6"] = "Control"
GSE171524_norm$group[GSE171524_norm$id == "1_7"] = "Control"
GSE171524_norm$group[GSE171524_norm$id == "1_8"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_9"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_10"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_11"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_12"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_13"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_14"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_15"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_16"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_17"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_18"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_19"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_20"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_21"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_22"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_23"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_24"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_25"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_26"] = "COVID-19"
GSE171524_norm$group[GSE171524_norm$id == "1_27"] = "COVID-19"

# This is the end of the Data Processing in GSE171524
#################################################################################################################
