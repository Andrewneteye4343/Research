# R version 4.2.2
# Load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(ggplot2)
library(ggpubr)

#################################################################################################################
# GSE149878 analysis
# Load GSE149878 data
data1 = Read10X_h5(filename = "path/to/GSE149878/GSM4516279_C166.h5")
data2 = Read10X_h5(filename = "path/to/GSE149878/GSM4516280_C168.h5")
data3 = Read10X_h5(filename = "path/to/GSE149878/GSM4516281_C170.h5")
data4 = Read10X_h5(filename = "path/to/GSE149878/GSM4516282_C172.h5")

# Create Seurat object with the setting of min.cell = 3 and min.features = 200
data1 = CreateSeuratObject(counts = data1, project = "C166", min.cells = 3, min.features = 200)
data2 = CreateSeuratObject(counts = data2, project = "C168", min.cells = 3, min.features = 200)
data3 = CreateSeuratObject(counts = data3, project = "C170", min.cells = 3, min.features = 200)
data4 = CreateSeuratObject(counts = data4, project = "C172", min.cells = 3, min.features = 200)

# Merge all the seurat objects
GSE149878_object = merge(x = data1, y = c(data2,data3,data4),
                         add.cell.ids = c("C166","C168","C170","C172"), project = "merge_data")

# Identify the mitochondria-related genes in GSE149878_object
GSE149878_object[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

# Visualization of nFeature_RNA, nCount_RNA and percent.mt for QC
VlnPlot(object = GSE149878_object, features = c("nCount_RNA","nFeature_RNA","percent.mt"))
GSE149878_object <- subset(GSE149878_object, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)

# Normalization following the LogNormalize method
GSE149878_norm <- NormalizeData(GSE149878_object, normalization.method = "LogNormalize", scale.factor = 10000)

# Scaling the expression with all genes in GSE155249_norm
GSE149878_norm <- FindVariableFeatures(GSE149878_norm, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GSE149878_norm)
GSE149878_norm <- ScaleData(GSE149878_norm, features = all.genes)

# PCA for dimensional reduction with the highly variable features in GSE149878_norm
GSE149878_norm <- RunPCA(GSE149878_norm, features = VariableFeatures(object = GSE149878_norm))

# Find neighbors with top20 PCs and conduct clustering with resolution 0.5
GSE149878_norm <- FindNeighbors(GSE149878_norm, dims = 1:20)
GSE149878_norm <- FindClusters(GSE149878_norm, resolution = 0.5)

# Conduct UMAP and visualization
GSE149878_norm <- RunUMAP(GSE149878_norm, dims = 1:20, seed.use = 1)
DimPlot(GSE149878_norm, reduction = "umap", shuffle = TRUE, raster = FALSE)

# Find markers in each clusters, and save the results for manually cell annotation
GSE149878_markers = FindAllMarkers(object = GSE149878_norm, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(x = GSE149878_markers, file = "path/to/YourDirectory/marker.csv")

# Manually cell annotation by searching cell markers
GSE149878_norm$label = ""
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "0"] = "Neutrophil"
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "1"] = "Neutrophil"
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "2"] = "CD8T" # CD8+ T cell
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "3"] = "PF-M" # Profibrotic macrophage
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "4"] = "CD14-Mo" # CD14+CD16- monocyte
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "5"] = "CD4T" # CD4+ T cell
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "6"] = "CD14CD16-Mo" # CD14+CD16+ monocyte
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "7"] = "PF-M-like" # Profibrotic macrophage-like cell
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "8"] = "AEP" # Alveolar epithelial cell
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "9"] = "NK" # Nature killer cell
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "10"] = "CE" # Ciliated epithelial cell
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "11"] = "Mast" # Mast cell
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "12"] = "RB" # Red blood cell
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "13"] = "AEP" # Alveolar epithelial cell
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "14"] = "AEP" # Alveolar epithelial cell
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "15"] = "PC" # Plasma cell
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "16"] = "Neutrophil"
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "17"] = "Myo-FB" # Myofibroblast
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "18"] = "Bcell"
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "19"] = "ED" # Endothelial cell
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "20"] = "Cycling" # Cycling cell
GSE149878_norm$label[GSE149878_norm$seurat_clusters == "21"] = "Mast" # Mast cell

# This is the end of the Data Processing in GSE149878
#################################################################################################################