# R version 4.2.2
# Load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(ggplot2)
library(ggpubr)

#################################################################################################################
# GSE155249 analysis
# Load GSE155249 data
data1 = Read10X_h5(filename = "path/to/GSE155249/GSM4698176.h5")
data2 = Read10X_h5(filename = "path/to/GSE155249/GSM4698177.h5")
data3 = Read10X_h5(filename = "path/to/GSE155249/GSM4698178.h5")
data4 = Read10X_h5(filename = "path/to/GSE155249/GSM4698179.h5")
data5 = Read10X_h5(filename = "path/to/GSE155249/GSM4698180.h5")
data6 = Read10X_h5(filename = "path/to/GSE155249/GSM4698181.h5")
data7 = Read10X_h5(filename = "path/to/GSE155249/GSM4698182.h5")
data8 = Read10X_h5(filename = "path/to/GSE155249/GSM4698183.h5")
data9 = Read10X_h5(filename = "path/to/GSE155249/GSM4698184.h5")
data10 = Read10X_h5(filename = "path/to/GSE155249/GSM4850111.h5")
data11 = Read10X_h5(filename = "path/to/GSE155249/GSM4850112.h5")
data12 = Read10X_h5(filename = "path/to/GSE155249/GSM4850113.h5")
data13 = Read10X_h5(filename = "path/to/GSE155249/GSM4850114.h5")
data14 = Read10X_h5(filename = "path/to/GSE155249/GSM4850115.h5")
data15 = Read10X_h5(filename = "path/to/GSE155249/GSM4850116.h5")
data16 = Read10X_h5(filename = "path/to/GSE155249/GSM4850118.h5")
data17 = Read10X_h5(filename = "path/to/GSE155249/GSM4850119.h5")
data18 = Read10X_h5(filename = "path/to/GSE155249/GSM4850120.h5")

# Create Seurat object with the setting of min.cell = 3 and min.features = 200
data1 = CreateSeuratObject(counts = data1, project = "8176", min.cells = 3, min.features = 200)
data2 = CreateSeuratObject(counts = data2, project = "8177", min.cells = 3, min.features = 200)
data3 = CreateSeuratObject(counts = data3, project = "8178", min.cells = 3, min.features = 200)
data4 = CreateSeuratObject(counts = data4, project = "8179", min.cells = 3, min.features = 200)
data5 = CreateSeuratObject(counts = data5, project = "8180", min.cells = 3, min.features = 200)
data6 = CreateSeuratObject(counts = data6, project = "8181", min.cells = 3, min.features = 200)
data7 = CreateSeuratObject(counts = data7, project = "8182", min.cells = 3, min.features = 200)
data8 = CreateSeuratObject(counts = data8, project = "8183", min.cells = 3, min.features = 200)
data9 = CreateSeuratObject(counts = data9, project = "8184", min.cells = 3, min.features = 200)
data10 = CreateSeuratObject(counts = data10, project = "0111", min.cells = 3, min.features = 200)
data11 = CreateSeuratObject(counts = data11, project = "0112", min.cells = 3, min.features = 200)
data12 = CreateSeuratObject(counts = data12, project = "0113", min.cells = 3, min.features = 200)
data13 = CreateSeuratObject(counts = data13, project = "0114", min.cells = 3, min.features = 200)
data14 = CreateSeuratObject(counts = data14, project = "0115", min.cells = 3, min.features = 200)
data15 = CreateSeuratObject(counts = data15, project = "0116", min.cells = 3, min.features = 200)
data16 = CreateSeuratObject(counts = data16, project = "0118", min.cells = 3, min.features = 200)
data17 = CreateSeuratObject(counts = data17, project = "0119", min.cells = 3, min.features = 200)
data18 = CreateSeuratObject(counts = data18, project = "0120", min.cells = 3, min.features = 200)

# Merge all the seurat objects
GSE155249_object = merge(x = data1, y = c(data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14,data15,data16,data17,data18),
                         add.cell.ids = c("8176","8177","8178","8179","8180","8181","8182","8183","8184","0111","0112","0113","0114","0115","0116","0118","0119","0120"), project = "merge_data")

# Rename gene name from Ensembl to Gene Symbol
newname = read.csv("path/to/GSE155249/new_GSE155249_genename.csv", header = TRUE, row.names = 1)
newname = as.data.frame(newname)
colnames(newname) = "gene"
GSE155249_object@assays$RNA@counts@Dimnames[[1]] <- newname$gene
GSE155249_object@assays$RNA@data@Dimnames[[1]] <- newname$gene
GSE155249_object[["RNA"]]@meta.features <- data.frame(row.names = rownames(GSE155249_object[["RNA"]]))

# Identify the mitochondria-related genes in GSE155249_object
GSE155249_object[["percent.mt"]] <- PercentageFeatureSet(GSE155249_object, pattern = "^MT-")

# Visualization of nFeature_RNA, nCount_RNA and percent.mt for QC
VlnPlot(object = GSE155249_object, features = c("nCount_RNA","nFeature_RNA","percent.mt"))
GSE155249_object <- subset(GSE155249_object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)

# Normalization following the LogNormalize method
GSE155249_norm <- NormalizeData(GSE155249_object, normalization.method = "LogNormalize", scale.factor = 10000)

# Scaling the expression with all genes in GSE155249_norm
GSE155249_norm <- FindVariableFeatures(GSE155249_norm, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GSE155249_norm)
GSE155249_norm <- ScaleData(GSE155249_norm, features = all.genes)

# PCA for dimensional reduction with the highly variable features in GSE155249_norm
GSE155249_norm <- RunPCA(GSE155249_norm, features = VariableFeatures(object = GSE155249_norm))

# Find neighbors with top20 PCs and conduct clustering with resolution 0.5
GSE155249_norm <- FindNeighbors(GSE155249_norm, dims = 1:20)
GSE155249_norm <- FindClusters(GSE155249_norm, resolution = 0.5)

# Conduct UMAP and visualization
GSE155249_norm <- RunUMAP(GSE155249_norm, dims = 1:20, seed.use = 1)
DimPlot(GSE155249_norm, reduction = "umap", shuffle = TRUE, raster = FALSE)

# Find markers in each clusters, and save the results for manually cell annotation
GSE155249_markers = FindAllMarkers(object = GSE155249_norm, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(x = GSE155249_markers, file = "path/to/YourDirectory/marker.csv")

# Manually cell annotation by searching cell markers
GSE155249_norm$label = ""
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "0"] = "PF-M" # Profibrotic macrophage
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "1"] = "CD14CD16-Mo1" # CD14+CD16+ monocyte population 1
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "2"] = "TR-M" # Tissue-resident macrophage
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "3"] = "CD8T" # CD8+ T cell
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "4"] = "CD14CD16-Mo2"# CD14+CD16+ monocyte population 2
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "5"] = "CD4T" # CD4+ T cell
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "6"] = "TR-M" # Tissue-resident macrophage
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "7"] = "TR-M" # Tissue-resident macrophage
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "8"] = "DC" # Dendritic cell
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "9"] = "CD14CD16-Mo3" # CD14+CD16+ monocyte population 3
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "10"] = "Pro-CD8T" # Proliferative CD8+ T cell
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "11"] = "PC" # Plasma cell
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "12"] = "CD4T" # CD4+ T cell
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "13"] = "pDC" # Plasmacytoid dendritic cell
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "14"] = "CD8T" # CD8+ T cell
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "15"] = "CE" # Ciliated epithelial cell
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "16"] = "CD8T" # CD8+ T cell
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "17"] = "Bcell" 
GSE155249_norm$label[GSE155249_norm$seurat_clusters == "18"] = "AEP" # Alveolar epithelial cell

# Label with control and disease group
GSE155249_norm$group = ""
GSE155249_norm$group[GSE155249_norm$orig.ident == "8183"] = "Control"
GSE155249_norm$group[GSE155249_norm$orig.ident == "8184"] = "Control"
GSE155249_norm$group[GSE155249_norm$orig.ident == "0111"] = "Control"
GSE155249_norm$group[GSE155249_norm$orig.ident == "0112"] = "Control"
GSE155249_norm$group[GSE155249_norm$orig.ident == "8176"] = "COVID-19"
GSE155249_norm$group[GSE155249_norm$orig.ident == "8177"] = "COVID-19"
GSE155249_norm$group[GSE155249_norm$orig.ident == "8178"] = "COVID-19"
GSE155249_norm$group[GSE155249_norm$orig.ident == "8179"] = "COVID-19"
GSE155249_norm$group[GSE155249_norm$orig.ident == "8180"] = "COVID-19"
GSE155249_norm$group[GSE155249_norm$orig.ident == "8181"] = "COVID-19"
GSE155249_norm$group[GSE155249_norm$orig.ident == "8182"] = "COVID-19"
GSE155249_norm$group[GSE155249_norm$orig.ident == "0113"] = "COVID-19"
GSE155249_norm$group[GSE155249_norm$orig.ident == "0114"] = "COVID-19"
GSE155249_norm$group[GSE155249_norm$orig.ident == "0115"] = "COVID-19"
GSE155249_norm$group[GSE155249_norm$orig.ident == "0116"] = "COVID-19"
GSE155249_norm$group[GSE155249_norm$orig.ident == "0117"] = "COVID-19"
GSE155249_norm$group[GSE155249_norm$orig.ident == "0118"] = "COVID-19"
GSE155249_norm$group[GSE155249_norm$orig.ident == "0119"] = "COVID-19"
GSE155249_norm$group[GSE155249_norm$orig.ident == "0120"] = "COVID-19"

# This is the end of the Data Processing in GSE155249
#################################################################################################################