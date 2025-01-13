# R version 4.2.2
# Load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(ggplot2)
library(ggpubr)

#################################################################################################################
# GSE192483 analysis
# Load GSE192483 data
data1 <- Read10X(data.dir = "path/to/GSE192483/GSM5747736/", gene.column = 1)
data2 <- Read10X(data.dir = "path/to/GSE192483/GSM5747739/", gene.column = 1)
data3 <- Read10X(data.dir = "path/to/GSE192483/GSM5747741/", gene.column = 1)
data4 <- Read10X(data.dir = "path/to/GSE192483/GSM5747743/", gene.column = 1)
data5 <- Read10X(data.dir = "path/to/GSE192483/GSM5747745/", gene.column = 1)
data6 <- Read10X(data.dir = "path/to/GSE192483/GSM5747737/", gene.column = 1)
data7 <- Read10X(data.dir = "path/to/GSE192483/GSM5747738/", gene.column = 1)
data8 <- Read10X(data.dir = "path/to/GSE192483/GSM5747740/", gene.column = 1)
data9 <- Read10X(data.dir = "path/to/GSE192483/GSM5747742/", gene.column = 1)
data10 <- Read10X(data.dir = "path/to/GSE192483/GSM5747744/", gene.column = 1)
data11 <- Read10X(data.dir = "path/to/GSE192483/GSM5747746/", gene.column = 1)

# Create Seurat object with the setting of min.cell = 3 and min.features = 200
data1 <- CreateSeuratObject(counts = data1, project = "7736", min.cells = 3, min.features = 200)
data2 <- CreateSeuratObject(counts = data2, project = "7739", min.cells = 3, min.features = 200)
data3 <- CreateSeuratObject(counts = data3, project = "7741", min.cells = 3, min.features = 200)
data4 <- CreateSeuratObject(counts = data4, project = "7743", min.cells = 3, min.features = 200)
data5 <- CreateSeuratObject(counts = data5, project = "7745", min.cells = 3, min.features = 200)
data6 <- CreateSeuratObject(counts = data6, project = "7737", min.cells = 3, min.features = 200)
data7 <- CreateSeuratObject(counts = data7, project = "7738", min.cells = 3, min.features = 200)
data8 <- CreateSeuratObject(counts = data8, project = "7740", min.cells = 3, min.features = 200)
data9 <- CreateSeuratObject(counts = data9, project = "7742", min.cells = 3, min.features = 200)
data10 <- CreateSeuratObject(counts = data10, project = "7744", min.cells = 3, min.features = 200)
data11 <- CreateSeuratObject(counts = data11, project = "7746", min.cells = 3, min.features = 200)

# Merge all the seurat objects
GSE192483_object = merge(x = data1, y = c(data2,data3,data4,data5,data6,data7,data8,data9,data10,data11),
                         add.cell.ids = c("7736","7739","7741","7743","7745","7737","7738","7740","7742","7744","7746"), project = "TB")

# Rename gene name from Ensembl to Gene Symbol
newname = read.csv("path/to/GSE192483/new_name.csv", header = TRUE, row.names = 1)
newname = as.data.frame(newname)
colnames(newname) = "gene"
GSE192483_object@assays$RNA@counts@Dimnames[[1]] <- newname$gene
GSE192483_object@assays$RNA@data@Dimnames[[1]] <- newname$gene
GSE192483_object[["RNA"]]@meta.features <- data.frame(row.names = rownames(GSE192483_object[["RNA"]]))

# Identify the mitochondria-related genes in GSE192483_object
GSE192483_object[["percent.mt"]] <- PercentageFeatureSet(GSE192483_object, pattern = "^MT-")

# Visualization of nFeature_RNA, nCount_RNA and percent.mt for QC
VlnPlot(object = GSE192483_object, features = c("nCount_RNA","nFeature_RNA","percent.mt"))
GSE192483_object <- subset(GSE192483_object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)

# Normalization following the LogNormalize method
GSE192483_norm <- NormalizeData(GSE192483_object, normalization.method = "LogNormalize", scale.factor = 10000)

# Scaling the expression with all genes in GSE192483_norm
GSE192483_norm <- FindVariableFeatures(GSE192483_norm, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GSE192483_norm)
GSE192483_norm <- ScaleData(GSE192483_norm, features = all.genes)

# PCA for dimensional reduction with the highly variable features in GSE192483_norm
GSE192483_norm <- RunPCA(GSE192483_norm, features = VariableFeatures(object = GSE192483_norm))

# Find neighbors with top20 PCs and conduct clustering with resolution 0.5
GSE192483_norm <- FindNeighbors(GSE192483_norm, dims = 1:20)
GSE192483_norm <- FindClusters(GSE192483_norm, resolution = 0.5)

# Conduct UMAP and visualization
GSE192483_norm <- RunUMAP(GSE192483_norm, dims = 1:20, seed.use = 1)
DimPlot(GSE192483_norm, reduction = "umap", shuffle = TRUE)

# Find markers in each clusters, and save the results for manually cell annotation
GSE192483_markers = FindAllMarkers(object = GSE192483_norm, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(x = GSE192483_markers, file = "path/to/YourDirectory/marker.csv")

# Manually cell annotation by searching cell markers
GSE192483_norm$labe = ""
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "0"] = "CD4T" # CD4+ T cell
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "1"] = "CD8T" # CD8+ T cell
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "2"] = "PC" # Plasma cell
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "3"] = "PF-M" # Profibrotic macrophage
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "4"] = "NK" # Nature killer cell
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "5"] = "NK" # Nature killer cell
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "6"] = "TR-M" # Tissue-resident macrophage
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "7"] = "CD14CD16-Mo" # CD14+CD16+ monocyte
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "8"] = "ExT" # Exhaustive T cell
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "9"] = "Bcell"
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "10"] = "TR-M" # Tissue-resident macrophage
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "11"] = "PC" # Plasma cell
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "12"] = "AEP" # Avleolar epithelial cell
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "13"] = "Cycling" # Cycling cell
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "14"] = "DC" # Dendritic cell
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "15"] = "BA" # Basophil
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "16"] = "pDC" # Plasmacytoid dendritic cell
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "17"] = "Cycling" # Cycling cell
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "18"] = "Myo-FB" # Myofibroblast
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "19"] = "PC" # Plasma cell
GSE192483_norm$label[GSE192483_norm$seurat_clusters == "20"] = "CE" # Ciliated epithelial cell

# Label with control and disease group
GSE192483_norm$group = ""
GSE192483_norm$group[GSE192483_norm$orig.ident == "7736"] = "Control"
GSE192483_norm$group[GSE192483_norm$orig.ident == "7739"] = "Control"
GSE192483_norm$group[GSE192483_norm$orig.ident == "7741"] = "Control"
GSE192483_norm$group[GSE192483_norm$orig.ident == "7743"] = "Control"
GSE192483_norm$group[GSE192483_norm$orig.ident == "7745"] = "Control"
GSE192483_norm$group[GSE192483_norm$orig.ident == "7737"] = "TB"
GSE192483_norm$group[GSE192483_norm$orig.ident == "7738"] = "TB"
GSE192483_norm$group[GSE192483_norm$orig.ident == "7740"] = "TB"
GSE192483_norm$group[GSE192483_norm$orig.ident == "7742"] = "TB"
GSE192483_norm$group[GSE192483_norm$orig.ident == "7744"] = "TB"
GSE192483_norm$group[GSE192483_norm$orig.ident == "7746"] = "TB"

# This is the end of the Data Processing in GSE192483
#################################################################################################################