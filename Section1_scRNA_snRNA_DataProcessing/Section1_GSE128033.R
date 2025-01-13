# R version 4.2.2
# Load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(ggplot2)
library(ggpubr)

#################################################################################################################
# GSE128033 analysis
# Load GSE128033 data
data1 = Read10X(data.dir = "path/to/GSE128033/GSM3660641/", gene.column = 1)
data2 = Read10X(data.dir = "path/to/GSE128033/GSM3660642/", gene.column = 1)
data3 = Read10X(data.dir = "path/to/GSE128033/GSM3660643/", gene.column = 1)
data4 = Read10X(data.dir = "path/to/GSE128033/GSM3660644/", gene.column = 1)
data5 = Read10X(data.dir = "path/to/GSE128033/GSM3660645/", gene.column = 1)
data6 = Read10X(data.dir = "path/to/GSE128033/GSM3660646/", gene.column = 1)
data7 = Read10X(data.dir = "path/to/GSE128033/GSM3660647/", gene.column = 1)
data8 = Read10X(data.dir = "path/to/GSE128033/GSM3660648/", gene.column = 1)
data9 = Read10X(data.dir = "path/to/GSE128033/GSM3660649/", gene.column = 1)
data10 = Read10X(data.dir = "path/to/GSE128033/GSM3660650/", gene.column = 1)
data11 = Read10X(data.dir = "path/to/GSE128033/GSM3660651/", gene.column = 1)
data12 = Read10X(data.dir = "path/to/GSE128033/GSM3660652/", gene.column = 1)
data13 = Read10X(data.dir = "path/to/GSE128033/GSM3660653/", gene.column = 1)
data14 = Read10X(data.dir = "path/to/GSE128033/GSM3660654/", gene.column = 1)
data15 = Read10X(data.dir = "path/to/GSE128033/GSM3660655/", gene.column = 1)
data16 = Read10X(data.dir = "path/to/GSE128033/GSM3660656/", gene.column = 1)
data17 = Read10X(data.dir = "path/to/GSE128033/GSM3660657/", gene.column = 1)
data18 = Read10X(data.dir = "path/to/GSE128033/GSM3660658/", gene.column = 1)

# Create Seurat object with the setting of min.cell = 3 and min.features = 200
data1 <- CreateSeuratObject(counts = data1, project = "0641", min.cells = 3, min.features = 200)
data2 <- CreateSeuratObject(counts = data2, project = "0642", min.cells = 3, min.features = 200)
data3 <- CreateSeuratObject(counts = data3, project = "0643", min.cells = 3, min.features = 200)
data4 <- CreateSeuratObject(counts = data4, project = "0644", min.cells = 3, min.features = 200)
data5 <- CreateSeuratObject(counts = data5, project = "0645", min.cells = 3, min.features = 200)
data6 <- CreateSeuratObject(counts = data6, project = "0646", min.cells = 3, min.features = 200)
data7 <- CreateSeuratObject(counts = data7, project = "0647", min.cells = 3, min.features = 200)
data8 <- CreateSeuratObject(counts = data8, project = "0648", min.cells = 3, min.features = 200)
data9 <- CreateSeuratObject(counts = data9, project = "0649", min.cells = 3, min.features = 200)
data10 <- CreateSeuratObject(counts = data10, project = "0650", min.cells = 3, min.features = 200)
data11 <- CreateSeuratObject(counts = data11, project = "0651", min.cells = 3, min.features = 200)
data12 <- CreateSeuratObject(counts = data12, project = "0652", min.cells = 3, min.features = 200)
data13 <- CreateSeuratObject(counts = data13, project = "0653", min.cells = 3, min.features = 200)
data14 <- CreateSeuratObject(counts = data14, project = "0654", min.cells = 3, min.features = 200)
data15 <- CreateSeuratObject(counts = data15, project = "0655", min.cells = 3, min.features = 200)
data16 <- CreateSeuratObject(counts = data16, project = "0656", min.cells = 3, min.features = 200)
data17 <- CreateSeuratObject(counts = data17, project = "0657", min.cells = 3, min.features = 200)
data18 <- CreateSeuratObject(counts = data18, project = "0658", min.cells = 3, min.features = 200)

# Merge all the seurat objects
GSE128033_object = merge(x = data1, y = c(data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14,data15,data16,data17,data18),
                         add.cell.ids = c("0641","0642","0643","0644","0645","0646","0647","0648","0649","0650","0651","0652","0653","0654","0655","0656","0657","0658"), project = "IPF")

# Rename gene name from Ensembl to Gene Symbol
newname = read.csv("path/to/GSE128033/new_name.csv", header = TRUE, row.names = 1)
newname = as.data.frame(newname)
colnames(newname) = "gene"
GSE128033_object@assays$RNA@counts@Dimnames[[1]] <- newname$gene
GSE128033_object@assays$RNA@data@Dimnames[[1]] <- newname$gene
GSE128033_object[["RNA"]]@meta.features <- data.frame(row.names = rownames(GSE128033_object[["RNA"]]))

# Identify the mitochondria-related genes in GSE128033_object
GSE128033_object[["percent.mt"]] <- PercentageFeatureSet(GSE128033_object, pattern = "^MT-")

# Visualization of nFeature_RNA, nCount_RNA and percent.mt for QC
VlnPlot(object = GSE128033_object, features = c("nCount_RNA","nFeature_RNA","percent.mt"))
GSE128033_object <- subset(GSE128033_object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)

# Normalization following the LogNormalize method
GSE128033_norm <- NormalizeData(GSE128033_object, normalization.method = "LogNormalize", scale.factor = 10000)

# Scaling the expression with all genes in GSE128033_norm
GSE128033_norm <- FindVariableFeatures(GSE128033_norm, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GSE128033_norm)
GSE128033_norm <- ScaleData(GSE128033_norm, features = all.genes)

# PCA for dimensional reduction with the highly variable features in GSE128033_norm
GSE128033_norm <- RunPCA(GSE128033_norm, features = VariableFeatures(object = GSE128033_norm))

# Find neighbors with top20 PCs and conduct clustering with resolution 0.5
GSE128033_norm <- FindNeighbors(GSE128033_norm, dims = 1:20)
GSE128033_norm <- FindClusters(GSE128033_norm, resolution = 0.5)

# Conduct UMAP and visualization
GSE128033_norm <- RunUMAP(GSE128033_norm, dims = 1:20, seed.use = 1)
DimPlot(GSE128033_norm, reduction = "umap", shuffle = TRUE)

# Due to the batch effect between samples, we decided to conduct integration
# Integration
GSE128033_norm$batch = "B"
GSE128033_norm$batch[GSE128033_norm$orig.ident == "0641"] = "A"
GSE128033_norm$batch[GSE128033_norm$orig.ident == "0642"] = "A"
GSE128033_norm$batch[GSE128033_norm$orig.ident == "0643"] = "A"

batch_A = subset(x = GSE128033_norm, subset = batch == "A")
batch_B = subset(x = GSE128033_norm, subset = batch == "B")

ifnb.list <- list(data1 = batch_A, data2 = batch_B)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 4000)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
GSE128033_integrated <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(GSE128033_integrated) <- "integrated"

GSE128033_integrated <- ScaleData(GSE128033_integrated, verbose = FALSE)
GSE128033_integrated <- RunPCA(GSE128033_integrated, verbose = FALSE, assay = "integrated")
GSE128033_integrated <- FindNeighbors(GSE128033_integrated, dims = 1:20)
GSE128033_integrated <- FindClusters(GSE128033_integrated, resolution = 0.5)
GSE128033_integrated <- RunUMAP(GSE128033_integrated, dims = 1:20, verbose = FALSE)

# Find markers in each clusters, and save the results for manually cell annotation
GSE128033_markers = FindAllMarkers(object = GSE128033_integrated, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(x = GSE128033_markers, file = "path/to/YourDirectory/marker.csv")

# Manually cell annotation by searching cell markers
GSE128033_integrated$label = ""
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "0"] = "PF-M" # Profibrotic macrophage
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "1"] = "TR-M" # Tissue-resident macrophage
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "2"] = "Tcell"
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "3"] = "CD8T" # CD8+ T cell
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "4"] = "TR-M" # Tissue-resident macrophage
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "5"] = "CD14CD16-Mo" # CD14+CD16+ monocyte
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "6"] = "Myo-FB" # Myofibroblast
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "7"] = "ED" # Endothelial cell
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "8"] = "ED" # Endothelial cell
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "9"] = "SMC" # Smooth muscle cell
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "10"] = "AEP" # Alveolar epithelial cell
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "11"] = "AEP" # Avleolar epithelial cell
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "12"] = "CE" # Ciliated epithelial cell
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "13"] = "Mast" # Mast cell
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "14"] = "AEP" # Alveolar epithelial cell
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "15"] = "PC" # Plasma cell
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "16"] = "L-ED" # Lymphatic endothelial cell
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "17"] = "DC" # Dendritic cell
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "18"] = "CD8T" # CD8+ T cell
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "19"] = "CE" # Ciliated cell
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "20"] = "Cycling" # Cycling cell
GSE128033_integrated$label[GSE128033_integrated$integrated_snn_res.0.5 == "21"] = "TR-M" # Tissue-resident cell

# Label with control and disease group
GSE128033_integrated$group = ""
GSE128033_integrated$group[GSE128033_integrated$orig.ident == "0641"] = "Control"
GSE128033_integrated$group[GSE128033_integrated$orig.ident == "0642"] = "Control"
GSE128033_integrated$group[GSE128033_integrated$orig.ident == "0643"] = "Control"
GSE128033_integrated$group[GSE128033_integrated$orig.ident == "0644"] = "Control"
GSE128033_integrated$group[GSE128033_integrated$orig.ident == "0645"] = "Control"
GSE128033_integrated$group[GSE128033_integrated$orig.ident == "0646"] = "Control"
GSE128033_integrated$group[GSE128033_integrated$orig.ident == "0647"] = "Control"
GSE128033_integrated$group[GSE128033_integrated$orig.ident == "0648"] = "Control"
GSE128033_integrated$group[GSE128033_integrated$orig.ident == "0649"] = "Control"
GSE128033_integrated$group[GSE128033_integrated$orig.ident == "0650"] = "Control"

GSE128033_integrated$group[GSE128033_integrated$orig.ident == "0651"] = "IPF"
GSE128033_integrated$group[GSE128033_integrated$orig.ident == "0652"] = "IPF"
GSE128033_integrated$group[GSE128033_integrated$orig.ident == "0653"] = "IPF"
GSE128033_integrated$group[GSE128033_integrated$orig.ident == "0654"] = "IPF"
GSE128033_integrated$group[GSE128033_integrated$orig.ident == "0655"] = "IPF"
GSE128033_integrated$group[GSE128033_integrated$orig.ident == "0656"] = "IPF"
GSE128033_integrated$group[GSE128033_integrated$orig.ident == "0657"] = "IPF"
GSE128033_integrated$group[GSE128033_integrated$orig.ident == "0658"] = "IPF"

# This is the end of the Data Processing in GSE128033
#################################################################################################################