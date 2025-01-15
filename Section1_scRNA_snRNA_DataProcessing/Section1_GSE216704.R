# R version 4.2.2
# Load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(ggplot2)
library(ggpubr)

#################################################################################################################
# GSE216704 analysis
data1 = Read10X_h5(filename = "path/to/GSE216704/GSM6689175.h5")
data2 = Read10X_h5(filename = "path/to/GSE216704/GSM6689176.h5")
data3 = Read10X_h5(filename = "path/to/GSE216704/GSM6689177.h5")
data4 = Read10X_h5(filename = "path/to/GSE216704/GSM6689178.h5")
data5 = Read10X_h5(filename = "path/to/GSE216704/GSM6689179.h5")
data6 = Read10X_h5(filename = "path/to/GSE216704/GSM6689180.h5")
data7 = Read10X_h5(filename = "path/to/GSE216704/GSM6689181.h5")
data8 = Read10X_h5(filename = "path/to/GSE216704/GSM6689182.h5")
data9 = Read10X_h5(filename = "path/to/GSE216704/GSM6689183.h5")
data10 = Read10X_h5(filename = "path/to/GSE216704/GSM6689184.h5")
data11 = Read10X_h5(filename = "path/to/GSE216704/GSM6689186.h5")
data12 = Read10X_h5(filename = "path/to/GSE216704/GSM6689187.h5")
data13 = Read10X_h5(filename = "path/to/GSE216704/GSM6689188.h5")
data14 = Read10X_h5(filename = "path/to/GSE216704/GSM6689189.h5")
data15 = Read10X_h5(filename = "path/to/GSE216704/GSM6689190.h5")
data16 = Read10X_h5(filename = "path/to/GSE216704/GSM6689191.h5")
data17 = Read10X_h5(filename = "path/to/GSE216704/GSM6689192.h5")
data18 = Read10X_h5(filename = "path/to/GSE216704/GSM6689193.h5")

# Create Seurat object with the setting of min.cell = 3 and min.features = 200
data1 = CreateSeuratObject(counts = data1$`Gene Expression`, project = "9175", min.cells = 3, min.features = 200)
data2 = CreateSeuratObject(counts = data2$`Gene Expression`, project = "9176", min.cells = 3, min.features = 200)
data3 = CreateSeuratObject(counts = data3, project = "9177", min.cells = 3, min.features = 200)
data4 = CreateSeuratObject(counts = data4, project = "9178", min.cells = 3, min.features = 200)
data5 = CreateSeuratObject(counts = data5$`Gene Expression`, project = "9179", min.cells = 3, min.features = 200)
data6 = CreateSeuratObject(counts = data6$`Gene Expression`, project = "9180", min.cells = 3, min.features = 200)
data7 = CreateSeuratObject(counts = data7$`Gene Expression`, project = "9181", min.cells = 3, min.features = 200)
data8 = CreateSeuratObject(counts = data8$`Gene Expression`, project = "9182", min.cells = 3, min.features = 200)
data9 = CreateSeuratObject(counts = data9$`Gene Expression`, project = "9183", min.cells = 3, min.features = 200)
data10 = CreateSeuratObject(counts = data10$`Gene Expression`, project = "9184", min.cells = 3, min.features = 200)
data11 = CreateSeuratObject(counts = data11$`Gene Expression`, project = "9186", min.cells = 3, min.features = 200)
data12 = CreateSeuratObject(counts = data12$`Gene Expression`, project = "9187", min.cells = 3, min.features = 200)
data13 = CreateSeuratObject(counts = data13$`Gene Expression`, project = "9188", min.cells = 3, min.features = 200)
data14 = CreateSeuratObject(counts = data14$`Gene Expression`, project = "9189", min.cells = 3, min.features = 200)
data15 = CreateSeuratObject(counts = data15$`Gene Expression`, project = "9190", min.cells = 3, min.features = 200)
data16 = CreateSeuratObject(counts = data16$`Gene Expression`, project = "9191", min.cells = 3, min.features = 200)
data17 = CreateSeuratObject(counts = data17$`Gene Expression`, project = "9192", min.cells = 3, min.features = 200)
data18 = CreateSeuratObject(counts = data18$`Gene Expression`, project = "9193", min.cells = 3, min.features = 200)

# Merge all the seurat objects
GSE216704_object = merge(x = data1, y = c(data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14,data15,data16,data17,data18),
                         add.cell.ids = c("9175","9176","9177","9178","9179","9180","9181","9182","9183","9184","9186","9187","9188","9189","9190","9191","9192","9193"), project = "merge_data")

# Rename gene name from Ensembl to Gene Symbol
newname = read.csv("path/to/GSE216704/new_genename.csv", header = TRUE, row.names = 1) # input re-name gene list
newname = as.data.frame(newname)
colnames(newname) = "gene"
GSE216704_object@assays$RNA@counts@Dimnames[[1]] <- newname$gene
GSE216704_object@assays$RNA@data@Dimnames[[1]] <- newname$gene
GSE216704_object[["RNA"]]@meta.features <- data.frame(row.names = rownames(GSE216704_object[["RNA"]]))

# Identify the mitochondria-related genes in GSE216704_object
GSE216704_object[["percent.mt"]] <- PercentageFeatureSet(GSE216704_object, pattern = "^MT-")

# Visualization of nFeature_RNA, nCount_RNA and percent.mt for QC
VlnPlot(object = GSE216704_object, features = c("nCount_RNA","nFeature_RNA","percent.mt"))
GSE216704_object <- subset(GSE216704_object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)

# Normalization following the LogNormalize method
GSE216704_norm <- NormalizeData(GSE216704_object, normalization.method = "LogNormalize", scale.factor = 10000)

# Scaling the expression in GSE216704_norm
GSE216704_norm <- FindVariableFeatures(GSE216704_norm, selection.method = "vst", nfeatures = 2000)
GSE216704_norm <- ScaleData(GSE216704_norm, verbose = FALSE)


# PCA for dimensional reduction
GSE216704_norm <- RunPCA(GSE216704_norm)

# Find neighbors with top20 PCs and conduct clustering with resolution 0.5
GSE216704_norm <- FindNeighbors(GSE216704_norm, dims = 1:20, reduction = "pca")
GSE216704_norm <- FindClusters(GSE216704_norm, resolution = 0.5)

# Conduct UMAP and visualization
GSE216704_norm <- RunUMAP(GSE216704_norm, dims = 1:20)

# Due to the batch effect between samples, we decided to conduct integration
# Integration

data1 = subset(x = GSE216704_norm, subset = orig.ident == "9175")
data2 = subset(x = GSE216704_norm, subset = orig.ident == "9176")
data3 = subset(x = GSE216704_norm, subset = orig.ident == "9177")
data4 = subset(x = GSE216704_norm, subset = orig.ident == "9178")
data5 = subset(x = GSE216704_norm, subset = orig.ident == "9179")
data6 = subset(x = GSE216704_norm, subset = orig.ident == "9180")
data7 = subset(x = GSE216704_norm, subset = orig.ident == "9181")
data8 = subset(x = GSE216704_norm, subset = orig.ident == "9182")
data9 = subset(x = GSE216704_norm, subset = orig.ident == "9183")
data10 = subset(x = GSE216704_norm, subset = orig.ident == "9184")
data11 = subset(x = GSE216704_norm, subset = orig.ident == "9186")
data12 = subset(x = GSE216704_norm, subset = orig.ident == "9187")
data13 = subset(x = GSE216704_norm, subset = orig.ident == "9188")
data14 = subset(x = GSE216704_norm, subset = orig.ident == "9189")
data15 = subset(x = GSE216704_norm, subset = orig.ident == "9190")
data16 = subset(x = GSE216704_norm, subset = orig.ident == "9191")
data17 = subset(x = GSE216704_norm, subset = orig.ident == "9192")
data18 = subset(x = GSE216704_norm, subset = orig.ident == "9193")

ifnb.list <- list(data1 = data1, data2 = data2, data3 = data3, data4 = data4, data5 = data5, data6 = data6, data7 = data7, data8 = data8, data9 = data9, data10 = data10, data11 = data11, data12 = data12, data13 = data13, data14 = data14, data15 = data15, data16 = data16, data17 = data17, data18 = data18)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 4000)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
GSE216704_integrated <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(GSE216704_integrated) <- "integrated"

GSE216704_integrated <- ScaleData(GSE216704_integrated, verbose = FALSE)
GSE216704_integrated <- RunPCA(GSE216704_integrated, verbose = FALSE, assay = "integrated")
GSE216704_integrated <- FindNeighbors(GSE216704_integrated, dims = 1:20)
GSE216704_integrated <- FindClusters(GSE216704_integrated, resolution = 0.5)
GSE216704_integrated <- RunUMAP(GSE216704_integrated, dims = 1:20, verbose = FALSE)

# Find markers in each clusters, and save the results for manually cell annotation
GSE216704_markers = FindAllMarkers(object = GSE216704_integrated, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(x = GSE216704_markers, file = "path/to/YourDirectory/marker.csv")

# Manually cell annotation by searching cell markers
GSE216704_integrated$label = ""
GSE216704_integrated$label[GSE216704_integrated$integrated_snn_res.0.5 == "0"] = "Neutrophil_1"
GSE216704_integrated$label[GSE216704_integrated$integrated_snn_res.0.5 == "1"] = "Neutrophil_2"
GSE216704_integrated$label[GSE216704_integrated$integrated_snn_res.0.5 == "2"] = "Neutrophil_3"
GSE216704_integrated$label[GSE216704_integrated$integrated_snn_res.0.5 == "3"] = "PF-M" # Profibrotic macrophage
GSE216704_integrated$label[GSE216704_integrated$integrated_snn_res.0.5 == "4"] = "Tcell"
GSE216704_integrated$label[GSE216704_integrated$integrated_snn_res.0.5 == "5"] = "Macrophage_1"
GSE216704_integrated$label[GSE216704_integrated$integrated_snn_res.0.5 == "6"] = "CD14CD16-Mo" # CD14+CD16+ monocyte
GSE216704_integrated$label[GSE216704_integrated$integrated_snn_res.0.5 == "7"] = "CD8T" # CD8+ T cell
GSE216704_integrated$label[GSE216704_integrated$integrated_snn_res.0.5 == "8"] = "Neutrophil_4"
GSE216704_integrated$label[GSE216704_integrated$integrated_snn_res.0.5 == "9"] = "Pro-CD8T" # Proliferative CD8+ T cell
GSE216704_integrated$label[GSE216704_integrated$integrated_snn_res.0.5 == "10"] = "Neutrophil_5"
GSE216704_integrated$label[GSE216704_integrated$integrated_snn_res.0.5 == "11"] = "NK" # Nature killer cell
GSE216704_integrated$label[GSE216704_integrated$integrated_snn_res.0.5 == "12"] = "Macrophage_2"
GSE216704_integrated$label[GSE216704_integrated$integrated_snn_res.0.5 == "13"] = "pDC" # Plasmacytoid dendritic cell
GSE216704_integrated$label[GSE216704_integrated$integrated_snn_res.0.5 == "14"] = "CE" # Ciliated epithelial cell
GSE216704_integrated$label[GSE216704_integrated$integrated_snn_res.0.5 == "15"] = "DC" # Dendritic cell
GSE216704_integrated$label[GSE216704_integrated$integrated_snn_res.0.5 == "16"] = "PC" # Plasma cell

# Label with control and disease group
GSE216704_integrated$group = ""
GSE216704_integrated$group[GSE216704_integrated$orig.ident == "9175"] = "COVID-19"
GSE216704_integrated$group[GSE216704_integrated$orig.ident == "9176"] = "Non-COVID-19"
GSE216704_integrated$group[GSE216704_integrated$orig.ident == "9177"] = "Non-COVID-19"
GSE216704_integrated$group[GSE216704_integrated$orig.ident == "9178"] = "Non-COVID-19"
GSE216704_integrated$group[GSE216704_integrated$orig.ident == "9179"] = "COVID-19"
GSE216704_integrated$group[GSE216704_integrated$orig.ident == "9180"] = "COVID-19"
GSE216704_integrated$group[GSE216704_integrated$orig.ident == "9181"] = "Non-COVID-19"
GSE216704_integrated$group[GSE216704_integrated$orig.ident == "9182"] = "COVID-19"
GSE216704_integrated$group[GSE216704_integrated$orig.ident == "9183"] = "Control"
GSE216704_integrated$group[GSE216704_integrated$orig.ident == "9184"] = "Control"
GSE216704_integrated$group[GSE216704_integrated$orig.ident == "9186"] = "Non-COVID-19"
GSE216704_integrated$group[GSE216704_integrated$orig.ident == "9187"] = "Non-COVID-19"
GSE216704_integrated$group[GSE216704_integrated$orig.ident == "9188"] = "COVID-19"
GSE216704_integrated$group[GSE216704_integrated$orig.ident == "9189"] = "Non-COVID-19"
GSE216704_integrated$group[GSE216704_integrated$orig.ident == "9190"] = "COVID-19"
GSE216704_integrated$group[GSE216704_integrated$orig.ident == "9191"] = "Non-COVID-19"
GSE216704_integrated$group[GSE216704_integrated$orig.ident == "9192"] = "COVID-19"
GSE216704_integrated$group[GSE216704_integrated$orig.ident == "9193"] = "COVID-19"

# This is the end of the Data Processing in GSE216704
#################################################################################################################
