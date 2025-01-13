# R version 4.2.2
# Load packages
library(Seurat)
library(ggplot2)
library(SeuratObject)
library(patchwork)
library(dplyr)
library(CellChat)
library(ggpubr)
library(scales)

#################################################################################################################
# Spatial transcriptomics from Mendeley Data (DOI:10.17632/xjtv62ncwr.1)
# Load C1 data
C1_img = Read10X_Image("path/to/V10F24-110_C1/spatial", image.name = "tissue_lowres_image.png")
C1_data = Load10X_Spatial(data.dir = "path/to/V10F24-110_C1/", filename = "filtered_feature_bc_matrix.h5", image = C1_img)

# Rename gene name from Ensembl to Gene Symbol
newname = read.csv("path/to/V10F24-110_C1/new_name.csv", header = TRUE, row.names = 1)
newname = as.data.frame(newname)
colnames(newname) = "gene"
C1_data@assays$Spatial@counts@Dimnames[[1]] <- newname$gene
C1_data@assays$Spatial@data@Dimnames[[1]] <- newname$gene
C1_data[["Spatial"]]@meta.features <- data.frame(row.names = rownames(C1_data[["Spatial"]]))

# Normalization following the SCTransform method
C1_norm <- SCTransform(C1_data, assay = "Spatial", verbose = FALSE)

# PCA for dimensional reduction with the highly variable features in GSE155249_norm
C1_norm <- RunPCA(C1_norm, assay = "SCT", verbose = FALSE)

# Find neighbors with top20 PCs and conduct clustering with resolution 0.5
C1_norm <- FindNeighbors(C1_norm, dims = 1:20)
C1_norm <- FindClusters(C1_norm, verbose = FALSE, resolution = 0.5)

# Conduct UMAP and visualization
C1_norm <- RunUMAP(C1_norm, dims = 1:20)
SpatialFeaturePlot(object = C1_norm, features = "GAS6", pt.size.factor = 6, stroke = 1, image.alpha = 0) + 
  scale_fill_continuous(low = "transparent", high = "red")

# Find markers in each clusters, and save the results for manually cell annotation
C1_marker = FindAllMarkers(object = C1_norm, logfc.threshold = 0, only.pos = FALSE)
write.csv(x = C1_marker, file = "path/to/YourDirectory/marker.csv")

# Manually cell annotation by searching cell markers
C1_norm$label[C1_norm$seurat_clusters == "0"] = "AEP" # Alveolar epithelial cell-enriched region
C1_norm$label[C1_norm$seurat_clusters == "1"] = "Immune" # Immune cell-enriched region
C1_norm$label[C1_norm$seurat_clusters == "2"] = "AEP" # Alveolar epithelial cell-enriched region
C1_norm$label[C1_norm$seurat_clusters == "3"] = "CE" # Ciliated epithelial cell-enriched region

# This is the end of the Data Processing in C1 data
#################################################################################################################
# Spatial transcriptomics from Mendeley Data (DOI:10.17632/xjtv62ncwr.1)
# Load D1 data
D1_img = Read10X_Image("path/to/V10F24-110_D1/spatial", image.name = "tissue_hires_image.png")
D1_data = Load10X_Spatial(data.dir = "path/to/V10F24-110_D1/", filename = "filtered_feature_bc_matrix.h5", image = D1_img)

# Rename gene name from Ensembl to Gene Symbol
newname = read.csv("path/to/V10F24-110_D1/new_name.csv", header = TRUE, row.names = 1)
newname = as.data.frame(newname)
colnames(newname) = "gene"
D1_data@assays$Spatial@counts@Dimnames[[1]] <- newname$gene
D1_data@assays$Spatial@data@Dimnames[[1]] <- newname$gene
D1_data[["Spatial"]]@meta.features <- data.frame(row.names = rownames(D1_data[["Spatial"]]))

# Normalization following the SCTransform method
D1_norm <- SCTransform(D1_data, assay = "Spatial", verbose = FALSE)

# PCA for dimensional reduction with the highly variable features in GSE155249_norm
D1_norm <- RunPCA(D1_norm, assay = "SCT", verbose = FALSE)

# Find neighbors with top20 PCs and conduct clustering with resolution 0.5
D1_norm <- FindNeighbors(D1_norm, dims = 1:20)
D1_norm <- FindClusters(D1_norm, verbose = FALSE, resolution = 0.5)

# Conduct UMAP and visualization
D1_norm <- RunUMAP(D1_norm, dims = 1:20)
SpatialFeaturePlot(object = D1_norm, features = "GAS6", pt.size.factor = 6, stroke = 1, image.alpha = 0) + 
  scale_fill_continuous(low = "transparent", high = "red")

# Find markers in each clusters, and save the results for manually cell annotation
D1_marker = FindAllMarkers(object = D1_norm, logfc.threshold = 0, only.pos = FALSE)
write.csv(x = D1_marker, file = "path/to/YourDirectory/marker.csv")

# Manually cell annotation by searching cell markers
D1_norm$label[D1_norm$seurat_clusters == "0"] = "AEP" # Alveolar epithelial cell-enriched region
D1_norm$label[D1_norm$seurat_clusters == "1"] = "AEP" # Alveolar epithelial cell-enriched region
D1_norm$label[D1_norm$seurat_clusters == "2"] = "Immune" # Immune cell-enriched region
D1_norm$label[D1_norm$seurat_clusters == "3"] = "CE" # Ciliated epithelial cell-enriched region

# This is the end of the Data Processing in D1 data
#################################################################################################################