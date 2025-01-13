# R version 4.2.2
# Load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(ggplot2)
library(ggpubr)

#################################################################################################################
# GSE188172 analysis
data1 = read.table("path/to/GSE188172/GSM5671235_CF1_all_counts.txt",sep="\t",header=TRUE)
data2 = read.table("path/to/GSE188172/GSM5671236_CF2_all_counts.txt",sep="\t",header=TRUE)
data3 = read.table("path/to/GSE188172/GSM5671237_CF3_all_counts.txt",sep="\t",header=TRUE)
data4 = read.table("path/to/GSE188172/GSM5671238_CF4_all_counts.txt",sep="\t",header=TRUE)
data5 = read.table("path/to/GSE188172/GSM5671239_CF5_all_counts.txt",sep="\t",header=TRUE)
data6 = read.table("path/to/GSE188172/GSM5671240_CF6_all_counts.txt",sep="\t",header=TRUE)
data7 = read.table("path/to/GSE188172/GSM5671241_CF7_all_counts.txt",sep="\t",header=TRUE)
data8 = read.table("path/to/GSE188172/GSM5671242_CM1_all_counts.txt",sep="\t",header=TRUE)
data9 = read.table("path/to/GSE188172/GSM5671243_CM2_all_counts.txt",sep="\t",header=TRUE)
data10 = read.table("path/to/GSE188172/GSM5671244_CM3_all_counts.txt",sep="\t",header=TRUE)
data11 = read.table("path/to/GSE188172/GSM5671245_CM4_all_counts.txt",sep="\t",header=TRUE)
data12 = read.table("path/to/GSE188172/GSM5671246_CM5_all_counts.txt",sep="\t",header=TRUE)

data13 = read.table("path/to/GSE188172/GSM5671304_PMF1_all_counts.txt",sep="\t",header=TRUE)
data14 = read.table("path/to/GSE188172/GSM5671305_PMF2_all_counts.txt",sep="\t",header=TRUE)
data15 = read.table("path/to/GSE188172/GSM5671306_PMF3_all_counts.txt",sep="\t",header=TRUE)
data16 = read.table("path/to/GSE188172/GSM5671307_PMF4_all_counts.txt",sep="\t",header=TRUE)
data17 = read.table("path/to/GSE188172/GSM5671308_PMF5_all_counts.txt",sep="\t",header=TRUE)
data18 = read.table("path/to/GSE188172/GSM5671309_PMF6_all_counts.txt",sep="\t",header=TRUE)
data19 = read.table("path/to/GSE188172/GSM5671310_PMF7_all_counts.txt",sep="\t",header=TRUE)
data20 = read.table("path/to/GSE188172/GSM5671311_PMF8_all_counts.txt",sep="\t",header=TRUE)
data21 = read.table("path/to/GSE188172/GSM5671312_PMM1_all_counts.txt",sep="\t",header=TRUE)
data22 = read.table("path/to/GSE188172/GSM5671313_PMM2_all_counts.txt",sep="\t",header=TRUE)
data23 = read.table("path/to/GSE188172/GSM5671314_PMM3_all_counts.txt",sep="\t",header=TRUE)
data24 = read.table("path/to/GSE188172/GSM5671315_PMM4_all_counts.txt",sep="\t",header=TRUE)
data25 = read.table("path/to/GSE188172/GSM5671316_PMM5_all_counts.txt",sep="\t",header=TRUE)
data26 = read.table("path/to/GSE188172/GSM5671317_PMM6_all_counts.txt",sep="\t",header=TRUE)
data27 = read.table("path/to/GSE188172/GSM5671318_PMM7_all_counts.txt",sep="\t",header=TRUE)
data28 = read.table("path/to/GSE188172/GSM5671319_PMM8_all_counts.txt",sep="\t",header=TRUE)

data29 = read.table("path/to/GSE188172/GSM5671320_PSF1a_all_counts.txt",sep="\t",header=TRUE)
data30 = read.table("path/to/GSE188172/GSM5671321_PSF1b_all_counts.txt",sep="\t",header=TRUE)
data31 = read.table("path/to/GSE188172/GSM5671322_PSF1c_all_counts.txt",sep="\t",header=TRUE)
data32 = read.table("path/to/GSE188172/GSM5671323_PSF2a_all_counts.txt",sep="\t",header=TRUE)
data33 = read.table("path/to/GSE188172/GSM5671324_PSF2b_all_counts.txt",sep="\t",header=TRUE)
data34 = read.table("path/to/GSE188172/GSM5671325_PSF2c_all_counts.txt",sep="\t",header=TRUE)
data35 = read.table("path/to/GSE188172/GSM5671326_PSF3a_all_counts.txt",sep="\t",header=TRUE)
data36 = read.table("path/to/GSE188172/GSM5671327_PSF3b_all_counts.txt",sep="\t",header=TRUE)
data37 = read.table("path/to/GSE188172/GSM5671328_PSF4a_all_counts.txt",sep="\t",header=TRUE)
data38 = read.table("path/to/GSE188172/GSM5671329_PSF4b_all_counts.txt",sep="\t",header=TRUE)
data39 = read.table("path/to/GSE188172/GSM5671330_PSF5a_all_counts.txt",sep="\t",header=TRUE)
data40 = read.table("path/to/GSE188172/GSM5671331_PSF5b_all_counts.txt",sep="\t",header=TRUE)
data41 = read.table("path/to/GSE188172/GSM5671332_PSF6a_all_counts.txt",sep="\t",header=TRUE)
data42 = read.table("path/to/GSE188172/GSM5671333_PSF6b_all_counts.txt",sep="\t",header=TRUE)
data43 = read.table("path/to/GSE188172/GSM5671334_PSM1a_all_counts.txt",sep="\t",header=TRUE)
data44 = read.table("path/to/GSE188172/GSM5671335_PSM1b_all_counts.txt",sep="\t",header=TRUE)
data45 = read.table("path/to/GSE188172/GSM5671336_PSM1c_all_counts.txt",sep="\t",header=TRUE)
data46 = read.table("path/to/GSE188172/GSM5671337_PSM2a_all_counts.txt",sep="\t",header=TRUE)
data47 = read.table("path/to/GSE188172/GSM5671338_PSM2b_all_counts.txt",sep="\t",header=TRUE)
data48 = read.table("path/to/GSE188172/GSM5671339_PSM2c_all_counts.txt",sep="\t",header=TRUE)
data49 = read.table("path/to/GSE188172/GSM5671340_PSM3a_all_counts.txt",sep="\t",header=TRUE)
data50 = read.table("path/to/GSE188172/GSM5671341_PSM3b_all_counts.txt",sep="\t",header=TRUE)
data51 = read.table("path/to/GSE188172/GSM5671342_PSM3c_all_counts.txt",sep="\t",header=TRUE)
data52 = read.table("path/to/GSE188172/GSM5671343_PSM4a_all_counts.txt",sep="\t",header=TRUE)
data53 = read.table("path/to/GSE188172/GSM5671344_PSM4b_all_counts.txt",sep="\t",header=TRUE)
data54 = read.table("path/to/GSE188172/GSM5671345_PSM4c_all_counts.txt",sep="\t",header=TRUE)
data55 = read.table("path/to/GSE188172/GSM5671346_PSM5a_all_counts.txt",sep="\t",header=TRUE)
data56 = read.table("path/to/GSE188172/GSM5671347_PSM5b_all_counts.txt",sep="\t",header=TRUE)
data57 = read.table("path/to/GSE188172/GSM5671348_PSM6a_all_counts.txt",sep="\t",header=TRUE)
data58 = read.table("path/to/GSE188172/GSM5671349_PSM6b_all_counts.txt",sep="\t",header=TRUE)
data59 = read.table("path/to/GSE188172/GSM5671350_PSM6c_all_counts.txt",sep="\t",header=TRUE)

# Create Seurat object with the setting of min.cell = 3 and min.features = 200
data1 = CreateSeuratObject(counts = data1, project = "CF1", min.cells = 3, min.features = 200, names.delim = "-")
data2 = CreateSeuratObject(counts = data2, project = "CF2", min.cells = 3, min.features = 200, names.delim = "-")
data3 = CreateSeuratObject(counts = data3, project = "CF3", min.cells = 3, min.features = 200, names.delim = "-")
data4 = CreateSeuratObject(counts = data4, project = "CF4", min.cells = 3, min.features = 200, names.delim = "-")
data5 = CreateSeuratObject(counts = data5, project = "CF5", min.cells = 3, min.features = 200, names.delim = "-")
data6 = CreateSeuratObject(counts = data6, project = "CF6", min.cells = 3, min.features = 200, names.delim = "-")
data7 = CreateSeuratObject(counts = data7, project = "CF7", min.cells = 3, min.features = 200, names.delim = "-")
data8 = CreateSeuratObject(counts = data8, project = "CM1", min.cells = 3, min.features = 200, names.delim = "-")
data9 = CreateSeuratObject(counts = data9, project = "CM2", min.cells = 3, min.features = 200, names.delim = "-")
data10 = CreateSeuratObject(counts = data10, project = "CM3", min.cells = 3, min.features = 200, names.delim = "-")
data11 = CreateSeuratObject(counts = data11, project = "CM4", min.cells = 3, min.features = 200, names.delim = "-")
data12 = CreateSeuratObject(counts = data12, project = "CM5", min.cells = 3, min.features = 200, names.delim = "-")
data13 = CreateSeuratObject(counts = data13, project = "PMF1", min.cells = 3, min.features = 200, names.delim = "-")
data14 = CreateSeuratObject(counts = data14, project = "PMF2", min.cells = 3, min.features = 200, names.delim = "-")
data15 = CreateSeuratObject(counts = data15, project = "PMF3", min.cells = 3, min.features = 200, names.delim = "-")
data16 = CreateSeuratObject(counts = data16, project = "PMF4", min.cells = 3, min.features = 200, names.delim = "-")
data17 = CreateSeuratObject(counts = data17, project = "PMF5", min.cells = 3, min.features = 200, names.delim = "-")
data18 = CreateSeuratObject(counts = data18, project = "PMF6", min.cells = 3, min.features = 200, names.delim = "-")
data19 = CreateSeuratObject(counts = data19, project = "PMF7", min.cells = 3, min.features = 200, names.delim = "-")
data20 = CreateSeuratObject(counts = data20, project = "PMF8", min.cells = 3, min.features = 200, names.delim = "-")
data21 = CreateSeuratObject(counts = data21, project = "PMM1", min.cells = 3, min.features = 200, names.delim = "-")
data22 = CreateSeuratObject(counts = data22, project = "PMM2", min.cells = 3, min.features = 200, names.delim = "-")
data23 = CreateSeuratObject(counts = data23, project = "PMM3", min.cells = 3, min.features = 200, names.delim = "-")
data24 = CreateSeuratObject(counts = data24, project = "PMM4", min.cells = 3, min.features = 200, names.delim = "-")
data25 = CreateSeuratObject(counts = data25, project = "PMM5", min.cells = 3, min.features = 200, names.delim = "-")
data26 = CreateSeuratObject(counts = data26, project = "PMM6", min.cells = 3, min.features = 200, names.delim = "-")
data27 = CreateSeuratObject(counts = data27, project = "PMM7", min.cells = 3, min.features = 200, names.delim = "-")
data28 = CreateSeuratObject(counts = data28, project = "PMM8", min.cells = 3, min.features = 200, names.delim = "-")
data29 = CreateSeuratObject(counts = data29, project = "PSF1a", min.cells = 3, min.features = 200, names.delim = "-")
data30 = CreateSeuratObject(counts = data30, project = "PSF1b", min.cells = 3, min.features = 200, names.delim = "-")
data31 = CreateSeuratObject(counts = data31, project = "PSF1c", min.cells = 3, min.features = 200, names.delim = "-")
data32 = CreateSeuratObject(counts = data32, project = "PSF2a", min.cells = 3, min.features = 200, names.delim = "-")
data33 = CreateSeuratObject(counts = data33, project = "PSF2b", min.cells = 3, min.features = 200, names.delim = "-")
data34 = CreateSeuratObject(counts = data34, project = "PSF2c", min.cells = 3, min.features = 200, names.delim = "-")
data35 = CreateSeuratObject(counts = data35, project = "PSF3a", min.cells = 3, min.features = 200, names.delim = "-")
data36 = CreateSeuratObject(counts = data36, project = "PSF3b", min.cells = 3, min.features = 200, names.delim = "-")
data37 = CreateSeuratObject(counts = data37, project = "PSF4a", min.cells = 3, min.features = 200, names.delim = "-")
data38 = CreateSeuratObject(counts = data38, project = "PSF4b", min.cells = 3, min.features = 200, names.delim = "-")
data39 = CreateSeuratObject(counts = data39, project = "PSF5a", min.cells = 3, min.features = 200, names.delim = "-")
data40 = CreateSeuratObject(counts = data40, project = "PSF5b", min.cells = 3, min.features = 200, names.delim = "-")
data41 = CreateSeuratObject(counts = data41, project = "PSF6a", min.cells = 3, min.features = 200, names.delim = "-")
data42 = CreateSeuratObject(counts = data42, project = "PSF6b", min.cells = 3, min.features = 200, names.delim = "-")
data43 = CreateSeuratObject(counts = data43, project = "PSM1a", min.cells = 3, min.features = 200, names.delim = "-")
data44 = CreateSeuratObject(counts = data44, project = "PSM1b", min.cells = 3, min.features = 200, names.delim = "-")
data45 = CreateSeuratObject(counts = data45, project = "PSM1c", min.cells = 3, min.features = 200, names.delim = "-")
data46 = CreateSeuratObject(counts = data46, project = "PSM2a", min.cells = 3, min.features = 200, names.delim = "-")
data47 = CreateSeuratObject(counts = data47, project = "PSM2b", min.cells = 3, min.features = 200, names.delim = "-")
data48 = CreateSeuratObject(counts = data48, project = "PSM2c", min.cells = 3, min.features = 200, names.delim = "-")
data49 = CreateSeuratObject(counts = data49, project = "PSM3a", min.cells = 3, min.features = 200, names.delim = "-")
data50 = CreateSeuratObject(counts = data50, project = "PSM3b", min.cells = 3, min.features = 200, names.delim = "-")
data51 = CreateSeuratObject(counts = data51, project = "PSM3c", min.cells = 3, min.features = 200, names.delim = "-")
data52 = CreateSeuratObject(counts = data52, project = "PSM4a", min.cells = 3, min.features = 200, names.delim = "-")
data53 = CreateSeuratObject(counts = data53, project = "PSM4b", min.cells = 3, min.features = 200, names.delim = "-")
data54 = CreateSeuratObject(counts = data54, project = "PSM4c", min.cells = 3, min.features = 200, names.delim = "-")
data55 = CreateSeuratObject(counts = data55, project = "PSM5a", min.cells = 3, min.features = 200, names.delim = "-")
data56 = CreateSeuratObject(counts = data56, project = "PSM5b", min.cells = 3, min.features = 200, names.delim = "-")
data57 = CreateSeuratObject(counts = data57, project = "PSM6a", min.cells = 3, min.features = 200, names.delim = "-")
data58 = CreateSeuratObject(counts = data58, project = "PSM6b", min.cells = 3, min.features = 200, names.delim = "-")
data59 = CreateSeuratObject(counts = data59, project = "PSM6c", min.cells = 3, min.features = 200, names.delim = "-")

# Merge all the seurat objects
GSE188172_object = merge(x = data1, y = c(data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14,data15,data16,data17,data18,data19,data20,
                                          data21,data22,data23,data24,data25,data26,data27,data28,data29,data30,data31,data32,data33,data34,data35,data36,data37,data38,data39,data40,
                                          data41,data42,data43,data44,data45,data46,data47,data48,data49,data50,data51,data52,data53,data54,data55,data56,data57,data58,data59),
                         add.cell.ids = c("CF1","CF2","CF3","CF4","CF5","CF6","CF7","CM1","CM2","CM3","CM4","CM5","PMF1","PMF2","PMF3","PMF4","PMF5","PMF6","PMF7","PMF8",
                                          "PMM1","PMM2","PMM3","PMM4","PMM5","PMM6","PMM7","PMM8","PSF1a","PSF1b","PSF1c","PSF2a","PSF2b","PSF2c","PSF3a","PSF3b","PSF4a","PSF4b",
                                          "PSF5a","PSF5b","PSF6a","PSF6b","PSM1a","PSM1b","PSM1c","PSM2a","PSM2b","PSM2c","PSM3a","PSM3b","PSM3c","PSM4a","PSM4b","PSM4c",
                                          "PSM5a","PSM5b","PSM6a","PSM6b","PSM6c"), project = "merge_data")

# Identify the mitochondria-related genes in GSE188172_object
GSE188172_object[["percent.mt"]] <- PercentageFeatureSet(GSE188172_object, pattern = "^MT-")

# Visualization of nFeature_RNA, nCount_RNA and percent.mt for QC
VlnPlot(object = GSE188172_object, features = c("nCount_RNA","nFeature_RNA","percent.mt"))
GSE188172_object <- subset(GSE188172_object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)

# Normalization following the LogNormalize method
GSE188172_norm <- NormalizeData(GSE188172_object, normalization.method = "LogNormalize", scale.factor = 10000)

# Scaling the expression with all genes in GSE188172_norm
GSE188172_norm <- FindVariableFeatures(GSE188172_norm, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GSE188172_norm)
GSE188172_norm <- ScaleData(GSE188172_norm, features = all.genes)

# PCA for dimensional reduction with the highly variable features in GSE188172_norm
GSE188172_norm <- RunPCA(GSE188172_norm, features = VariableFeatures(object = GSE188172_norm))

# Find neighbors with top20 PCs and conduct clustering with resolution 0.5
GSE188172_norm <- FindNeighbors(GSE188172_norm, dims = 1:20)
GSE188172_norm <- FindClusters(GSE188172_norm, resolution = 0.5)

# Conduct UMAP and visualization
GSE188172_norm <- RunUMAP(GSE188172_norm, dims = 1:20, seed.use = 1)
DimPlot(GSE188172_norm, reduction = "umap", shuffle = TRUE, raster = FALSE)

# Find markers in each clusters, and save the results for manually cell annotation
GSE188172_markers = FindAllMarkers(object = GSE188172_norm, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(x = GSE188172_markers, file = "path/to/YourDirectory/marker.csv")

# Manually cell annotation by searching cell markers
GSE188172_norm$label = ""
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "0"] = "NK" # Nature killer cell
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "1"] = "CD4T" # CD4+ T cell
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "2"] = "CD8T" # CD8+ T cell
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "3"] = "Bcell"
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "4"] = "CD14-Mo" # CD14+CD16- monocyte
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "5"] = "CD14CD16-Mo" # CD14+CD16+ monocyte
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "6"] = "CD4T"# CD4+ T cell
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "7"] = "CD14-Mo" # CD14+CD16- monocyte
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "8"] = "CD14-Mo" # CD14+CD16- monocyte
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "9"] = "CD4T" # CD4+ T cell
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "10"] = "CD16-Mo" # CD14-CD16+ monocyte
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "11"] = "Platelet"
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "12"] = "Pro-CD8T" # Proliferative CD8+ T cell
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "13"] = "Neutrophil"
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "14"] = "RB" # Red blood cell
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "15"] = "Bcell"
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "16"] = "DC" # Dendritic cell
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "17"] = "CD14-Mo" # CD14+CD16- monocyte
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "18"] = "PC" # Plasma cell
GSE188172_norm$label[GSE188172_norm$seurat_clusters == "19"] = "PC" # Plasma cell

# Label with control and disease group
GSE188172_norm$group = ""
GSE188172_norm$group[GSE188172_norm$orig.ident == "CF1"] = "Control"
GSE188172_norm$group[GSE188172_norm$orig.ident == "CF2"] = "Control"
GSE188172_norm$group[GSE188172_norm$orig.ident == "CF3"] = "Control"
GSE188172_norm$group[GSE188172_norm$orig.ident == "CF4"] = "Control"
GSE188172_norm$group[GSE188172_norm$orig.ident == "CF5"] = "Control"
GSE188172_norm$group[GSE188172_norm$orig.ident == "CF6"] = "Control"
GSE188172_norm$group[GSE188172_norm$orig.ident == "CF7"] = "Control"
GSE188172_norm$group[GSE188172_norm$orig.ident == "CM1"] = "Control"
GSE188172_norm$group[GSE188172_norm$orig.ident == "CM2"] = "Control"
GSE188172_norm$group[GSE188172_norm$orig.ident == "CM3"] = "Control"
GSE188172_norm$group[GSE188172_norm$orig.ident == "CM4"] = "Control"
GSE188172_norm$group[GSE188172_norm$orig.ident == "CM5"] = "Control"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PMF1"] = "Mild/Moderate"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PMF2"] = "Mild/Moderate"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PMF3"] = "Mild/Moderate"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PMF4"] = "Mild/Moderate"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PMF5"] = "Mild/Moderate"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PMF6"] = "Mild/Moderate"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PMF7"] = "Mild/Moderate"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PMF8"] = "Mild/Moderate"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PMM1"] = "Mild/Moderate"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PMM2"] = "Mild/Moderate"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PMM3"] = "Mild/Moderate"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PMM4"] = "Mild/Moderate"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PMM5"] = "Mild/Moderate"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PMM6"] = "Mild/Moderate"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PMM7"] = "Mild/Moderate"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PMM8"] = "Mild/Moderate"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSF1a"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSF1b"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSF1c"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSF2a"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSF2b"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSF2c"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSF3a"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSF3b"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSF4a"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSF4b"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSF5a"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSF5b"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSF6a"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSF6b"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSM1a"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSM1b"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSM1c"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSM2a"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSM2b"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSM2c"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSM3a"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSM3b"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSM3c"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSM4a"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSM4b"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSM4c"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSM5a"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSM5b"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSM6a"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSM6b"] = "Critical"
GSE188172_norm$group[GSE188172_norm$orig.ident == "PSM6c"] = "Critical"

# This is the end of the Data Processing in GSE188172
#################################################################################################################