# R version 4.2.2
# Load packages
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggpubr)

#################################################################################################################
# GSE149878 analysis

# Load the dataset
# GSE149878_norm.rds can be acquired by following the code in Section1_scRNA_snRNA_DataProcessing.R
GSE149878_norm = readRDS(file = "path/to/GSE149878/GSE149878_norm.rds")

# Select plasma cells
PC = subset(GSE149878_norm, subset = label == "PC")

# Check the expression of GAS6 in PC data
# We can find the mean value of GAS6 in PC data is 0.2301
summary(FetchData(object = PC, vars = "GAS6"))

# Define cells with lower GAS6 expression than 0.2301 in PC data are GAS6.neg
PC$GAS6.group <- "GAS6.pos"
PC$GAS6.group[WhichCells(PC, expression = GAS6 < 0.2301)] <- "GAS6.neg"
PC = subset(PC, subset = GAS6.group == "GAS6.pos")

# Add the GAS6.group id back to GSE149878_norm data
allcellname = colnames(GSE149878_norm)
for (i in 1:length(allcellname)){
  if (allcellname[i] %in% colnames(PC)){
    GSE149878_norm$label[i] = "PC_GAS6.pos"
  }
}

# Find markers in PC_GAS6.pos population
GSE149878_norm@active.ident = as.factor(x = GSE149878_norm$label)
markers = FindAllMarkers(object = GSE149878_norm, logfc.threshold = 0, only.pos = FALSE)

# Save the markers results for further analysis in Section10_GeneSetEnrichmentAnalysis.R
write.csv(x = markers, file = "C:/Users/user/Desktop/marker.csv")

# In the final part of this section, we overlap the top200 markers (ranked by log2 fold change) in PC_GAS6.pos population in each dataset (GSE171524, GSE155249, GSE149878)
# Then, get the 54 common markers in COVID-19 GAS6+ plasma cells, and regarded as  COVID-19 GAS6+ plasma cells signature

# This is the end of the GAS6 Plasma Cell in GSE149878
#################################################################################################################