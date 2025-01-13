# R version 4.2.2
# Load packages
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggpubr)

#################################################################################################################
# GSE155249 analysis

# Load the dataset
# GSE155249_norm.rds can be acquired by following the code in Section1_scRNA_snRNA_DataProcessing.R
GSE155249_norm = readRDS(file = "path/to/GSE155249/GSE155249_norm.rds")

# Select plasma cells
PC = subset(GSE155249_norm, subset = label == "PC")

# Select plasma cells in COVID-19
COVID19_PC = subset(PC, subset = group == "COVID-19")

# Check the expression of GAS6 in COVID19_PC data
# We can find the mean value of GAS6 in COVID19_PC data is 0.48
summary(FetchData(object = COVID19_PC, vars = "GAS6"))

# Define cells with lower GAS6 expression than 0.48 in COVID19_PC data are GAS6.neg
COVID19_PC$GAS6.group <- "GAS6.pos"
COVID19_PC$GAS6.group[WhichCells(COVID19_PC, expression = GAS6 < 0.48)] <- "GAS6.neg"
COVID19_PC_GAS6 = subset(COVID19_PC, subset = GAS6.group == "GAS6.pos")

# Add the GAS6.group id back to GSE155249_norm data
allcellname = colnames(GSE155249_norm)
for (i in 1:length(allcellname)){
  if (allcellname[i] %in% colnames(COVID19_PC_GAS6)){
    GSE155249_norm$label[i] = "PC_GAS6.pos"
  }
}

# Find markers in PC_GAS6.pos population
GSE155249_norm@active.ident = as.factor(x = GSE155249_norm$label)
markers = FindAllMarkers(object = GSE155249_norm, logfc.threshold = 0, only.pos = FALSE)

# Save the markers results for further analysis in Section10_GeneSetEnrichmentAnalysis.R
write.csv(x = markers, file = "C:/Users/user/Desktop/marker.csv")

# In the final part of this section, we overlap the top200 markers (ranked by log2 fold change) in PC_GAS6.pos population in each dataset (GSE171524, GSE155249, GSE149878)
# Then, get the 54 common markers in COVID-19 GAS6+ plasma cells, and regarded as  COVID-19 GAS6+ plasma cells signature

# This is the end of the GAS6 Plasma Cell in GSE155249
#################################################################################################################