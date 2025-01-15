# R version 4.1.2
# Load packages
library(data.table)
library(limma)
library(edgeR)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(ggcorrplot)

#################################################################################################################
# GSE164013 data

# Gene Expression Correlation
# Because GSE164013 is Q3 normalized data, we can directly use it
GSE164013 = read.csv(file = "path/to/GSE164013/GSE164013.csv", header = TRUE, row.names = 1)

# Transpose the matrix
GSE164013 = as.data.frame(t(GSE164013))

# list the interested gene members
names = c("GAS6","APOBEC3A","S100A8","S100A9","CD68","CD36","TFRC","FN1","CD9","CD63","CD163","LGMN","GPNMB",
          "SPP1","FABP4","INHBA","COL1A1","COL1A2","ACTA2")
gene_members = GSE164013[names]

# list the gene name for these gene members
colnames(gene_members) = c("GAS6","APOBEC3A","S100A8","S100A9","CD68","CD36","TFRC","FN1","CD9","CD63","CD163","LGMN","GPNMB",
                           "SPP1","FABP4","INHBA","COL1A1","COL1A2","ACTA2")

# Calculate the Pearson correlation
correlation_results = cor(x =  gene_members, method = "pearson")

# rounded to one decimal point
corr <- round(correlation_results, 1)
p = ggcorrplot(corr, hc.order = FALSE, 
               outline.color = "black", 
               colors = c("#6D9EC1", "white", "#E46726"),
               lab = TRUE)

# This is the end of the Data Processing in GSE164013
#################################################################################################################
