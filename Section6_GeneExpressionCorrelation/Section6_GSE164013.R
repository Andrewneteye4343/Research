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
names = c("")
gene_members = GSE164013[names]

# list the gene name for these gene members
colnames(gene_members) = c("")

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
