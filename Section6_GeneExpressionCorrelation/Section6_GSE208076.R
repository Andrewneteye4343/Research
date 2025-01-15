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
# GSE208076 data

# Because GSE208076 is raw counts data, we need to normalize it first
# Load GSE208076 raw data
GSE208076_raw = read.csv(file = "path/to/GSE208076/GSE208076_raw.csv")

# Define group types
sampletypes = factor(c(rep("Control",3), rep("COVID-19",7)))

# check the sample and gene numbers
n_samples = ncol(GSE208076_raw)
n_genes = nrow(GSE208076_raw)

# Create DGEList
dge = DGEList(GSE208076_raw)

# Filter low expressed genes
keep.exprs = filterByExpr(dge, group = sampletypes)
dge = dge[keep.exprs,]

# TMM normalization
# norm_dge can be saved as GSE208076_norm.csv
norm_dge = calcNormFactors(dge, method = "TMM")

# Gene Expression Correlation
# Load gene expression matrix
GSE208076_norm = read.csv(file = "path/to/GSE208076/GSE208076_norm.csv", header = TRUE, row.names = 1)

# Transpose the matrix
GSE208076_norm = as.data.frame(t(GSE208076_norm))

# list the interested gene members
names = c("ENSG00000183087","ENSG00000128383","ENSG00000143546",
          "ENSG00000163220","ENSG00000135218","ENSG00000072274","ENSG00000115414","ENSG00000010278",
          "ENSG00000135404","ENSG00000177575","ENSG00000100600","ENSG00000136235","ENSG00000118785",
          "ENSG00000170323","ENSG00000122641","ENSG00000108821","ENSG00000164692","ENSG00000107796")
gene_members = GSE208076_norm[,names]

# list the gene name for these gene members
colnames(gene_members) = c("GAS6","APOBEC3A","S100A8","S100A9","CD36","TFRC","FN1","CD9","CD63","CD163","LGMN","GPNMB",
                           "SPP1","FABP4","INHBA","COL1A1","COL1A2","ACTA2")

# Calculate the Pearson correlation
correlation_results = cor(x =  gene_members, method = "pearson")

# rounded to one decimal point
corr <- round(test_cor, 1)
p = ggcorrplot(corr, hc.order = FALSE, 
               outline.color = "black", 
               colors = c("#6D9EC1", "white", "#E46726"),
               lab = TRUE)

# This is the end of the Data Processing in GSE208076
#################################################################################################################
