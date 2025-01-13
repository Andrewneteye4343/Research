# R version 4.2.2
# Load packages
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggpubr)

#################################################################################################################
# GSE171524 analysis

# Load the dataset
# GSE171524_norm.rds can be acquired by following the code in Section1_scRNA_snRNA_DataProcessing.R
GSE171524_norm = readRDS(file = "path/to/GSE171524/GSE171524_norm.rds")

# Select cells (ex: CD14CD16-Mo, PF-M,...)
interested_cells = subset(x = GSE171524_norm, subset = label == "CD14CD16-Mo" | label == "PF-M" | label == "NA-M" | label == "HYP-M" | label == "VCAN-M" | label == "TR-M")

# Define gene member that ypu want to evaluate
# For instance, these genes represent the profibrotic macrophage population
cd_features <- list(c("CD9","TREM2","SPP1","LGMN","GPNMB","FABP5","CD63","CD163"))

# Use AddModuleScore() function to calcualte the module score
# Use random 500 features as control gene expression
# Set seed as 123
interested_cells <- AddModuleScore(
  object = interested_cells,
  features = cd_features,
  ctrl = 500,
  name = 'CD_Features',
  seed = 123
)

# Define a function that can draw Violin plot and add statistics on it
vp_case1 <- function(gene_signature, file_name, test_sign, y_max){
  plot_case1 <- function(signature){
    VlnPlot(interested_cells, features = signature,
            pt.size = 0.1, 
            group.by = "label",  # Choose the comparison group. Here, we compare it between cell populations in interested_cells
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + 
      stat_compare_means(comparisons = test_sign, method = "t.test") + 
      xlab(label = "") +
      labs(title = "Profibrotic macrophage signature") +
      stat_summary(fun=mean, geom="point", size=1, color = "red")
  }
  purrr::map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = .)
  file_name <- paste0(file_name, "_r.png")
  ggsave(file_name, width = 7, height = 7, dpi = 800, bg = "white")
}

# Define the gene that your are going to visualize
# CD_Features1 is the default expression id for AddModuleScore()
gene_sig <- c("CD_Features1")

# Choose which comparison will be add the statistics label
comparisons <- list(c("CD14CD16-Mo","PF-M"), c("NA-M","PF-M"), c("HYP-M","PF-M"), c("VCAN-M","PF-M"), c("TR-M","PF-M"))

# Conduct the visualization function
vp_case1(gene_signature = gene_sig, file_name = "path/to/YourDirectory/gene_sig", 
         test_sign = comparisons, y_max = 5.1)

# This is the end of the Module Score Calculation in GSE171524
#################################################################################################################
