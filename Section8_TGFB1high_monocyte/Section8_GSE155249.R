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

# Select CD14+CD16+ monocyte
monocyte = subset(GSE155249_norm, subset = label == "CD14CD16-Mo2")

# Check the expression of TGFB1 in monocyte data
# We can find the mean value of TGFB1 in monocyte data is 0.733
summary(FetchData(object = monocyte, vars = "TGFB1"))

# Define cells with lower TGFB1 expression than 0.733 in monocyte data are TGFB1.neg
monocyte$TGFB1.group <- "TGFB1.pos"
monocyte$TGFB1.group[WhichCells(monocyte, expression = TGFB1 < 0.733)] <- "TGFB1.neg"

# To select TGFB1 highly expressing monocytes (considered as differentiated monocytes), we further select the population with TGFB1 expression
# Select monocyte expresses higher TGFB1
monocyte_2 = subset(monocyte, TGFB1.group == "TGFB1.pos")

# Check the expression of TGFB1 in monocyte_2 data
# We can find the median value of TGFB1 in monocyte_2 data is 1.4640
summary(FetchData(object = monocyte_2, vars = "TGFB1"))

# Define cells with lower TGFB1 expression than 1.4640 in monocyte data are TGFB1.neg
monocyte_2$TGFB1.group <- "TGFB1.pos"
monocyte_2$TGFB1.group[WhichCells(monocyte_2, expression = TGFB1 < 1.4640)] <- "TGFB1.neg"

# Add the "TGFB1.group" label back to GSE155249_norm
group1 = subset(monocyte_2, subset = TGFB1.group == "TGFB1.pos")
group2 = subset(monocyte_2, subset = TGFB1.group == "TGFB1.neg")
group3 = subset(monocyte, subset = TGFB1.group == "TGFB1.neg")

allcellname = colnames(GSE155249_norm)
GSE155249_norm$TGFB_label = ""

for (i in 1:length(allcellname)){
  if (allcellname[i] %in% colnames(group1)){
    GSE155249_norm$TGFB_label[i] = "TGFB1high_CD14CD16-Mo"
  }
  else if (allcellname[i] %in% colnames(group2)){
    GSE155249_norm$TGFB_label[i] = "TGFB1low_CD14CD16-Mo"
  }
  else if (allcellname[i] %in% colnames(group3)){
    GSE155249_norm$TGFB_label[i] = "TGFB1low_CD14CD16-Mo"
  }
  else {
    GSE155249_norm$TGFB_label[i] = GSE155249_norm$label[i]
  }
}

# Select the TGFB1high_CD14CD16-Mo and TGFB1low_CD14CD16-Mo cells
monocyte_3 = subset(x = GSE155249_norm, subset = TGFB_label == "TGFB1high_CD14CD16-Mo" | TGFB_label == "TGFB1low_CD14CD16-Mo")

# Reorder the cells in monocyte_3 data
monocyte_3$TGFB_label <- factor(
  monocyte_3$TGFB_label,
  levels = c("TGFB1low_CD14CD16-Mo", "TGFB1high_CD14CD16-Mo")
)

# Define a function that can draw Violin plot and add statistics on it
vp_case1 <- function(gene_signature, file_name, test_sign, y_max){
  plot_case1 <- function(signature){
    VlnPlot(monocyte_3, features = signature,
            pt.size = 0.1, 
            group.by = "TGFB_label", 
            y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
            cols = c("#b3cde3","#fbb4ae")
    ) + stat_compare_means(comparisons = test_sign, method = "t.test") + 
      stat_summary(fun.y=median, geom="point", size=1, color = "white") +
      
      xlab(label = "")
  }
  purrr::map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = .)
  file_name <- paste0(file_name, "_r.png")
  ggsave(file_name, width = 7, height = 7, dpi = 800, bg = "white")
}

# Define the gene that your are going to visualize
gene_sig <- c("TGFB1")

# Choose which comparison will be add the statistics label
comparisons <- list(c("TGFB1low_CD14CD16-Mo", "TGFB1high_CD14CD16-Mo"))

# Conduct the visualization function
vp_case1(gene_signature = gene_sig, file_name = "path/to/YourDirectory/gene_sig", 
         test_sign = comparisons, y_max = 4)

# Calculate the differential expression genes (DEGs) between TGFB1high CD14+CD16+ monocyte and TGFB1low CD14+CD16+ monocyte
TGFB1_monocyte_markers = FindMarkers(object = "monocyte_3", group = "TGFB_label", ident.1 = "TGFB1high_CD14CD16-Mo")

# Save DEGs for further analysis in Section10_GeneSetEnrichmentAnalysis.R
write.csv(x = TGFB1_monocyte_markers, file = "path/to/YourDirectory/markers.csv")

# This is the end of the TGFB1high monocyte in GSE155249
#################################################################################################################