# Cell proportion calculation and statistics
# Here, we use GSE171524 dataset as an example

# R version 4.2.2
# Load packages
library(Seurat)

# Load the dataset
# GSE171524_norm.rds can be acquired by following the code in Section1_scRNA_snRNA_DataProcessing.R
GSE171524_norm = readRDS(file = "path/to/GSE171524/GSE171524_norm.rds")

# Use table() function to gain the cell number in each sample
# For instance, here's the sample 1_1
sample1 = subset(x = GSE171524_norm, subset = id == "1_1")
table(sample1$label)

# Record the cell number in each sample and calculate the proportion of each cell type in each sample
# The cell proportion (per sample) is calculated  by (cell number of specific cell type) / (total cell number in the sample) 
# If there are 27 samples in the dataset, this process require to be conducted for 27 times in different samples (sample1 to sample27)

# R version 4.3.1
# Load packages
library(ggplot2)
library(rstatix)
library(ggpubr)
library(scales)

# Based on the cell ratio that we have gain from the above process, we then evaluated the cell proportion with the following codes
# For instance, this is the example of AEP1 cell proportion from all samples in GSE171524
Control <- c(0.132731747,	0.202747618,	0.038623596, 0.146347032,	0.091468777,	0.003422983,	0.251618292)
COVID19 <- c(0.156638326,	0.00318218,	0.010192837,	0.005466871,	0.013761468,	0.023344764,	0.021096215,
             0.013503738,	0.00332871,	0.158709677,	0.003034901,	0.000775394,	0.017071164,	0.124784217, 
             0.008232711,	0.022580645,	0.001065341,	0.042247908,	0.073868613,	0.00458059)

# Create a data frame
Cell_proportion <- data.frame( 
  group = c(rep("Control",7), rep("COVID-19",20)),
  Proportion = c(Control, COVID19)
)

# Conduct statistics with Wilcoxon rank-sum test for evaluating whether there is significant difference between control and disease group
stat.test <- wilcox_test(data = Cell_proportion, Proportion ~ group, exact = TRUE)
stat.test <- add_xy_position(stat.test, x = "group")

# Draw box plot and input the 
# Adjust the max_value according to the y axis in each cell type proportion result
# For example, we used 0.3 in AEP1 cell proportion for visualization
max_value = 0.3
plot = ggboxplot(stat.test, x = "group", y = "Proportion", fill = "group", ylim = c(-0.01,max_value)) +
  xlab("") +
  ylab("Cell proportion") +
  geom_jitter(shape = 16, position = position_jitter(0.1)) +
  scale_fill_manual(values = c("#b3cde3","#fbb4ae"))+
  theme(axis.text = element_text(size = 30), 
        axis.title = element_text(size = 30),
        plot.title = element_text(hjust = 0.5, size = 30)) +
  stat_pvalue_manual(stat.test, label = "p = {p}", vjust = -0.2, y.position = c(max_value-0.02), size = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  ggtitle("Alveolar type 1 epithelial cell")

# Save the box plot
ggsave(filename = "path/to/YourDirectory/plot.png", plot = plot, width = 7, height = 7, dpi = 800)

# This is the end of the Cell Proportion Calculation in GSE171524
#################################################################################################################
