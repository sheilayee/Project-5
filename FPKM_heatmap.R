# author: Sheila Yee
# project: Transcriptional Profile of Mammalian Cardiac Regeneration with mRNA-Seq¶

library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(pheatmap)
library(patchwork)
library(tibble)
require(grid)
require(gridExtra)
library(gtable)
library(cowplot)
  
##### Create a heatmap of the top 1000 differentially expressed genes ----
# Read FPKM matrix
FPKM_matrix <- read.csv('/project/bf528/project_2/data/fpkm_matrix.csv', sep =  "\t")

# Read file containing differentially expressed genes (an output of cuffdiff)
diff_expr_file <-read.table("/usr4/bf528/scyee/project_5/project_5/cuffdiff_out/gene_exp.diff", header=TRUE, sep= "\t") 

# Read FPKM data for P0_1
FPKM_0 <- read.csv("/usr4/bf528/scyee/project_5/project_5/P0_1_cufflinks/genes.fpkm_tracking", header=TRUE, sep="\t") 
FPKM_0 <- FPKM_0  %>%  select(FPKM, tracking_id, gene_short_name) # select the columns of interest

##### Top 1000 genes ----
# Obtain the top 1000 significant genes from P0 to adult
sig_de_genes  <- diff_expr_file[diff_expr_file$significant =='yes',]
x <- sig_de_genes %>% filter(sig_de_genes$log2.fold_change. > 0 ) # 1084 genes that are up-regulated
y <- sig_de_genes %>% filter(sig_de_genes$log2.fold_change. < 0 ) # 1055 genes that are down-regulated
sig_de_genes <- diff_expr_file %>% filter(diff_expr_file$significant =='yes') %>% arrange(q_value) %>% slice_head(n=1000) # change to 10 when looking at top 10 genes

# Rename the columns in both FPKM tables
colnames(FPKM_0) <- c("P0_1", "tracking_id", "gene_short_name") # rename tracking_id 
colnames(FPKM_matrix) <- c("tracking_id", "Ad_1", "Ad_2", "P0_2", "P4_1", "P4_2", "P7_1", "P7_2")

# Create a super FPKM matrix of all 8 samples by combining P0_1 with the other 7 samples based on the tracking_id
FPKM_8_samples <- merge(FPKM_0, FPKM_matrix, by = 'tracking_id')
FPKM_8_samples <- FPKM_8_samples %>% relocate(P0_1, .after = Ad_2) # move the P0_1 column to be after Ad_2 and before P0_2
FPKM_8_samples$tracking_id <- NULL # remove the tracking_id column (this is the first column)
FPKM_8_samples <- FPKM_8_samples %>%
  group_by(gene_short_name) %>%
  summarise_all(mean)


# Filter the super FPKM for these 1000 genes
FPKM_1000 <- FPKM_8_samples %>% filter(gene_short_name %in% sig_de_genes$gene)

# Label each row in the FPKM as the gene name
rownames(FPKM_1000) <- NULL
FPKM_1000 <- data.frame(column_to_rownames(FPKM_1000, var = "gene_short_name"))

# Convert the FPKM dataframe to a to matrix in preparation for the heatmap
FPKM_matrix_1000 <- as.matrix(FPKM_1000[,])

# Set the colors for the heatmap 
my_colors = brewer.pal(n = 11, name = "BrBG")
my_colors = colorRampPalette(my_colors)(100)

# Create the heatmap
my_heatmap <- pheatmap(FPKM_matrix_1000, 
                       scale = "row", 
                       fontsize_row = 5,
                       border_color = NA,
                       color = my_colors,
                       show_rownames = FALSE, # hide gene names for rows b/c it looks too busy. Unhide for top 10 genes heatmap
                       clustering_distance_rows="euclidean",
                       clustering_distance_cols="euclidean")


