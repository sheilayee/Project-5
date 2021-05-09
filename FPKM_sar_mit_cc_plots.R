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

##### Plots of FPKM fold changes from Figure 1D ----
# Goal: determine how sarcomere, mitochondrial, and cell cycle genes significantly differentially expressed during in vivo maturation

##### Load files
# Read FPKM tables for P0, P4, P7, and Ad (8 samples total)
# P0_1 is processed by student, the remaining 7 samples were provided by instructor
P0_1 <- read.csv("/usr4/bf528/scyee/project_5/project_5/P0_1_cufflinks/genes.fpkm_tracking", header=TRUE, sep="\t") 
P0_2 <- read.csv("/project/bf528/project_2/data/samples/P0_2/genes.fpkm_tracking", header=TRUE, sep="\t")

P4_1 <- read.csv("/project/bf528/project_2/data/samples/P4_1/genes.fpkm_tracking", header=TRUE, sep="\t")
P4_2 <- read.csv("/project/bf528/project_2/data/samples/P4_2/genes.fpkm_tracking", header=TRUE, sep="\t")

P7_1 <- read.csv("/project/bf528/project_2/data/samples//P7_1/genes.fpkm_tracking", header=TRUE, sep="\t")
P7_2 <- read.csv("/project/bf528/project_2/data/samples/P7_2/genes.fpkm_tracking", header=TRUE, sep="\t")

Ad_1 <-read.csv("/project/bf528/project_2/data/samples/Ad_1/genes.fpkm_tracking", header=TRUE, sep="\t")
Ad_2 <- read.csv("/project/bf528/project_2/data/samples/Ad_2/genes.fpkm_tracking", header=TRUE, sep="\t")

# Read FPKM matrix
fpkm_matrix <- read.csv("/project/bf528/project_2/data/fpkm_matrix.csv", header=TRUE, sep="\t")

# Read file containing differentially expressed genes (an output of cuffdiff )
# This file is titled gene_diff_expr

##### Obtain gene names used by O'Meara et al. in Figure 1D ----
sar_genes <- c("Pdlim5", "Pygm", "Myoz2", "Des", "Csrp3", "Tcap", "Cryab")
mit_genes <- c("Mpc1", "Prdx3", "Acat1", "Echs1", "Slc25a11", "Phyh")
cell_cycle_genes <- c("Cdc7", "E2f8", "Cdk7", "Cdc26", "Cdc6", "Cdc27", "Bora", "Cdc45", "Rad51", "Aurkb", "Cdc23")

# Filtering function ----
# Create function to filter out genes of interest in each sample. Further filter each of these gene's respective mean FPKM value 
# Then, find the mean FPKM for this gene across both samples for each age
# Repeat for all ages and combine into one super dataset

filter_genes_FPKM <- function(sample_1, sample_2, gene_names, age){
  sample_1_filtered <- sample_1 %>% filter(gene_short_name %in% gene_names) %>% # filter out genes of interest
    select(gene_short_name, FPKM) %>% # select for genes of interest and their associated FPKM value
    column_to_rownames("gene_short_name") # convert column of gene names into name for each row of FPKM value
  sample_2_filtered <- sample_2 %>% filter(gene_short_name %in% gene_names) %>%  # repeat the same steps for the second sample for each age
    select(gene_short_name, FPKM) %>% 
    column_to_rownames("gene_short_name")
  combined_samples <- bind_cols(sample_1_filtered,sample_2_filtered) %>% # create new dataframes by combining the filtered FPKM values from both samples
    rownames_to_column("Genes") %>% # convert rows of gene names into column titled "Genes" 
    rowwise() %>% 
    mutate(FPKM = mean(c(FPKM...1, FPKM...2))) %>% # find the mean FPKM value for each gene of interest across both samples
    select(Genes, FPKM) %>%
    add_column(Age=age) # add column to specify the age associated with each FPKM value (either P0, P4, P7, or Ad)
  return(combined_samples)
}

##### Sarcomere ----
sar_P0_12 <- filter_genes_FPKM(sample_1=P0_1, sample_2=P0_2, gene_names=sar_genes, age='P0')
sar_P4_12 <- filter_genes_FPKM(sample_1=P4_1, sample_2=P4_2, gene_names=sar_genes, age='P4')
sar_P7_12 <- filter_genes_FPKM(sample_1=P7_1, sample_2=P7_2, gene_names=sar_genes, age='P7')
sar_AD_12 <- filter_genes_FPKM(sample_1=Ad_1, sample_2=Ad_2, gene_names=sar_genes, age='Ad')

# Combine all of the sarcomere FPKM values across all ages
sarcomere_genes_final_FPKM <- bind_rows(sar_P0_12, sar_P4_12, sar_P7_12, sar_AD_12)

##### Mitochondria ----
mit_P0_12 <- filter_genes_FPKM(sample_1=P0_1, sample_2=P0_2, gene_names=mit_genes, age='P0')
mit_P4_12 <- filter_genes_FPKM(sample_1=P4_1, sample_2=P4_2, gene_names=mit_genes, age='P4')
mit_P7_12 <- filter_genes_FPKM(sample_1=P7_1, sample_2=P7_2, gene_names=mit_genes, age='P7')
mit_AD_12 <- filter_genes_FPKM(sample_1=Ad_1, sample_2=Ad_2, gene_names=mit_genes, age='Ad')

# Combine all of the mitochondrial FPKM values across all ages
mitochondrial_genes_final_FPKM <- bind_rows(mit_P0_12, mit_P4_12, mit_P7_12, mit_AD_12)

##### Cell Cycle ----
cc_P0_12 <- filter_genes_FPKM(sample_1=P0_1, sample_2=P0_2, gene_names=cell_cycle_genes, age='P0')
cc_P4_12 <- filter_genes_FPKM(sample_1=P4_1, sample_2=P4_2, gene_names=cell_cycle_genes, age='P4')
cc_P7_12 <- filter_genes_FPKM(sample_1=P7_1, sample_2=P7_2, gene_names=cell_cycle_genes, age='P7')
cc_AD_12 <- filter_genes_FPKM(sample_1=Ad_1, sample_2=Ad_2, gene_names=cell_cycle_genes, age='Ad')

# Combine all of the cell cycle FPKM values across all ages
cell_cycle_genes_final_FPKM <- bind_rows(cc_P0_12, cc_P4_12, cc_P7_12, cc_AD_12)

###### Plots -----

# Sarcomere plot -----
sarcomere_genes_final_FPKM$Age = as.factor(sarcomere_genes_final_FPKM$Age) # categorize data and store it as levels
sarcomere_genes_final_FPKM$Age <- factor(sarcomere_genes_final_FPKM$Age, levels=c('P0','P4','P7','Ad'))
sarcomere_plot <- ggplot(sarcomere_genes_final_FPKM, aes(x=Age, y=FPKM, col=Genes, group=Genes)) + 
  geom_point() + 
  geom_path() + 
  labs(x='', y="FPKM") + 
  ylim(0,1500) + # set same y-axis scale as figure from original paper 
  guides(x = guide_axis(angle = 45)) + # tilt x-axis tick marks
  ggtitle('Sarcomere') + 
  theme(plot.title = element_text(hjust=0.5, size=18), 
        axis.title=element_text(size=14),  
        axis.text.x = element_text(color="black", size=12, angle=45), # angle x-axis tick mark labels
        axis.text.y = element_text(color="black", size=12), 
        axis.ticks.length =unit(2,"mm"), # define the space between tick mark labels and tick marks
        axis.line = element_line(colour = "black"), # set the x and y axis to be black
        panel.background = element_blank(), # remove gray background on plot
        panel.grid.major = element_blank(), # remove grid lines
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), # remove gray border
        legend.title = element_blank(), # remove legend title
        legend.key = element_rect(fill = NA, color = NA), # remove gray background around legend
        aspect.ratio = 0.8) + # set aspect ratio
  labs(tag = "A") + # add figure label (outside of the plot)
  scale_color_manual(values=c("black", "turquoise3", "hotpink2", "red2", "blue3", "chartreuse3", "yellow")) # assign color to lines in the order of the legend, use similar color scheme as authors
(sarcomere_plot)

# Mitochondria plot ----
mitochondrial_genes_final_FPKM$Age = as.factor(mitochondrial_genes_final_FPKM$Age)
mitochondrial_genes_final_FPKM$Age <- factor(mitochondrial_genes_final_FPKM$Age, levels=c("P0", "P4", "P7","Ad"))
mitochondrial_plot <- ggplot(mitochondrial_genes_final_FPKM, aes(x=Age, y=FPKM, col=Genes, group=Genes)) + 
  geom_point() + 
  geom_path() +
  labs(x='', y="FPKM") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Mitochondria") +
  theme(plot.title = element_text(hjust=0.5, size=18), 
        axis.title=element_text(size=14),  
        axis.text.x = element_text(color="black", size=12, angle=45),
        axis.text.y = element_text(color="black", size=12), 
        axis.ticks.length =unit(2,"mm"), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA, color = NA), 
        aspect.ratio = 0.8) +
  labs(tag = "B") +
  scale_color_manual(values=c("red2", "hotpink2", "black", "chartreuse3", "turquoise3")) 
(mitochondrial_plot)

# Cell cycle plot ----
cell_cycle_genes_final_FPKM$Age = as.factor(cell_cycle_genes_final_FPKM$Age)
cell_cycle_genes_final_FPKM$Age <- factor(cell_cycle_genes_final_FPKM$Age, levels=c("P0", "P4", "P7","Ad"))
cell_cycle_plot <- ggplot(cell_cycle_genes_final_FPKM, aes(x=Age, y=FPKM, col=Genes, group=Genes)) + 
  geom_point() + 
  geom_path() +
  scale_y_continuous(breaks=seq(0,50,10)) +
  # ylim(0,80) + 
  labs(x='', y="FPKM") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Cell Cycle") +
  theme(plot.title = element_text(hjust=0.5, size=18), 
        axis.title=element_text(size=14),  
        axis.text.x = element_text(color="black", size=12, angle=45),
        axis.text.y = element_text(color="black", size=12), 
        axis.ticks.length =unit(2,"mm"), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA, color = NA), 
        aspect.ratio = 0.8) +
  labs(tag = "C") +
  scale_color_manual(values=c("black", "bisque3", "deeppink4", "hotpink2", "pink1", "turquoise3", "blue3", "red2", "yellow1", "gray71")) 
(cell_cycle_plot)

# Combine all 3 plots
plot_grid(sarcomere_plot, mitochondrial_plot, cell_cycle_plot, ncol = 1) # width=2024, height=1160, scale=1

