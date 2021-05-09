# Project Description

This repository contains the scripts used to process and analyze the data.

O’Meara, C. C., Wamstad, J. A., Gladstone, R. A., Fomovsky, G. M., Butty, V. L., Shrikumar, A., Gannon, J. B., Boyer, L. A., & Lee, R. T. (2015). Transcriptional Reversion of Cardiac Myocyte Fate During Mammalian Cardiac Regeneration. Circulation Research, 116(5), 804–815. https://doi.org/10.1161/CIRCRESAHA.116.304269


# Contributor

Sheila Yee scyee@bu.edu

# Repository Contents

#### run_tophat.qsub ####
Implements TopHat, a splice aligner, to align RNA-seq reads from P0_1 sample to the mouse reference genome (mm9). Other utilities - Bowtie2, Boost, Python, and SAMtools - are also implemented. 

To run, submit as a job on the cluster: 
```
qsub run_tophat.qsub
```
This outputs *accepted_hits.bam*.

#### Determine alignment statistics ####

In the directory where the TopHat results are located, run on the command line: 
```
samtools flagstat P0_1_tophat/accepted_hits.bam 
```

#### Index BAM file ####
In the directory where the TopHat results are located, run on the command line: 
```
samtools index accepted_hits.bam
```
This outputs *accepted_hits.bam.bai*. 

#### rseqc.qsub ####
Performs quality control analysis on RNA-seq data and alignment.

To run, submit as a job on the cluster: 
```
qsub fastqc.sh
```
This outputs *gene_body.geneBodyCoverage.curves.pdf* and *inner_distance.inner_distance_plot.pdf* plots. 

#### run_cufflinks.qsub ####
Determine FPKM values of reads mapped to genomic regions. Run Cufflinks on *P0_1_tophat/accepted_hits.bam*.

To run, submit as a job on the cluster: 
```
qsub run_cufflinks.qsub
```
This outputs *genes.fpkm_tracking* file which contains FPKM values for all genes. 

#### run_cuffdiff.qsub ####
Identify differentially expressed genes in the P0_1 sample. 

To run, submit as a job on the cluster: 
```
qsub run_cuffdiff.qsub
```
This outputs  "*genes_exp.diff*. 

#### FPKM_histogram.R ####
Creates a histogram of FPKM values from the P0_1 sample. Requires *genes.fpkm_tracking* to run. Run using RStudio. 

#### FPKM_sar_mit_cc_plots.R ####
Creates plot to visualize FPKM values of representative sarcomere, mitochondrial, and cell cycle genes of interest during in vivo maturation across all neonatal and adult stages. Recreates Figure 1D from O'Meara et al. Requires *P0_1_cufflinks/genes.fpkm_tracking* and *genes.fpkm_tracking* for the remaining 7 samples (provided by instructor) to run. Run using RStudio. 

#### FPKM_heatmap.R ####
Creates a heatmap of the top 1000 and top 10 differentially expressed genes across all neonatal and adult stages (8 samples total). Recreates Figure 2A from O'Meara et al. Requires *fpkm_matrix.csv* (provided by instructor), *gene_exp.diff*, and *genes.fpkm_tracking* to run. Run using RStudio. 
