# author: Sheila Yee 
# project: Transcriptional Profile of Mammalian Cardiac Regeneration with mRNA-Seq¶

# Read FPKM data for sample P0_1
FPKM_P0_1 <- read.csv("/usr4/bf528/scyee/project_5/project_5/P0_1_cufflinks/genes.fpkm_tracking", header=TRUE, sep="\t") 

# Convert data to numeric value ----
FPKM_P0_1$FPKM <- as.numeric(FPKM_P0_1$FPKM) 

# Histogram ----
fpkm <- FPKM_P0_1$FPKM # there are 37469 FPKM values

# Select only FPKM values that is greater than 0
fpkm <- fpkm[which(fpkm>=1)] # after this filtering step, there are 14,205 FPKM values greater than 0

# Create histogram
# Use log10 scale to better visualize the large range of FPKM values
hist(log10(fpkm),
     main = "Log10-adjusted FPKM", 
     cex.main = 2,
     xlab = "FPKM (log10)",
     ylim = c(0, 1000),
     breaks=50, 
     col="powderblue")



