# Loading data into R -----------------------------------------------------
setwd("/projectnb2/bf528/users/hedgehog/project_2/programmer/testing/cufflinks/P0_1_cufflinks")
genes.fpkm <- read.table("genes.fpkm_tracking")


# Data cleaning -----------------------------------------------------------
colnames(genes.fpkm) <- genes.fpkm[1,] # column names are in the first row
genes.fpkm <- genes.fpkm[-1,] # get rid of this row
genes.fpkm$FPKM <- as.numeric(genes.fpkm$FPKM) # convert to numeric

# Histogram ---------------------------------------------------------------
fpkm <- genes.fpkm$FPKM
# Only inlcuding FPKM values >=1, as specified by paper
fpkm <- fpkm[which(fpkm>=1)] 
# reduces length to ~38% of original
hist(log10(fpkm),
     main = "Histogram of log-adjusted FPKM values")
