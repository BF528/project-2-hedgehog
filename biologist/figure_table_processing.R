# Taylor Falk
# tfalk@bu.edu
# biologist, BF528 prj2

library(gplots)
library(gridExtra)
library(reshape2)
library(tidyverse)


setwd("/projectnb/bf528/users/hedgehog/project_2/biologist/")

# criminal that someone would provide me with a csv that is not comma separated
fpkm <- read.csv("/project/bf528/project_2/data/fpkm_matrix.csv", 
                 header = TRUE, sep = "\t")

#### figure 1D recreate ####
genes <- c("Pdlim5", "Pygm", "Myoz2", "Des", "Csrp3", "Tcap", "Cryab", # sarc
           "Mpc1", "Prdx3", "Acat1", "Echs1", "Slc25a11", "Phyh", # mitochondria
           "Cdc7", "E2f8", "Cdk7", "Cdc26", "Cdc6", "E2f1", "Cdc27", # cellcycle
           "Bora", "Cdc45", "Rad51", "Aurkb", "Cdc23") # cellcycle

# find all available FPKM files to iterate through
files <- Sys.glob("*fpkm_tracking")
headerName <- sapply(files, function(x) strsplit(x, ".", fixed = TRUE)[[1]])[1,]

specificGenes <- data.frame(row.names = genes)

# load in FPKM data from all files
for (file in files) {
  table <- read.table(file, header = T)
  fpkm_vec <- c()
  for (gene in genes) {
    currentVal <- table$FPKM[which(table$gene_short_name == gene)]
    if (!(length(currentVal)) > 0) {
      currentVal <- NA
    }
    fpkm_vec <- c(fpkm_vec, currentVal)
  }
  specificGenes <- cbind(specificGenes, fpkm_vec)
}
names(specificGenes) <- headerName

#### plot details ####

plotFPKM <- function(genes, data, palette, samples, title) {
  geneData <- data[genes, samples]
  geneData <- melt(rownames_to_column(geneData))
  names(geneData) <- c("Gene", "mouse", "fpkm")
  
  plot <- ggplot(data=geneData, aes(x=mouse, y=fpkm, color = Gene)) +
    geom_point() + 
    geom_line(aes(group = Gene)) + 
    ggtitle(title) + 
    theme_bw() + 
    xlab("Mouse sample") +
    ylab("FPKM") +
    scale_color_brewer(palette = palette)
  return(plot)
}

sarc1 <- plotFPKM(genes[1:7], specificGenes, "Accent",
                  c("P0_1", "P4_1", "P7_1", "Ad_1"), "Sarcomere 1")
mito1 <- plotFPKM(genes[8:13], specificGenes, "Set2",
                  c("P0_1", "P4_1", "P7_1", "Ad_1"), "Mitochondria 1")
cc1   <- plotFPKM(genes[14:25], specificGenes, "Set3",
                  c("P0_1", "P4_1", "P7_1", "Ad_1"), "Cell Cycle 1")

sarc2 <- plotFPKM(genes[1:7], specificGenes, "Accent",
                  c("P0_2", "P4_2", "P7_2", "Ad_2"), "Sarcomere 2")
mito2 <- plotFPKM(genes[8:13], specificGenes, "Set2",
                  c("P0_2", "P4_2", "P7_2", "Ad_2"), "Mitochondria 2")
cc2 <- plotFPKM(genes[14:25], specificGenes, "Set3",
                c("P0_2", "P4_2", "P7_2", "Ad_2"), "Cell Cycle 2")

# combine all plots
grid.arrange(grobs = list(sarc1, sarc2, mito1, mito2, cc1, cc2))

#### subsetting + heatmaping ####

files <- Sys.glob("*fpkm_tracking")
headerName <- sapply(files, function(x) paste0(".", strsplit(x, ".", fixed = TRUE)[[1]]))[1,]

# R doesn't let you place objects into a list??? so we have to copy and paste 
# like neanderthals
t1 <- read.table(files[1], header = T)[,c("gene_id", "gene_short_name", "FPKM")]
t2 <- read.table(files[2], header = T)[,c("gene_id", "FPKM")]
t3 <- read.table(files[3], header = T)[,c("gene_id", "FPKM")]
t4 <- read.table(files[4], header = T)[,c("gene_id", "FPKM")]
t5 <- read.table(files[5], header = T)[,c("gene_id", "FPKM")]
t6 <- read.table(files[6], header = T)[,c("gene_id", "FPKM")]
t7 <- read.table(files[7], header = T)[,c("gene_id", "FPKM")]
t8 <- read.table(files[8], header = T)[,c("gene_id", "FPKM")]

# and merging these by gene_id gets messy but produces a nice table
merge1 <- merge(t1, t2, by="gene_id", all=T, suffixes = c(".Ad_1", ".Ad_2"))
merge2 <- merge(merge1, t3, by="gene_id", all=T)
merge3 <- merge(merge2, t4, by="gene_id", all=T, suffixes = c(".P0_1", ".P0_2"))
merge4 <- merge(merge3, t5, by="gene_id", all=T)
merge5 <- merge(merge4, t6, by="gene_id", all=T, suffixes = c(".P4_1", ".P4_2"))
merge6 <- merge(merge5, t7, by="gene_id", all=T)
merge_final <- merge(merge6, t8, by="gene_id", all=T, suffixes = c(".P7_1", ".P7_2"))

gene_exp <- read.table("gene_exp.diff", header = T)
top1000genes <- gene_exp[order(gene_exp$q_value)[1:100], c("gene")]
topGenes <- merge_final[merge_final$gene_short_name %in% top1000genes,]

# the original paper's blue to yellow range is narrower but there isn't quite 
# any easy way to arrange that magnitude
colfunc <- colorRampPalette(c("blue", "black", "yellow"), interpolate = "spline")
names(topGenes) <- c("gene_id", "gene_short_name", "Ad_1", "Ad_2", "P0_1", "P0_2",
                     "P4_1", "P4_2", "P7_1", "P7_2")
# plot two heatmaps to compare hierarchical clustering of samples
heatmap.2(data.matrix(topGenes[,3:10]), scale = "row", trace = "none", 
          na.color = "cyan", 
          col = colfunc(10),
          margins = c(5,2),
          cexCol = 0.9,
          xlab = "Sample name", ylab = "Gene ID",
          labRow = "",
          # key.xlab = "",
          lwid=c(1, 4),
          lhei=c(1, 4), 
          main = "Samples Clustered")

heatmap.2(data.matrix(topGenes[,3:10]), scale = "row", trace = "none", 
          na.color = "cyan", dendrogram = "row", Colv = F,
          col = colfunc(10),
          margins = c(5,2),
          cexCol = 0.9,
          xlab = "Sample name", ylab = "Gene ID",
          labRow = "",
          # key.xlab = "",
          lwid=c(1, 4),
          lhei=c(1, 4),
          main = "Samples Unclustered")


#### David table ####
# down <- read.csv("david_v1/downreg.csv")
# up <- read.csv("david_v1/upreg.csv")
down <- read.csv("downreg.csv")
up <- read.csv("upreg.csv")
down$Regulation <- "down"
up$Regulation <- "up"
top50 <- rbind(up, down)
dat1 <- data.frame(do.call(rbind, strsplit(as.vector(top50$Term), split = "~")))
names(dat1) <- c("GO_term", "GO_name")
top50 <- cbind(top50, dat1)
top50$lower <- tolower(top50$GO_name)

# manual curated online table iia
onlineTII <- read.csv("online_table_iia.csv")
onlineTII$lower1 <- tolower(onlineTII$term1)
onlineTII$lower1[which(onlineTII$lower1 == "i-band")] <- "i band"

onlineTII$lower1[which(onlineTII$lower1 == "z-disk")] <- "z disc"
# merge based on name
bigTable <- merge(top50, onlineTII, by.x = "lower", by.y = "lower1", all.x = T)

outTable <- bigTable[,c("GO_term", "GO_name", "Category", "PValue", "Fold.Enrichment",
                        "Bonferroni", "Benjamini", "Regulation", 
                        "term1", "FE", "BH")]
outTable <- outTable[order(outTable$Benjamini, decreasing = F),]
write.csv(outTable[1:25,], "regulated_table.csv", quote = T)
