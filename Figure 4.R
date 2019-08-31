library(ggplot2)
library(GSEABase)
library(GSVA)
library(gplots)
library(gdata)
library(geneplotter)

cor.dist <- function (x) 
{
  as.dist(1 - cor(t(x), use="pairwise.complete.obs"), method="kendall")
}
cor.dist.spearman <- function (x) 
{
  as.dist(1 - cor(t(x), use="pairwise.complete.obs"), method="spearman")
}
dist.manhattan <- function(x){dist(x, method="manhattan")}
hclust.ward.d <- function(x){hclust(x, method="ward.D")}
hclust.ward.d2 <- function(x){hclust(x, method="ward.D2")}

# Read variation-stablized RNA-seq data from a text file
vitro <- read.table("8cc-in-vitro.txt", header=TRUE, sep="\t", row.names=1)
vitro.ann <- data.frame(Cell=c(rep("2208L", 3), rep("4T1", 3), rep("67NR", 3), 
                               rep("AT3", 3), rep("E0771", 3), rep("PyMTE", 3), 
                               rep("PyMTM", 3), rep("T11",3)))

# Re-arrange sample order for presentation
vitro.ann$fig.3.plot.ind <- c(1,2,3,7,8,9,10,11,12,19,20,21,22,23,24,13,14,15,16,17,18,4,5,6)

# EMT-related genes curated from the literature
emt.cl.e <- c("Cdh1", "Epcam", "Cldn2", "Cldn4", 
              "Krt7", "Krt8", "Krt18", "Krt19", "Krt20", "Esrp1", "Mc1")
emt.cl.m <- c("Pdgfrb", "Zeb1","Fap", "Cdh2", "Cdh11", "Col1a1", 
              "Fn1", "Twist1", "Snai1", "Snai2")

vitro.e <- vitro[which(is.element(rownames(vitro), emt.cl.e)),]
vitro.m <- vitro[which(is.element(rownames(vitro), emt.cl.m)),]
vitro.em <- rbind(vitro.e, vitro.m)

# Generate Figure 4b
heatmap.2(vitro.em, col=c(rep("green", 30), greenred(100), rep("red",30)), dist=cor.dist,
          trace="none", scale="row", labCol=vitro.ann$Cell,  hclust=hclust.ward.d,
          cexRow=1.4, cexCol=1.4, margins=c(7,7), density.info="none")




# Read Hallmark Gene Sets, please obtain it from original authors' website: MSigDB
hm.sets <- getGmt("GSEA/Gene Sets/h.all.v6.0.symbols.gmt")

# Load TCGA dataset
load("TCGA-RNA-seq.RData")

# Isolate triple negative breast cancer
tcga.tn.ann <- tcga.ann[which(tcga.ann$er_status_by_ihc == "Negative" 
                              & tcga.ann$pr_status_by_ihc == "Negative" 
                              & tcga.ann$her2_status_by_ihc == "Negative"),]

tcga.tn <- tcga[,which(tcga.ann$er_status_by_ihc == "Negative" 
                       & tcga.ann$pr_status_by_ihc == "Negative" 
                       & tcga.ann$her2_status_by_ihc == "Negative")]

# Run GSVA on TCGA TN Breast Cancer data using Hallmark Gene Sets
tcga.tn.hm <- gsva(as.matrix(tcga.tn), hm.sets, method="gsva", rnaseq=TRUE)

# Focusing on mTOR and EMT pathways
tcga.tn.ann$mtor <- tcga.tn.hm$es.obs["HALLMARK_MTORC1_SIGNALING",]
tcga.tn.ann$mtor2 <- tcga.tn.hm$es.obs["HALLMARK_PI3K_AKT_MTOR_SIGNALING",]
tcga.tn.ann$emt.hm <- tcga.tn.hm$es.obs["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",]
tcga.tn <- as.matrix(tcga.tn)

# Read signatures from Taube et al., 2010, which are also provided in Supplemenatary Table.
taube.emt.up <- scan("emt-up.txt", what="")
taube.emt.dw <- scan("emt-dw.txt", what="")
tcga.tn.ann$emt.taube <- apply(tcga.tn[which(is.element(rownames(tcga.tn), emt.up)),], 2, sum) -
                         apply(tcga.tn[which(is.element(rownames(tcga.tn), emt.dw)),], 2, sum)
# Examine specific cytokines
tcga.tn.ann$IL8 <- tcga.tn["IL8",]
tcga.tn.ann$CXCL1 <- tcga.tn["CXCL1",]
tcga.tn.ann$CSF1 <- tcga.tn["CSF1",]
tcga.tn.ann$CCL2 <- tcga.tn["CCL2",]
tcga.tn.ann$TNFAIP6 <- tcga.tn["TNFAIP6",]

# Read TIMER output of TCGA dataset downloaded from the authors' website
timer <- read.table("Timer-TCGA-output.txt", header=TRUE, sep="\t")
# Reformat sample ID
timer$id2 <- substr(timer$ID, 1, 12)
timer$id2 <- sub("-", ".", timer$id2)
timer$id2 <- sub("-", ".", timer$id2)

# Isolate triple negative breast cancer
timer.tn <- as.matrix(timer[match(colnames(tcga.tn), timer$id2),2:7])
# Z-transformation
timer.tn.z <- apply(timer.tn, 2, function(x){(x-mean(x))/sd(x)})
# Timer scores of neutrophils and macrophages
timer.tn.ann$Neu <- timer.tn.z[,"Neutrophil"]
timer.tn.ann$Mac <- timer.tn.z[,"Macrophage"]

# Generate Figure 4g
fig.4f <- heatmap.2(t(as.matrix(tcga.tn.ann[,c(29:38)])), 
                    col=c(rep("blue",20), dChip.colors(100), rep("red",20)), trace="none",
                    scale=("row"), dist=cor.dist, hclust=hclust.ward.d)$colInd

heatmap.2(t(as.matrix(tcga.tn.ann[fig.4f,c(39:40)])), 
          col=c(dChip.colors(50)[1:50]), trace="none",
          dist=cor.dist.spearman, hclust=hclust.ward.d2, Colv=FALSE, margins=c(5,10))

