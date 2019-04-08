library(ggplot2)
library(GSEABase)
library(GSVA)
library(gplots)
library(gdata)

############################################################
# Load RNAseq data after variance stabilizing transformation
############################################################

# Load data
load("Macrophage-neutrophil-RNAseq.RData")
# Mapping human and mouse orthologous genes
human.ortho <- convertMouseGeneList(rownames(mye))
# Collapse data to the same genes
mye$human <- human.ortho[match(rownames(mye), human.ortho$MGI.symbol),2]
mye.human <- max.redundant.probes(mye, 26, 1:24)
# Log transformation and quantile normalization
mye.human <- log2(mye.human)
mye.human.qn <- normalize.quantiles(mye.human)
rownames(mye.human.qn) <- rownames(mye.human)
colnames(mye.human.qn) <- colnames(mye.human)

###########################################################
#Figure 5b: GSVA of Hallmark genesets on macrophage dataset
###########################################################

# Read gene sets
hm.sets <- getGmt("h.all.v6.0.symbols.gmt") 
# GSVA implementation on macrophage profiles
mac.hm <- gsva(as.matrix(mye.human[,c(1:3, 10:15, 19:21)]), hm.sets, method="gsva") 
# Data presentation
heatmap.2(mac.hm$es.obs, col=c(rep("green", 1), greenred(100), rep("red", 1)), 
          trace="none",  margins=c(10,25), density.info="none", Colv=FALSE,
          dist=cor.dist, hclust=hclust)

###########################################################
#Figure 6b: GSVA of Hallmark genesets on neutrophil dataset
###########################################################


neu.hm <- gsva(as.matrix(mye.human[,c(4:9, 16:18, 22:24)]), hm.sets, method="gsva")
heatmap.2(neu.hm$es.obs, col=c(rep("green", 1), greenred(100), rep("red", 1)), 
          trace="none", margins=c(10,25), density.info="none", Colv=FALSE,
          dist=cor.dist, hclust=hclust.ward.d)

###########################################################
#Figure 6a: Hiearchical clustering of tumor-associated neu-
#trophils with published MDSC and TAN signatures.
###########################################################

# Load data downloaded from Fridlender, Z. G. et al.2012
load("TAN.RData")
# Collapse probes of the same gene into one value - take the maximum
tan.df <- as.data.frame(tan.dt)
tan.df$Symbol <- tan.probes$Symbol
tan.dt.clps <- max.redundant.probes(tan.df, 16, 1:15)

# Use SAM to analyze the data and identify genes up- and down-regulated in 
# MDSCs and TANs. Cutoff is set to be FC > 10 and FDR < 0.05.
tan.sam <- sam(tan.dt.clps, cl=c(rep(0,7), rep(1,8)))
tan.up <- names(tan.sam@fold)[which(tan.sam@fold > 10 & tan.sam@q.value < 0.05)]
tan.dw <- names(tan.sam@fold)[which(tan.sam@fold < 0.1 & tan.sam@q.value < 0.05)]

# Use the identified genes to cluster neutrophil data.

tan.up.cm <- intersect(tan.up, rownames(mye))
tan.dw.cm <- intersect(tan.dw, rownames(mye))
tan.diff <- c(tan.up.cm, tan.dw.cm)
neu.tan <- mye.qn[match(tan.diff, rownames(mye.qn)), c(4:9, 16:18, 22:24)]
heatmap.2(t(as.matrix(neu.tan)), scale="col", trace="none", 
          col=c(rep("green", 1), greenred(100), rep("red", 1)), Rowv=TRUE,
          ColSideColors=ifelse(is.element(rownames(neu.tan), tan.up.cm), "orange", "yellow"),
          margins=c(10,20), Colv=TRUE, hclust=hclust.ward.d, dist=cor.dist.kendall)