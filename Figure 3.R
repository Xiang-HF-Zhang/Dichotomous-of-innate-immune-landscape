library(gplots)
library(geneplotter)


# Read variation-stablized RNA-seq data from a text file
load("Cancer cells of 8 models.RData")
vitro.ann <- data.frame(Cell=c(rep("2208L", 3), rep("4T1", 3), rep("67NR", 3), 
                               rep("AT3", 3), rep("E0771", 3), rep("PyMTE", 3), 
                               rep("PyMTM", 3), rep("T11",3)))


# Cytokine genes curated from the literature
cytokines <- c("Csf3", "Cxcl1", "Cxcl2","Cxcl5","Il6", "Csf2", "Csf1", "Ccl2", "Tnfaip6")

vitro.ann$type <- c(rep("Neu", 6), rep("Mac", 3), rep("MN-hybrid", 3), rep("Mac", 3),
                    rep("Neu", 3), rep("Mac", 6))

vitro.cyto <- vitro[match(cytokines, rownames(vitro)),]

# Re-arrange sample order for presentation
vitro.ann$fig.3.plot.ind <- c(1,2,3,7,8,9,10,11,12,19,20,21,22,23,24,13,14,15,16,17,18,4,5,6)

# Generate Figure 3f
heatmap.2(vitro.cyto[,order(vitro.ann$fig.3.plot.ind)], col=dChip.colors(50), 
          trace="none", scale="row", labCol=vitro.ann$Cell[order(vitro.ann$fig.3.plot.ind)], 
          Rowv=FALSE, Colv=FALSE, colsep=c(6,12,18), rowsep=c(4,7))