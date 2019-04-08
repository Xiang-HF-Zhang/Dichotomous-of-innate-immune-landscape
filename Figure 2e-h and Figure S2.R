library(ggplot2)
library(gplots)
library(gdata)
library(Rtsne)

# Define functions for later hierarchical clustering
cor.dist <- function (x) 
{
  as.dist(1 - cor(t(x), use="pairwise.complete.obs"))
}
dist.manhattan <- function(x){dist(x, method="manhattan")}
hclust.ward.d <- function(x){hclust(x, method="ward.D")}
hclust.ward.d2 <- function(x){hclust(x, method="ward.D2")}

# TCGA breast cancer dataset
load("TCGA-RNA-seq.RData")

# Isolate triple negative breast cancer
tcga.tn.ann <- tcga.ann[which(tcga.ann$er_status_by_ihc == "Negative" 
                              & tcga.ann$pr_status_by_ihc == "Negative" 
                              & tcga.ann$her2_status_by_ihc == "Negative"),]

tcga.tn <- tcga[,which(tcga.ann$er_status_by_ihc == "Negative" 
                       & tcga.ann$pr_status_by_ihc == "Negative" 
                       & tcga.ann$her2_status_by_ihc == "Negative")]

# Load METABRIC Dataset. Please contact the original author to obtain the data
load("METABRIC.RData")

# Isolate triple negative breast cancer

mb.d.tn <- metabric.d.exp.clps[,which(metabric.d.ann$ER.Expr=="-" 
                                      & metabric.d.ann$PR.Expr=="-" 
                                      & metabric.d.ann$Her2.Expr=="-")]
mb.v.tn <- metabric.v.exp.clps[,which(metabric.v.ann$ER.Expr=="-" 
                                      & metabric.v.ann$PR.Expr=="-" 
                                      & metabric.v.ann$Her2.Expr=="-")]
mb.d.tn.ann <- metabric.d.ann[which(metabric.d.ann$ER.Expr=="-" 
                                    & metabric.d.ann$PR.Expr=="-" 
                                    & metabric.d.ann$Her2.Expr=="-"),]
mb.v.tn.ann <- metabric.v.ann[which(metabric.v.ann$ER.Expr=="-" 
                                    & metabric.v.ann$PR.Expr=="-" 
                                    & metabric.v.ann$Her2.Expr=="-"),]
# Combining discovery and validation subsets
mb.tn <- cbind(mb.d.tn, mb.v.tn)
mb.tn.ann <- rbind(mb.d.tn.ann, mb.v.tn.ann)


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
# Generating Figure 2e
timer.tn.hclust <- heatmap.2(t(timer.tn), col=c(dChip.colors(50)), dist=dist.manhattan,
                             hclust=hclust.ward.d2, scale="row", margins=c(5,10), trace="row", tracecol="black")
# Isolate cold and hot clusters
timer.tn.cutree <- cutree(hclust(dist(timer.tn, method="manhattan"), method="ward.D2"), 2)
# Calculate "ratio" between neutrophils and macrophages
timer.tn.mnratio <- timer.tn.z[,"Neutrophil"] - timer.tn.z[, "Macrophage"]
# Set the ratio to "NA" if the tumors are cold.
timer.tn.mnratio <- ifelse(timer.tn.cutree==1, NA, timer.tn.mnratio)
# Execute t-SNE
timer.tsne <- Rtsne(as.matrix(timer.tn.z), initial_dims = 6, perplexity = 5)
# Generating Figure 2f
timer.plot <- data.frame(X=timer.tsne$Y[,1], Y=timer.tsne$Y[,2], mnratio=timer.tn.mnratio)
p7 <- ggplot(timer.plot, aes(X, Y, mnratio)) 
p7  + geom_point(size=6, aes(colour=mnratio)) + scale_colour_gradient2(low="blue", mid="light grey", high="orange") + theme(panel.grid.major = element_blank(), axis.line.x = element_line(size=0.5, colour="black"), axis.line.y = element_line(size=0.5, colour="black"), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.text.x=element_text(colour="black", size=14), axis.text.y=element_text(colour="black", size=14), text=element_text(size=16))

# Generating Figure 2g and 2h for non-triple negative breast cancers in similar procedures as above
timer.nontn <- timer.allbc[-which(is.element(timer.allbc[,1], timer.tn[,1])),]
timer.nontn.z <- apply(timer.nontn, 2, function(x){(x-mean(x))/sd(x)})
heatmap.2(t(timer.nontn), col=c(rep("blue", 50), dChip.colors(50), rep("red", 50)), dist=dist.manhattan,
          hclust=hclust.ward.d2, scale="row", margins=c(5,10), trace="row", tracecol="black")
timer.nontn.tsne <- Rtsne(as.matrix(timer.nontn.z), initial_dims = 6, perplexity = 5)
timer.nontn.mnratio <- timer.nontn.z[,"Neutrophil"] - timer.nontn.z[, "Macrophage"]
timer.nontn.hclust <- heatmap.2(t(timer.nontn), col=c(rep("blue", 50), dChip.colors(50), rep("red", 50)), dist=dist.manhattan,
                                hclust=hclust.ward.d2, scale="row", margins=c(5,10), trace="row", tracecol="black")timer.tn.cutree <- cutree(hclust(dist(timer.tn, method="manhattan"), method="ward.D2"), 2)
timer.nontn.cutree <- cutree(hclust(dist(timer.nontn, method="manhattan"), method="ward.D2"), 2)
timer.nontn.mnratio <- ifelse(timer.nontn.cutree==1, NA, timer.nontn.mnratio)

timer.nontn.plot <- data.frame(X=timer.nontn.tsne$Y[,1], Y=timer.nontn.tsne$Y[,2], mnratio=timer.nontn.mnratio)

p7 <- ggplot(timer.nontn.plot, aes(X, Y, mnratio)) 
p7  + geom_point(size=3, aes(colour=mnratio)) + scale_colour_gradient2(low="blue", mid="light grey", high="orange") + theme(panel.grid.major = element_blank(), axis.line.x = element_line(size=0.5, colour="black"), axis.line.y = element_line(size=0.5, colour="black"), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.text.x=element_text(colour="black", size=14), axis.text.y=element_text(colour="black", size=14), text=element_text(size=16))


# Repeat similar procedures on METABRIC and EMC-MSK datasets to generate Supplementary Figure 2c,d.
mb.timer.output <- read.table("metabric-timer-output", sep="\t", row.names=1, header=TRUE)
mb.tn.timer <- as.matrix(mb.timer.output[,1:6])
mb.tn.timer.z <- apply(mb.timer.output, 2, function(x){(x-mean(x))/sd(x)})
mb.tn.mnratio <- mb.tn.timer.z[,"Neutrophil"] - mb.tn.timer.z[,"Macrophage"]
mb.tn.total <- apply(mb.tn.timer.z, 1, sum)
# Here a different approach was used to define "cold" cluster.
mb.tn.mnratio <- ifelse(mb.tn.total < 0, NA, mb.tn.mnratio)
mb.tn.rtsne <- Rtsne(as.matrix(mb.tn.timer.z), initial_dims = 6, perplexity = 5, theta=0, max_iter=10000)

mb.tn.plot <- data.frame(X=mb.tn.rtsne$Y[,1], Y=mb.tn.rtsne$Y[,2], mnratio=mb.tn.mnratio)

p7 <- ggplot(mb.tn.plot, aes(X, Y, mnratio)) 
p7  + geom_point(size=3, aes(colour=mnratio)) + scale_colour_gradient2(low="blue", mid="light grey", high="orange") + theme(panel.grid.major = element_blank(), axis.line.x = element_line(size=0.5, colour="black"), axis.line.y = element_line(size=0.5, colour="black"), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.text.x=element_text(colour="black", size=14), axis.text.y=element_text(colour="black", size=14), text=element_text(size=16))

