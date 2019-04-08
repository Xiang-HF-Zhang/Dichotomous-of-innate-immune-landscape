library(gplots)
library(ggplot2)

# Define functions for later hierarchical clustering
cor.dist <- function (x) 
{
  as.dist(1 - cor(t(x), use="pairwise.complete.obs"))
}
dist.manhattan <- function(x){dist(x, method="manhattan")}
hclust.ward.d <- function(x){hclust(x, method="ward.D")}
hclust.ward.d2 <- function(x){hclust(x, method="ward.D2")}


# Read annotations of GSE78220 with timer output
mel.ann <- read.table("../Melanoma/gse78220-ann-wTIMER.txt", header=TRUE, sep="\t", row.names=1)

# Generate Figure 8h
heatmap.2(t(as.matrix(mel.ann[order(mel.ann$Neutrophil),c(10,11,7:9, 12)])), col=dChip.colors(100), 
          trace="none", scale="row",dist=dist.manhattan, hclust=hclust, Colv=FALSE, Rowv=FALSE,
          ColSideColors=ifelse(mel.ann$Response[order(mel.ann$Neutrophil)]=="Progressive Disease", "red", ifelse(mel.ann$Response[order(mel.ann$Neutrophil)]=="Partial Response", "orange", "green")), 
          margins=c(5,10))

# Analyze Riaz et al., 2016 data, please obtain original data and annotation from original authors
# Corrected some of the typos in original data, and group PR and SD together
riaz.ann <- read.table("riaz et al-ann.txt", header=TRUE, sep="\t", row.names=1)
riaz.ann$BOR[c(23:24,85)] <- c("PD", "PD", "PD")
riaz.ann$BOR <- as.vector(riaz.ann$BOR)
riaz.ann$myBOR <- ifelse(riaz.ann$BOR=="CR", "CR", 
                         ifelse(riaz.ann$BOR=="PD", "PD", "PR/SD"))
# Read timer output. 
# Note that Timer was run on viration-stablized and regularized log-transformed RNAseq data

riaz.timer <- read.table("riaz et al-timer.txt", header=TRUE, sep="\t", row.names=1)

# There is one extreme outlier in the dataset, which was removed. This did not affect any conclusion
# shown in Figure 7.
riaz.timer <- riaz.timer[-39,]
riaz.timer <- apply(riaz.timer, 2, function(x){(x-median(x))/sd(x)})
riaz.ann <- cbind(riaz.ann[-39,], as.data.frame(riaz.timer))

# Restrict to pre-treatment samples and generate Figure 8f
riaz.ann.pre <- riaz.ann[which(riaz.ann$Treatment=="Pre"),]
riaz.ann.pre$color <- ifelse(riaz.ann.pre$myBOR=="CR", "green", ifelse(riaz.ann.pre$myBOR=="PD", "red", "orange"))
heatmap.2(t(as.matrix(riaz.ann.pre[order(riaz.ann.pre$Neutrophil),c(15,16,12:14,17)])), trace="row",
          scale="none", Rowv=FALSE, Colv=FALSE,
          ColSideColors=riaz.ann.pre$color[order(riaz.ann.pre$Neutrophil)],
          col=dChip.colors(100))

# Generate Figure 8g
p <- ggplot(riaz.ann[which(riaz.ann$Treatment=="Pre"),], aes(myBOR, Neutrophil))
p + geom_boxplot(outlier.colour="white", fill=c("green", "red", "orange")) + geom_jitter(width=0.4, size=4, col="dark grey") + theme(panel.grid.major = element_blank(), axis.line.x = element_line(size=0.5, colour="black"), axis.line.y = element_line(size=0.5, colour="black"), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.text.x=element_text(colour="black", size=14), axis.text.y=element_text(colour="black", size=14), text=element_text(size=16))


# Combine the two melanoma datasets
composite.data <- data.frame(OS=c(mel.ann$OS, riaz.ann.pre$OS), 
                             Death=c(mel.ann$Death, riaz.ann.pre$OS_SOR), 
                             Neutrophil=c(mel.ann$Neutrophil.z, riaz.ann.pre$Neutrophil),
                             Response=c(mel.ann$Response, riaz.ann.pre$myBOR),
                             Dataset=c(rep("GSE78220", 28), rep("Riaz et al", 56)))

# Use top 20% as a cutoff to identify patients with extraordinarily high neutrophil content
# Data shown in Figure 8i
composite.data$Neu.bi <- ifelse(composite.data$Neutrophil > 0.44620026, "High", "Low")
fisher.test(composite.data$Neu.bi, composite.data$Dataset)
table(composite.data$Neu.bi, composite.data$Response)