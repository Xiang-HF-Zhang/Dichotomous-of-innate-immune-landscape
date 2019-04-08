library(siggenes)
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


# Data loading from BioGPS primary cell atlas. Please obtain the data from original authors.
imm.atlas <- load("BioGPS Atlas.RData")

# Identify neutrophils and macrophages cell indices
neu.ind <- c(699:716)
mac.ind <- c(356:505)

# Perform SAM to identify differnetially expressed genes as schematically shown in Figure 2a.
mvn <- sam(imm.atlas[, c(mac.ind, neu.ind)], cl=c(rep(1, length(mac.ind)), rep(0, length(neu.ind))))
fc.cut <- 8
q.cut <- 0.1
mac.vs.neu <- names(mvn@fold)[which(mvn@fold > fc.cut & mvn@q.value < q.cut)]
neu.vs.mac <- names(mvn@fold)[which(mvn@fold < 1/fc.cut & mvn@q.value < q.cut)]

# Load Nanostring data of 72 TNBC
load("Nanostring-of-72-TNBC-with-IHC-info.RData")

# Hiearchical clustering shown in Figure 2b
saw.mac <- saw[which(is.element(rownames(saw), mac.vs.neu)),]
saw.neu <- saw[which(is.element(rownames(saw), neu.vs.mac)),]
saw.mn <- rbind(saw.mac, saw.neu)

heatmap.2(saw.mn, col=c(rep(greenred(100)[1],30), greenred(100), rep(greenred(100)[100],30)), trace="none",
          scale="row", RowSideColors=c(rep("blue", dim(saw.mac)[1]), rep("yellow", dim(saw.neu)[1])),
          Rowv=TRUE, hclust=hclust.ward.d, dist=dist.manhattan)

# To generate boxplots on top of the heatmap in Figure 2b
saw.hclust <- hclust.ward.d(dist.manhattan(t(saw.mn)))
cluster.names <- c("Mac", "Int-1", "Neu", "Int-2")
saw.ann$mn.cutree <- cluster.names[cutree(saw.hclust,4)]

p <- ggplot(saw.ann, aes(mn.cutree, log2(TIL.CD68.count+1), fill=mn.cutree))
p + geom_boxplot(outlier.colour="white") + geom_jitter(width=0.4, size=3) + theme_minimal()

# To generate the separte heatmaps for selected markers shown at the bottom of Figure 2b
ind.markers <- c("ELANE", "CSF3", "CD68")
heatmap.2(saw[match(ind.markers, rownames(saw)), mn.heatmap$colInd], 
          col=c(rep(greenred(100)[1],30), greenred(100), rep(greenred(100)[100],30)), trace="none",
          scale="row", Colv=FALSE, dist=dist.manhattan, hclust=hclust.ward.d, margins=c(10,5))

# To generate boxplots in Figure 2c
cd68.temp <- saw.ann$TIL.CD68.count
saw.ann$TIL.CD68.cat <- ifelse(cd68.temp < 10, "Negative or Weak", "Strong")
t.test(mac~TIL.CD68.cat, saw.ann)
p <- ggplot(saw.ann, aes(TIL.CD68.cat, mac, fill=TIL.CD68.cat))
p + geom_boxplot(outlier.colour="white") + geom_jitter(width=0.4, size=3) + theme_minimal()

# To generate the histogram shown in Figure 2d
hist(apply(saw.neu, 2, sum) - apply(saw.mac, 2, sum), breaks=25)
