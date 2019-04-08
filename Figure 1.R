library(gplots)
library(ggplot2)

# Load data
imm.log <- read.table("immune-profiles-in-18-models.txt")

# Generating Figure 1e
heatmap.2(t(as.matrix(imm.log[,3:7])), col=c(rep("green", 1), greenred(100), rep("red",1)), 
          trace="none", scale="row", dist=dist.manhattan)

# Generating Figure 1f
p7 <- ggplot(imm.log, aes(Tu.Neu-Tu.Mac, Tu.CD45))
p7  + geom_point(size=6, aes(col=Line, shape=Background)) + geom_smooth(method=lm) + theme(panel.grid.major = element_blank(), axis.line.x = element_line(size=0.5, colour="black"), axis.line.y = element_line(size=0.5, colour="black"), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.text.x=element_text(colour="black", size=14), axis.text.y=element_text(colour="black", size=14), text=element_text(size=16))

# Generating Figure S2d
imm.hclust <- hclust(dist(as.matrix(imm.log[,3:7]), method="manhattan"))
imm.cutree <- cutree(imm.hclust, 4)
imm.log$cluster=imm.cutree
p <- ggplot(imm.log, aes(cluster, Tu.CD3.CD4+Tu.CD3.CD8))
p + geom_boxplot(outlier.colour="white", aes(group=cluster)) + geom_jitter(width=0.2, size=4, col="dark grey") + theme(panel.grid.major = element_blank(), axis.line.x = element_line(size=0.5, colour="black"), axis.line.y = element_line(size=0.5, colour="black"), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.text.x=element_text(colour="black", size=14), axis.text.y=element_text(colour="black", size=14), text=element_text(size=16))

# Generate Figure 1g-1j
p7 <- ggplot(imm.log, aes(Tu.Neu, Bl.Neu)) 
p7  + geom_point(size=6, aes(col=Line, shape=Background)) + geom_smooth(method=lm) + theme(panel.grid.major = element_blank(), axis.line.x = element_line(size=0.5, colour="black"), axis.line.y = element_line(size=0.5, colour="black"), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.text.x=element_text(colour="black", size=14), axis.text.y=element_text(colour="black", size=14), text=element_text(size=16))

p7 <- ggplot(imm.log, aes(Tu.Mac, Bl.Mono)) 
p7  + geom_point(size=6, aes(col=Line, shape=Background)) + geom_smooth(method=lm) + theme(panel.grid.major = element_blank(), axis.line.x = element_line(size=0.5, colour="black"), axis.line.y = element_line(size=0.5, colour="black"), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.text.x=element_text(colour="black", size=14), axis.text.y=element_text(colour="black", size=14), text=element_text(size=16))

p7 <- ggplot(imm.log, aes(Tu.mono, Bl.Mono)) 
p7  + geom_point(size=6, aes(col=Line, shape=Background)) + geom_smooth(method=lm) + theme(panel.grid.major = element_blank(), axis.line.x = element_line(size=0.5, colour="black"), axis.line.y = element_line(size=0.5, colour="black"), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.text.x=element_text(colour="black", size=14), axis.text.y=element_text(colour="black", size=14), text=element_text(size=16))

p7 <- ggplot(imm.log, aes(Tu.Mac, Tu.Mono)) 
p7  + geom_point(size=6, aes(col=Line, shape=Background)) + geom_smooth(method=lm) + theme(panel.grid.major = element_blank(), axis.line.x = element_line(size=0.5, colour="black"), axis.line.y = element_line(size=0.5, colour="black"), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.text.x=element_text(colour="black", size=14), axis.text.y=element_text(colour="black", size=14), text=element_text(size=16))
