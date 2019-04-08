setwd("Z:/Manuscripts in Preparation/Dichotomous of innate immune landscape/Nature Cell Biology/New analyses for revision")
tcga.m2h <- read.table("TCGA-mir200-data.txt", header=TRUE, sep="\t", row.names=1)
tcga.m2h.z <- apply(as.matrix(tcga.m2h[,4:8]), 2, function(x){(x-mean(x))/sd(x)})

tcga.tn.ann$m200a <- tcga.m2h.z[match(rownames(tcga.tn.ann), rownames(tcga.m2h)),1]
tcga.tn.ann$m200b <- tcga.m2h.z[match(rownames(tcga.tn.ann), rownames(tcga.m2h)),2]
tcga.tn.ann$m200c <- tcga.m2h.z[match(rownames(tcga.tn.ann), rownames(tcga.m2h)),3]
tcga.tn.ann$m141 <- tcga.m2h.z[match(rownames(tcga.tn.ann), rownames(tcga.m2h)),4]
tcga.tn.ann$m429 <- tcga.m2h.z[match(rownames(tcga.tn.ann), rownames(tcga.m2h)),5]
tcga.tn.ann$m200fam <- tcga.tn.ann$m200a + tcga.tn.ann$m200b + tcga.tn.ann$m200c + tcga.tn.ann$m141 +
                      tcga.tn.ann$m429

cor.test(tcga.tn.ann$m200fam, tcga.tn.ann$timer.mac, use="pairwise.complete.obs")

timer.tcga <- timer[match(rownames(tcga.ann), timer$id2), 2:7]
tcga.ann <- cbind(tcga.ann, timer.tcga)
m2h.tcga <- tcga.m2h.z[match(rownames(tcga.ann), rownames(tcga.m2h.z)),]
tcga.ann <- cbind(tcga.ann, m2h.tcga)

# Plot
library(ggplot2)
p5 <- ggplot(tcga.ann, aes(Macrophage, hsa.mir.200c)) 
p5  + geom_point(size=6, col="light grey") + theme(panel.grid.major = element_blank(), axis.line.x = element_line(size=0.5, colour="black"), axis.line.y = element_line(size=0.5, colour="black"), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.text.x=element_text(colour="black", size=14), axis.text.y=element_text(colour="black", size=14), text=element_text(size=16))

tcga.ann$Macrophag.cat <- ifelse(tcga.ann$Macrophage == 0, "None",
                                 ifelse(tcga.ann$Macrophage < 0.3, "Low",
                                        "High"))
tcga.ann$TN <- ifelse(tcga.ann$er_status_by_ihc=="Negative" & tcga.ann$pr_status_by_ihc=="Negative" &
                        tcga.ann$her2_status_by_ihc=="Negative", "TN", "Non-TN")


boxplot(hsa.mir.200c~Macrophag.cat, tcga.ann[which(tcga.ann$TN=="TN"),])
mir200c.aov <- aov(hsa.mir.200c~Macrophag.cat, tcga.ann[which(tcga.ann$TN=="TN"),])
summary(mir200c.aov)

table(tcga.ann$Macrophag.cat,tcga.ann$TN)
