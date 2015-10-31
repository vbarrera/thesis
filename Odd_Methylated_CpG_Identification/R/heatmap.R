# TF Analysis- Heatmap

setwd("~/projects/CGProfile/results/selectedPosAnalysis/integratedAnalysis/TFAnalysis/")
setwd("/databases/vbarrera/newSNP_check/TFAnalysis/s1t0/tfFiles/")
overlapMat<-read.table("overlapMatrixClean.csv",header=T,sep=";",row.names=1)
countTF<-apply(overlapMat,2,sum)
countTFSubSet<-countTF[countTF>=5]

write.table(sort(countTFSubSet,decreasing=T),file="countTFsubset_table.csv",sep=";",quote=F,col.names=F)
png("barplotTFOverlap.png")
barplot(sort(countTFSubSet,decreasing=T),cex.names=0.7,las=2)
dev.off()

overlapMatCleanTF<-overlapMat[,apply(overlapMat,2,sum)>0]

d.row = dist(overlapMatCleanTF, method = "binary")
hc.row = hclust(d.row, method="ward")

d.col = dist(t(overlapMatCleanTF), method = "binary")
hc.col = hclust(d.col, method="ward")

#heatmap

png("pngs/heatmap.png",height=1080,width=1920)
#pdf("pngs/heatmap.pdf")
heatmap(as.matrix(overlapMatCleanTF),Rowv=as.dendrogram(hc.row), Colv=as.dendrogram(hc.col), 
        scale='none',col=c("lightblue","red"),cexRow=0.4,cexCol=0.4)
dev.off()

png("heatmap_noClusterTF.jpeg",height=1080,width=1920)
heatmap(as.matrix(overlapMatCleanTF),Rowv=as.dendrogram(hc.row), Colv=F, 
        scale='none',col=c("lightblue","red"),cexRow=1,cexCol=1)
dev.off()