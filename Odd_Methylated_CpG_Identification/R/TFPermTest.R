# load package
library(regioneR)

# Region Sets
hg19_genome=getGenomeAndMask("hg19")
interestingCG<-toGRanges("input/interesting1t0PointsHg19.bed_forRegioneR")
allCG_Changing<-toGRanges("input/changingCG_H1_IMR90_NoX_Y.bed_forRegioneR")
tfAll<-toGRanges("input/wgEncodeRegTfbsClustered.bed_forRegioneR")
#tfMerged<-joinRegions(tfAll)


##
listTF<-unique(tfAll$TF)
resultsDF<-data.frame(TF=character(),pval=numeric(),z_score=numeric(),alternative=character())
for (A in listTF){
  tempGR<-subset(tfAll,tfAll$TF==A)
  joinedTempGR<-joinRegions(tempGR)
  pt<-permTest(A=interestingCG,ntimes=5000,randomize.function=resampleRegions,universe=allCG_Changing,
                   evaluate.function=numOverlaps,B=joinedTempGR,genome=hg19_genome$genome,mask=hg19_genome$mask)
  save(pt,file=paste("TFAnalysis/tfFiles/",A,".Rdata",sep=""))
  jpeg(paste("TFAnalysis/tfFiles/jpegs/",A,".jpeg",sep=""))
  plot(pt)
  dev.off()
  tempV<-data.frame(A,pt$pval,pt$zscore,pt$alternative)
  colnames(tempV)=c("TF","pval","z_score","alternative")
  resultsDF<-rbind(resultsDF,tempV)
  colnames(resultsDF)=c("TF","pval","z_score","alternative")
  write.table(resultsDF,file="TFAnalysis/tfFiles/resultsDF",sep=";",row.names=FALSE)
}

lz<-localZScore(A=interestingCG,pt=pt,window=10*mean(width(interestingCG)),step=mean(width(interestingCG))/2,
                B=joinedTempGR)
