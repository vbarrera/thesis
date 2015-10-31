# load package
library(regioneR)

'''# Region Sets
interestingCG<-toGRanges("interesting1t0PointsHg19.bed_tmp2")
allCG<-toGRanges("methList_CGi_hg19_H1_CG_allGenome.bed_tmp")
allCG_Changing<-toGRanges("changingCG_H1_IMR90.bed")
dnaseSites<-toGRanges("Simple_wgEncodeRegDnaseClusteredV3.bed")
tfAll<-toGRanges("wgEncodeRegTfbsClustered.bed")
tfMerged<-joinRegions(tfAll)
hg19_genome=getGenomeAndMask("hg19")
##'''


# Region Sets
hg19_genome=getGenomeAndMask("hg19")
allCG_Changing<-toGRanges("input/changingCG_H1_IMR90_NoX_Y.bed_forRegioneR")
tfAll<-toGRanges("input/wgEncodeRegTfbsClustered.bed_forRegioneR")
tfMerged<-joinRegions(tfAll)
######
s1t0<-toGRanges("input/bedFiles_hg19/11111111111_11111011111_H1_IMR90_NoSNP.bed_hg19_forRegioneR")
s0t0<-toGRanges("input/bedFiles_hg19/11111011111_11111011111_H1_IMR90_NoSNP.bed_hg19_forRegioneR")
s0t1<-toGRanges("input/bedFiles_hg19/11111011111_11111111111_H1_IMR90_NoSNP.bed_hg19_forRegioneR")
s1t0_HC<-toGRanges("input/bedFiles_hg19/11111111111_11111011111_H1_IMR90_NoSNP_HighConf.bed_hg19_forRegioneR")
#####

#s1t0
interestingCG=s1t0_HC
folder="s1t0_HC"

##

listTF<-unique(tfAll$TF)
resultsDF<-data.frame(ID=character())
countsTF<-data.frame(TF=character(),counts=integer())
for (A in listTF){
  tempGR<-subset(tfAll,tfAll$TF==A)
  joinedTempGR<-joinRegions(tempGR)
  overlapDF<-overlapRegions(interestingCG,joinedTempGR)
  numOL<-numOverlaps(interestingCG,joinedTempGR)
  ID<-paste(overlapDF$chr,overlapDF$startA,sep="_")
  TF<-rep("1",numOL)
  tempTFOverlapDF<-cbind(ID,TF)
  colnames(tempTFOverlapDF)<-c("ID",A)
  resultsDF<-merge(resultsDF,tempTFOverlapDF,all=TRUE)
  tempCountsDF<-data.frame(TF=character(),counts=integer())
  tempCountsDF<-rbind(tempCountsDF,c(A,numOL))
  colnames(tempCountsDF)<-c("TF","counts")
  countsTF<-rbind(countsTF,tempCountsDF)
}

countsTF<-data.frame(TF=as.character(countsTF$TF),counts=as.numeric(as.vector(countsTF$counts)))
countsTF<-countsTF[order(-countsTF$counts),]
countsTF<-data.frame(TF=as.character(countsTF$TF),counts=as.numeric(as.vector(countsTF$counts)))
significativeTF<-subset(countsTF,countsTF$counts>5)

barplot(significativeTF$counts,names.arg=significativeTF$TF,cex.names=0.7)


write.table(resultsDF,file=paste0("TFAnalysis/",folder,"/tfFiles/overlapMatrix.csv"),na="0",sep=";")
write.table(resultsDF,file=paste0("TFAnalysis/",folder,"/tfFiles/overlapMatrixClean.csv"),na="0",sep=";",quote=F,row.names=F)
write.table(t(resultsDF),file=paste0("TFAnalysis/",folder,"/tfFiles/overlapMatrixTransposed.csv"),na="0",sep=";")
write.table(countsTF,file=paste0("TFAnalysis/",folder,"/tfFiles/countsTF.csv"),sep=";")

countsTF<-read.table(paste0("TFAnalysis/",folder,"/tfFiles/overlapMatrix.csv"),sep=";")

vectorNames<-character()
for (i in 1:nrow(countsTF)){
  overlapString=paste(as.character(countsTF[i,1]),";",sep="")
  for (j in 2:ncol(countsTF)){
    if (countsTF[i,j]==1){
      overlapString=paste(overlapString,colnames(countsTF)[j],sep=",")
    }
  }
  vectorNames[i]<-overlapString
}

write.table(vectorNames,file=paste0("TFAnalysis/",folder,"/tfFiles/OverlapVectorNames.txt"),row.names=F,quote=F,col.names=F)

######