# load package
library(regioneR)

# Region Sets
hg19_genome=getGenomeAndMask("hg19")
#interestingCG<-toGRanges("input/interesting1t0PointsHg19.bed_forRegioneR")
allCG_Changing<-toGRanges("input/changingCG_H1_IMR90_NoX_Y.bed_forRegioneR")
#H1_HMM<-toGRanges("input/H1_HMM_States.bed_forRegioneR")
H1_HMM<-toGRanges("input/H1_HMM_States_WithGroups.bed_forRegioneR")

######
s1t0<-toGRanges("input/bedFiles_hg19/11111111111_11111011111_H1_IMR90_NoSNP.bed_hg19_forRegioneR")
s0t0<-toGRanges("input/bedFiles_hg19/11111011111_11111011111_H1_IMR90_NoSNP.bed_hg19_forRegioneR")
s0t1<-toGRanges("input/bedFiles_hg19/11111011111_11111111111_H1_IMR90_NoSNP.bed_hg19_forRegioneR")
s1t0_HC<-toGRanges("input/bedFiles_hg19/11111111111_11111011111_H1_IMR90_NoSNP_HighConf.bed_hg19_forRegioneR")
#####

#s1t0
interestingCG=s1t0
folder="s1t0"

##
listStates<-unique(H1_HMM$minorGroup)
resultsDF<-data.frame(State=character(),pval=numeric(),z_score=numeric(),alternative=character())
for (A in listStates){
  tempGR<-subset(H1_HMM,H1_HMM$minorGroup==A)
  joinedTempGR<-joinRegions(tempGR)
  pt<-permTest(A=interestingCG,ntimes=5000,randomize.function=resampleRegions,universe=allCG_Changing,
               evaluate.function=numOverlaps,B=joinedTempGR,genome=hg19_genome$genome,mask=hg19_genome$mask)
  save(pt,file=paste("chromatinStates/minorGroup/",folder,"/statesFiles/",A,".Rdata",sep=""))
  png(paste("chromatinStates/minorGroup/",folder,"/statesFiles/pngs/",A,".png",sep=""))
  plot(pt)
  dev.off()
  tempV<-data.frame(A,pt$pval,pt$zscore,pt$alternative)
  colnames(tempV)=c("State","pval","z_score","alternative")
  resultsDF<-rbind(resultsDF,tempV)
  colnames(resultsDF)=c("State","pval","z_score","alternative")
  write.table(resultsDF,file=paste("chromatinStates/minorGroup/",folder,"/statesFiles/resultsDF",sep=""),sep=";",row.names=FALSE)
}

#s0t0
interestingCG=s0t0
folder="s0t0"

##
listStates<-unique(H1_HMM$minorGroup)
resultsDF<-data.frame(State=character(),pval=numeric(),z_score=numeric(),alternative=character())
for (A in listStates){
  tempGR<-subset(H1_HMM,H1_HMM$minorGroup==A)
  joinedTempGR<-joinRegions(tempGR)
  pt<-permTest(A=interestingCG,ntimes=5000,randomize.function=resampleRegions,universe=allCG_Changing,
               evaluate.function=numOverlaps,B=joinedTempGR,genome=hg19_genome$genome,mask=hg19_genome$mask)
  save(pt,file=paste("chromatinStates/minorGroup/",folder,"/statesFiles/",A,".Rdata",sep=""))
  png(paste("chromatinStates/minorGroup/",folder,"/statesFiles/pngs/",A,".png",sep=""))
  plot(pt)
  dev.off()
  tempV<-data.frame(A,pt$pval,pt$zscore,pt$alternative)
  colnames(tempV)=c("State","pval","z_score","alternative")
  resultsDF<-rbind(resultsDF,tempV)
  colnames(resultsDF)=c("State","pval","z_score","alternative")
  write.table(resultsDF,file=paste("chromatinStates/minorGroup/",folder,"/statesFiles/resultsDF",sep=""),sep=";",row.names=FALSE)
}

#s0t1
interestingCG=s0t1
folder="s0t1"

##
listStates<-unique(H1_HMM$minorGroup)
resultsDF<-data.frame(State=character(),pval=numeric(),z_score=numeric(),alternative=character())
for (A in listStates){
  tempGR<-subset(H1_HMM,H1_HMM$minorGroup==A)
  joinedTempGR<-joinRegions(tempGR)
  pt<-permTest(A=interestingCG,ntimes=5000,randomize.function=resampleRegions,universe=allCG_Changing,
               evaluate.function=numOverlaps,B=joinedTempGR,genome=hg19_genome$genome,mask=hg19_genome$mask)
  save(pt,file=paste("chromatinStates/minorGroup/",folder,"/statesFiles/",A,".Rdata",sep=""))
  png(paste("chromatinStates/minorGroup/",folder,"/statesFiles/pngs/",A,".png",sep=""))
  plot(pt)
  dev.off()
  tempV<-data.frame(A,pt$pval,pt$zscore,pt$alternative)
  colnames(tempV)=c("State","pval","z_score","alternative")
  resultsDF<-rbind(resultsDF,tempV)
  colnames(resultsDF)=c("State","pval","z_score","alternative")
  write.table(resultsDF,file=paste("chromatinStates/minorGroup/",folder,"/statesFiles/resultsDF",sep=""),sep=";",row.names=FALSE)
}

#s1t0_HC
interestingCG=s1t0_HC
folder="s1t0_HC"

##
listStates<-unique(H1_HMM$minorGroup)
resultsDF<-data.frame(State=character(),pval=numeric(),z_score=numeric(),alternative=character())
for (A in listStates){
  tempGR<-subset(H1_HMM,H1_HMM$minorGroup==A)
  joinedTempGR<-joinRegions(tempGR)
  pt<-permTest(A=interestingCG,ntimes=5000,randomize.function=resampleRegions,universe=allCG_Changing,
               evaluate.function=numOverlaps,B=joinedTempGR,genome=hg19_genome$genome,mask=hg19_genome$mask)
  save(pt,file=paste("chromatinStates/minorGroup/",folder,"/statesFiles/",A,".Rdata",sep=""))
  png(paste("chromatinStates/minorGroup/",folder,"/statesFiles/pngs/",A,".png",sep=""))
  plot(pt)
  dev.off()
  tempV<-data.frame(A,pt$pval,pt$zscore,pt$alternative)
  colnames(tempV)=c("State","pval","z_score","alternative")
  resultsDF<-rbind(resultsDF,tempV)
  colnames(resultsDF)=c("State","pval","z_score","alternative")
  write.table(resultsDF,file=paste("chromatinStates/minorGroup/",folder,"/statesFiles/resultsDF",sep=""),sep=";",row.names=FALSE)
}