# load package
library(regioneR)

# Region Sets
hg19_genome=getGenomeAndMask("hg19")
#interestingCG<-toGRanges("input/interesting1t0PointsHg19.bed_forRegioneR")
allCG_Changing<-toGRanges("input/changingCG_H1_IMR90_NoX_Y.bed_forRegioneR")
dnaseSites<-toGRanges("input/Simple_wgEncodeRegDnaseClusteredV3.bed_forRegioneR")
dnaseSitesMerged<-joinRegions(dnaseSites)

######
s1t0<-toGRanges("input/bedFiles_hg19/11111111111_11111011111_H1_IMR90_NoSNP.bed_hg19_forRegioneR")
s0t0<-toGRanges("input//bedFiles_hg19/11111011111_11111011111_H1_IMR90_NoSNP.bed_hg19_forRegioneR")
s0t1<-toGRanges("input/bedFiles_hg19/11111011111_11111111111_H1_IMR90_NoSNP.bed_hg19_forRegioneR")
s1t0_HC<-toGRanges("input//bedFiles_hg19/11111111111_11111011111_H1_IMR90_NoSNP_HighConf.bed_hg19_forRegioneR")
#####

#s1t0
interestingCG=s1t0
folder="s1t0"

pt<-permTest(A=interestingCG,ntimes=5000,randomize.function=resampleRegions,universe=allCG_Changing,
             evaluate.function=numOverlaps,B=dnaseSitesMerged,genome=hg19_genome$genome,mask=hg19_genome$mask)

save(pt,file=paste("dnaseAnalysis/",folder,"/dnasePermTest_changingNoX_Y_masked.Rdata",sep=""))
pt
png(paste("dnaseAnalysis/",folder,"/permutationTest.png",sep=""))
plot(pt)
dev.off()

#s0t0
interestingCG=s0t0
folder="s0t0"

pt<-permTest(A=interestingCG,ntimes=5000,randomize.function=resampleRegions,universe=allCG_Changing,
             evaluate.function=numOverlaps,B=dnaseSitesMerged,genome=hg19_genome$genome,mask=hg19_genome$mask)

save(pt,file=paste("dnaseAnalysis/",folder,"/dnasePermTest_changingNoX_Y_masked.Rdata",sep=""))
pt
png(paste("dnaseAnalysis/",folder,"/permutationTest.png",sep=""))
plot(pt)
dev.off()

#s0t1
interestingCG=s0t1
folder="s0t1"

pt<-permTest(A=interestingCG,ntimes=5000,randomize.function=resampleRegions,universe=allCG_Changing,
             evaluate.function=numOverlaps,B=dnaseSitesMerged,genome=hg19_genome$genome,mask=hg19_genome$mask)

save(pt,file=paste("dnaseAnalysis/",folder,"/dnasePermTest_changingNoX_Y_masked.Rdata",sep=""))
pt
png(paste("dnaseAnalysis/",folder,"/permutationTest.png",sep=""))
plot(pt)
dev.off()

#s1t0_HC
interestingCG=s1t0_HC
folder="s1t0_HC"

pt<-permTest(A=interestingCG,ntimes=5000,randomize.function=resampleRegions,universe=allCG_Changing,
             evaluate.function=numOverlaps,B=dnaseSitesMerged,genome=hg19_genome$genome,mask=hg19_genome$mask)

save(pt,file=paste("dnaseAnalysis/",folder,"/dnasePermTest_changingNoX_Y_masked.Rdata",sep=""))
pt
png(paste("dnaseAnalysis/",folder,"/permutationTest.png",sep=""))
plot(pt)
dev.off()
