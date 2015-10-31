################
### Figure 31 ###
################

#################################
## With Informativeness filter ##
#################################
selected_reDinucleotide_Element_ADS<-subset(reDinucleotide_Element_ADS,
                                           reDinucleotide_Element_ADS$meanMethCoef_ADS!=-1 & 
                                             reDinucleotide_Element_ADS$E_methCoef!=-1 & ((reDinucleotide_Element_ADS$E_posInf/2)/reDinucleotide_Element_ADS$E_nCG)>=0.25)
selected_reDinucleotide_Element_ADS_Adipose<-subset(reDinucleotide_Element_ADS_Adipose,
                                              reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose!=-1 & 
                                                reDinucleotide_Element_ADS_Adipose$E_methCoef!=-1 & ((reDinucleotide_Element_ADS_Adipose$E_posInf/2)/reDinucleotide_Element_ADS_Adipose$E_nCG)>=0.25)
selected_reDinucleotide_Element_ADS_IPSC<-subset(reDinucleotide_Element_ADS_IPSC,
                                                    reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC!=-1 & 
                                                      reDinucleotide_Element_ADS_IPSC$E_methCoef!=-1 & ((reDinucleotide_Element_ADS_IPSC$E_posInf/2)/reDinucleotide_Element_ADS_IPSC$E_nCG)>=0.25)


##########################

#library(RColorBrewer)
#cols<-brewer.pal(3,"Dark2")
reducedDataADS<-selected_reDinucleotide_Element_ADS[,c(1,15,7,12)]
reducedDataADS_Adipose<-selected_reDinucleotide_Element_ADS_Adipose[,c(1,15,7,12)]
reducedDataADS_IPSC<-selected_reDinucleotide_Element_ADS_IPSC[,c(1,15,7,12)]



meanValuesDIR<-paste(Sys.getenv("HOME"),"/projects/Evaluation_of_single_CpG_sites_as_proxies_of_CpG_island_methylation_states_at_the_genome_scale/results/20120802/",sep="")

meanValuesADS<-read.table(paste(meanValuesDIR,"meanValueADS_Filtered.txt",sep=""),header=FALSE,sep="\t")
meanValuesADS_Adipose<-read.table(paste(meanValuesDIR,"meanValueADS_Adipose_Filtered.txt",sep=""),header=FALSE,sep="\t")
meanValuesADS_IPSC<-read.table(paste(meanValuesDIR,"meanValueADS_IPSC_Filtered.txt",sep=""),header=FALSE,sep="\t")


colnames(meanValuesADS)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")
colnames(meanValuesADS_Adipose)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")
colnames(meanValuesADS_IPSC)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")



mergedReducedMeanADS<-merge(reducedDataADS,meanValuesADS,by="id_In_Type",all=TRUE)
mergedReducedMeanADS_Adipose<-merge(reducedDataADS_Adipose,meanValuesADS_Adipose,by="id_In_Type",all=TRUE)
mergedReducedMeanADS_IPSC<-merge(reducedDataADS_IPSC,meanValuesADS_IPSC,by="id_In_Type",all=TRUE)



xR<-subset(mergedReducedMeanADS,mergedReducedMeanADS$nREs>1)
xJ<-subset(mergedReducedMeanADS_Adipose,mergedReducedMeanADS_Adipose$nREs>1)
xL<-subset(mergedReducedMeanADS_IPSC,mergedReducedMeanADS_IPSC$nREs>1)

### ADS
######
png(paste(resultsDIR,"figure31_ADSFilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist((xR$E_methCoef-xR$meanMethCoef_ADS), breaks=seq(-1,1,0.1), plot=FALSE)
yhist <- hist((xR$E_methCoef-xR$meanREMethCoef), breaks=seq(-1,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))

par(lwd=3)
par(cex.axis=1.2)
plot((xR$E_methCoef-xR$meanMethCoef_ADS),(xR$E_methCoef-xR$meanREMethCoef),pch=16,xlim=c(-1,1),ylim=c(-1,1),col=rgb(0,0,0,0.2),
     xlab="",ylab="",cex=0.8)
abline(h=0,lty=2,lwd=1.5)
abline(v=0,lty=2,lwd=1.5)
par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(-1, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(-1, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()


## ADS_Adipose
########
png(paste(resultsDIR,"figure31_ADS_AdiposeFilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)                                        

def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist((xJ$E_methCoef-xJ$meanMethCoef_ADS_Adipose), breaks=seq(-1,1,0.1), plot=FALSE)
yhist <- hist((xJ$E_methCoef-xJ$meanREMethCoef), breaks=seq(-1,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))

par(lwd=3)
par(cex.axis=1.2)
plot((xJ$E_methCoef-xJ$meanMethCoef_ADS_Adipose),(xJ$E_methCoef-xJ$meanREMethCoef),pch=16,xlim=c(-1,1),ylim=c(-1,1),col=rgb(0,0,0,0.2),
     xlab="",ylab="",cex=0.8)
abline(h=0,lty=2,lwd=1.5)
abline(v=0,lty=2,lwd=1.5)
par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(-1, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(-1, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()


## ADS_IPSC
########
png(paste(resultsDIR,"figure31_ADS_IPSCFilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)                                        

def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist((xL$E_methCoef-xL$meanMethCoef_ADS_IPSC), breaks=seq(-1,1,0.1), plot=FALSE)
yhist <- hist((xL$E_methCoef-xL$meanREMethCoef), breaks=seq(-1,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))

par(lwd=3)
par(cex.axis=1.2)
plot((xL$E_methCoef-xL$meanMethCoef_ADS_IPSC),(xL$E_methCoef-xL$meanREMethCoef),pch=16,xlim=c(-1,1),ylim=c(-1,1),col=rgb(0,0,0,0.2),
     xlab="",ylab="",cex=0.8)
abline(h=0,lty=2,lwd=1.5)
abline(v=0,lty=2,lwd=1.5)
par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(-1, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(-1, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()


####################################
## Without Informativeness filter ##
####################################
selected_reDinucleotide_Element_ADS<-subset(reDinucleotide_Element_ADS,
                                           reDinucleotide_Element_ADS$meanMethCoef_ADS!=-1 & 
                                             reDinucleotide_Element_ADS$E_methCoef!=-1)
selected_reDinucleotide_Element_ADS_Adipose<-subset(reDinucleotide_Element_ADS_Adipose,
                                              reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose!=-1 & 
                                                reDinucleotide_Element_ADS_Adipose$E_methCoef!=-1)
selected_reDinucleotide_Element_ADS_IPSC<-subset(reDinucleotide_Element_ADS_IPSC,
                                                    reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC!=-1 & 
                                                      reDinucleotide_Element_ADS_IPSC$E_methCoef!=-1)



#######################
#library(RColorBrewer)
#cols<-brewer.pal(3,"Dark2")
reducedDataADS<-selected_reDinucleotide_Element_ADS[,c(1,15,7,12)]
reducedDataADS_Adipose<-selected_reDinucleotide_Element_ADS_Adipose[,c(1,15,7,12)]
reducedDataADS_IPSC<-selected_reDinucleotide_Element_ADS_IPSC[,c(1,15,7,12)]



meanValuesDIR<-paste(Sys.getenv("HOME"),"/projects/Evaluation_of_single_CpG_sites_as_proxies_of_CpG_island_methylation_states_at_the_genome_scale/results/20120802/",sep="")

meanValuesADS<-read.table(paste(meanValuesDIR,"meanValueADS_All.txt",sep=""),header=FALSE,sep="\t")
meanValuesADS_Adipose<-read.table(paste(meanValuesDIR,"meanValueADS_Adipose_All.txt",sep=""),header=FALSE,sep="\t")
meanValuesADS_IPSC<-read.table(paste(meanValuesDIR,"meanValueADS_IPSC_All.txt",sep=""),header=FALSE,sep="\t")

colnames(meanValuesADS)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")
colnames(meanValuesADS_Adipose)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")
colnames(meanValuesADS_IPSC)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")

mergedReducedMeanADS<-merge(reducedDataADS,meanValuesADS,by="id_In_Type",all=TRUE)
mergedReducedMeanADS_Adipose<-merge(reducedDataADS_Adipose,meanValuesADS_Adipose,by="id_In_Type",all=TRUE)
mergedReducedMeanADS_IPSC<-merge(reducedDataADS_IPSC,meanValuesADS_IPSC,by="id_In_Type",all=TRUE)


xR<-subset(mergedReducedMeanADS,mergedReducedMeanADS$nREs>1)
xJ<-subset(mergedReducedMeanADS_Adipose,mergedReducedMeanADS_Adipose$nREs>1)
xL<-subset(mergedReducedMeanADS_IPSC,mergedReducedMeanADS_IPSC$nREs>1)

## ADS
#####
png(paste(resultsDIR,"figure31_ADSAllBN.png",sep=""),height=15,width=15,units="cm",res=300)

def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist((xR$E_methCoef-xR$meanMethCoef_ADS), breaks=seq(-1,1,0.1), plot=FALSE)
yhist <- hist((xR$E_methCoef-xR$meanREMethCoef), breaks=seq(-1,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))


par(lwd=3)
par(cex.axis=1.2)
plot((xR$E_methCoef-xR$meanMethCoef_ADS),(xR$E_methCoef-xR$meanREMethCoef),pch=16,xlim=c(-1,1),ylim=c(-1,1),col=rgb(0,0,0,0.2),
     xlab="",ylab="",cex=0.8)
abline(h=0,lty=2,lwd=1.5)
abline(v=0,lty=2,lwd=1.5)
par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(-1, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(-1, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()




## ADS_Adipose
########
png(paste(resultsDIR,"figure31_ADS_AdiposeAllBN.png",sep=""),height=15,width=15,units="cm",res=300)                                        

def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist((xJ$E_methCoef-xJ$meanMethCoef_ADS_Adipose), breaks=seq(-1,1,0.1), plot=FALSE)
yhist <- hist((xJ$E_methCoef-xJ$meanREMethCoef), breaks=seq(-1,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))

par(lwd=3)
par(cex.axis=1.2)
plot((xJ$E_methCoef-xJ$meanMethCoef_ADS_Adipose),(xJ$E_methCoef-xJ$meanREMethCoef),pch=16,xlim=c(-1,1),ylim=c(-1,1),col=rgb(0,0,0,0.2),
     xlab="",ylab="",cex=1)
abline(h=0,lty=2,lwd=1.5)
abline(v=0,lty=2,lwd=1.5)
par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(-1, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(-1, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()

## ADS_IPSC
########
png(paste(resultsDIR,"figure31_ADS_IPSCAllBN.png",sep=""),height=15,width=15,units="cm",res=300)                                        

def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist((xL$E_methCoef-xL$meanMethCoef_ADS_IPSC), breaks=seq(-1,1,0.1), plot=FALSE)
yhist <- hist((xL$E_methCoef-xL$meanREMethCoef), breaks=seq(-1,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))

par(lwd=3)
par(cex.axis=1.2)
plot((xL$E_methCoef-xL$meanMethCoef_ADS_IPSC),(xL$E_methCoef-xL$meanREMethCoef),pch=16,xlim=c(-1,1),ylim=c(-1,1),col=rgb(0,0,0,0.2),
     xlab="",ylab="",cex=1)
abline(h=0,lty=2,lwd=1.5)
abline(v=0,lty=2,lwd=1.5)
par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(-1, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(-1, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()
