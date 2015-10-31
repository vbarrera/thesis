################
### Figure 31 ###
################

#################################
## With Informativeness filter ##
#################################
selected_reDinucleotide_Element_H1<-subset(reDinucleotide_Element_H1,
                                           reDinucleotide_Element_H1$meanMethCoef_H1!=-1 & 
                                             reDinucleotide_Element_H1$E_methCoef!=-1 & ((reDinucleotide_Element_H1$E_posInf/2)/reDinucleotide_Element_H1$E_nCG)>=0.25)
selected_reDinucleotide_Element_IMR90<-subset(reDinucleotide_Element_IMR90,
                                              reDinucleotide_Element_IMR90$meanMethCoef_IMR90!=-1 & 
                                                reDinucleotide_Element_IMR90$E_methCoef!=-1 & ((reDinucleotide_Element_IMR90$E_posInf/2)/reDinucleotide_Element_IMR90$E_nCG)>=0.25)
##########################

#library(RColorBrewer)
#cols<-brewer.pal(3,"Dark2")
reducedDataH1<-selected_reDinucleotide_Element_H1[,c(1,15,7,12)]
reducedDataIMR90<-selected_reDinucleotide_Element_IMR90[,c(1,15,7,12)]

meanValuesDIR<-paste(Sys.getenv("HOME"),"/projects/Evaluation_of_single_CpG_sites_as_proxies_of_CpG_island_methylation_states_at_the_genome_scale/results/20120308/",sep="")

meanValuesH1<-read.table(paste(meanValuesDIR,"meanValueH1.txt",sep=""),header=FALSE,sep="\t")
meanValuesIMR90<-read.table(paste(meanValuesDIR,"meanValueIMR90.txt",sep=""),header=FALSE,sep="\t")

colnames(meanValuesH1)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")
colnames(meanValuesIMR90)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")

mergedReducedMeanH1<-merge(reducedDataH1,meanValuesH1,by="id_In_Type",all=TRUE)
mergedReducedMeanIMR90<-merge(reducedDataIMR90,meanValuesIMR90,by="id_In_Type",all=TRUE)

xR<-subset(mergedReducedMeanH1,mergedReducedMeanH1$nREs>1)
xJ<-subset(mergedReducedMeanIMR90,mergedReducedMeanIMR90$nREs>1)

### H1
######
png(paste(resultsDIR,"figure31_H1FilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist((xR$E_methCoef-xR$meanMethCoef_H1), breaks=seq(-1,1,0.1), plot=FALSE)
yhist <- hist((xR$E_methCoef-xR$meanREMethCoef), breaks=seq(-1,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))

par(lwd=3)
par(cex.axis=1.2)
plot((xR$E_methCoef-xR$meanMethCoef_H1),(xR$E_methCoef-xR$meanREMethCoef),pch=16,xlim=c(-1,1),ylim=c(-1,1),col=rgb(0,0,0,0.2),
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


## IMR90
########
png(paste(resultsDIR,"figure31_IMR90FilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)                                        

def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist((xJ$E_methCoef-xJ$meanMethCoef_IMR90), breaks=seq(-1,1,0.1), plot=FALSE)
yhist <- hist((xJ$E_methCoef-xJ$meanREMethCoef), breaks=seq(-1,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))

par(lwd=3)
par(cex.axis=1.2)
plot((xJ$E_methCoef-xJ$meanMethCoef_IMR90),(xJ$E_methCoef-xJ$meanREMethCoef),pch=16,xlim=c(-1,1),ylim=c(-1,1),col=rgb(0,0,0,0.2),
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
selected_reDinucleotide_Element_H1<-subset(reDinucleotide_Element_H1,
                                           reDinucleotide_Element_H1$meanMethCoef_H1!=-1 & 
                                             reDinucleotide_Element_H1$E_methCoef!=-1)
selected_reDinucleotide_Element_IMR90<-subset(reDinucleotide_Element_IMR90,
                                              reDinucleotide_Element_IMR90$meanMethCoef_IMR90!=-1 & 
                                                reDinucleotide_Element_IMR90$E_methCoef!=-1)

#######################
#library(RColorBrewer)
#cols<-brewer.pal(3,"Dark2")
reducedDataH1<-selected_reDinucleotide_Element_H1[,c(1,15,7,12)]
reducedDataIMR90<-selected_reDinucleotide_Element_IMR90[,c(1,15,7,12)]

meanValuesDIR<-paste(Sys.getenv("HOME"),"/projects/Evaluation_of_single_CpG_sites_as_proxies_of_CpG_island_methylation_states_at_the_genome_scale/results/20120423/",sep="")

meanValuesH1<-read.table(paste(meanValuesDIR,"meanValueH1_v2.txt",sep=""),header=FALSE,sep="\t")
meanValuesIMR90<-read.table(paste(meanValuesDIR,"meanValueIMR90_v2.txt",sep=""),header=FALSE,sep="\t")

colnames(meanValuesH1)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")
colnames(meanValuesIMR90)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")

mergedReducedMeanH1<-merge(reducedDataH1,meanValuesH1,by="id_In_Type",all=TRUE)
mergedReducedMeanIMR90<-merge(reducedDataIMR90,meanValuesIMR90,by="id_In_Type",all=TRUE)

xR<-subset(mergedReducedMeanH1,mergedReducedMeanH1$nREs>1)
xJ<-subset(mergedReducedMeanIMR90,mergedReducedMeanIMR90$nREs>1)

## H1
#####
png(paste(resultsDIR,"figure31_H1AllBN.png",sep=""),height=15,width=15,units="cm",res=300)

def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist((xR$E_methCoef-xR$meanMethCoef_H1), breaks=seq(-1,1,0.1), plot=FALSE)
yhist <- hist((xR$E_methCoef-xR$meanREMethCoef), breaks=seq(-1,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))


par(lwd=3)
par(cex.axis=1.2)
plot((xR$E_methCoef-xR$meanMethCoef_H1),(xR$E_methCoef-xR$meanREMethCoef),pch=16,xlim=c(-1,1),ylim=c(-1,1),col=rgb(0,0,0,0.2),
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




## IMR90
########
png(paste(resultsDIR,"figure31_IMR90AllBN.png",sep=""),height=15,width=15,units="cm",res=300)                                        

def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist((xJ$E_methCoef-xJ$meanMethCoef_IMR90), breaks=seq(-1,1,0.1), plot=FALSE)
yhist <- hist((xJ$E_methCoef-xJ$meanREMethCoef), breaks=seq(-1,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))

par(lwd=3)
par(cex.axis=1.2)
plot((xJ$E_methCoef-xJ$meanMethCoef_IMR90),(xJ$E_methCoef-xJ$meanREMethCoef),pch=16,xlim=c(-1,1),ylim=c(-1,1),col=rgb(0,0,0,0.2),
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
