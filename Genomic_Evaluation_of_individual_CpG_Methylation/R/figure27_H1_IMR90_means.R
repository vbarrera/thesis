##################
### Figure 27 ###
##################
#################################
## With Informativeness filter ##
#################################

meanValuesDIR<-paste(Sys.getenv("HOME"),"/projects/Evaluation_of_single_CpG_sites_as_proxies_of_CpG_island_methylation_states_at_the_genome_scale/results/20120308/",sep="")

meanValuesH1<-read.table(paste(meanValuesDIR,"meanValueH1.txt",sep=""),header=FALSE,sep="\t")
meanValuesIMR90<-read.table(paste(meanValuesDIR,"meanValueIMR90.txt",sep=""),header=FALSE,sep="\t")

colnames(meanValuesH1)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")
colnames(meanValuesIMR90)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")

#H1


meanValuesH1Jitter<-as.data.frame(cbind(MethCoef=jitter(meanValuesH1$MethCoef,factor=25),meanREMethCoef=jitter(meanValuesH1$meanREMethCoef,factor=25)))
incorrectPoints_H1_plus<-subset(meanValuesH1,meanValuesH1$MethCoef-meanValuesH1$meanREMethCoef>0.25)
#x_H1_plus<-subset(meanValuesH1Jitter,meanValuesH1Jitter$MethCoef-meanValuesH1Jitter$meanREMethCoef>0.25)
dim(incorrectPoints_H1_plus)
x_H1_plus<-subset(meanValuesH1Jitter,meanValuesH1Jitter$MethCoef>=0.75)

incorrectPoints_H1_minus<-subset(meanValuesH1,meanValuesH1$MethCoef-meanValuesH1$meanREMethCoef<(-0.25))
#x_H1_minus<-subset(meanValuesH1Jitter,meanValuesH1Jitter$MethCoef-meanValuesH1Jitter$meanREMethCoef<(-0.25))
dim(incorrectPoints_H1_minus)
x_H1_minus<-subset(meanValuesH1Jitter,meanValuesH1Jitter$MethCoef<=0.25)

correctPoints_H1<-subset(meanValuesH1,abs(meanValuesH1$MethCoef-meanValuesH1$meanREMethCoef)<=0.25)
#y_H1<-subset(meanValuesH1Jitter,abs(meanValuesH1Jitter$MethCoef-meanValuesH1Jitter$meanREMethCoef)<=0.25)
dim(correctPoints_H1)
y_H1<-subset(meanValuesH1Jitter,meanValuesH1Jitter$MethCoef>0.25 & meanValuesH1Jitter$MethCoef<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure27_meansH1FilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure27_meansH1FilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)

xhist <- hist(meanValuesH1$meanREMethCoef, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(meanValuesH1$MethCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (meanValuesH1$MethCoef,meanValuesH1$meanREMethCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n")


points(x_H1_plus$meanREMethCoef,x_H1_plus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_H1_minus$meanREMethCoef,x_H1_minus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_H1$meanREMethCoef,y_H1$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# points(incorrectPoints_H1_plus$meanREMethCoef,incorrectPoints_H1_plus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(incorrectPoints_H1_minus$meanREMethCoef,incorrectPoints_H1_minus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(correctPoints_H1$meanREMethCoef,correctPoints_H1$MethCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if (numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_H1_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_H1_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_H1)[1],cex=1.2,col="red")
}else{}

abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()

# IMR90

meanValuesIMR90Jitter<-as.data.frame(cbind(MethCoef=jitter(meanValuesIMR90$MethCoef,factor=25),meanREMethCoef=jitter(meanValuesIMR90$meanREMethCoef,factor=25)))
incorrectPoints_IMR90_plus<-subset(meanValuesIMR90,meanValuesIMR90$MethCoef-meanValuesIMR90$meanREMethCoef>0.25)
x_IMR90_plus<-subset(meanValuesIMR90Jitter,meanValuesIMR90Jitter$MethCoef-meanValuesIMR90Jitter$meanREMethCoef>0.25)
dim(incorrectPoints_IMR90_plus)
x_IMR90_plus<-subset(meanValuesIMR90Jitter,meanValuesIMR90Jitter$MethCoef>=0.75)


incorrectPoints_IMR90_minus<-subset(meanValuesIMR90,meanValuesIMR90$MethCoef-meanValuesIMR90$meanREMethCoef<(-0.25))
x_IMR90_minus<-subset(meanValuesIMR90Jitter,meanValuesIMR90Jitter$MethCoef-meanValuesIMR90Jitter$meanREMethCoef<(-0.25))
dim(incorrectPoints_IMR90_minus)
x_IMR90_minus<-subset(meanValuesIMR90Jitter,meanValuesIMR90Jitter$MethCoef<=0.25)


correctPoints_IMR90<-subset(meanValuesIMR90,abs(meanValuesIMR90$MethCoef-meanValuesIMR90$meanREMethCoef)<=0.25)
y_IMR90<-subset(meanValuesIMR90Jitter,abs(meanValuesIMR90Jitter$MethCoef-meanValuesIMR90Jitter$meanREMethCoef)<=0.25)
dim(correctPoints_IMR90)
y_IMR90<-subset(meanValuesIMR90Jitter,meanValuesIMR90Jitter$MethCoef>0.25 & meanValuesIMR90Jitter$MethCoef<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure27_meansIMR90FilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure27_meansIMR90FilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}

def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)

xhist <- hist(meanValuesIMR90$meanREMethCoef, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(meanValuesIMR90$MethCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (meanValuesIMR90$MethCoef,meanValuesIMR90$meanREMethCoef,pch=16,col=rgb(0,0,0,0.2),xlab="",
      ylab="",main="",type="n")


points(x_IMR90_plus$meanREMethCoef,x_IMR90_plus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_IMR90_minus$meanREMethCoef,x_IMR90_minus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_IMR90$meanREMethCoef,y_IMR90$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# points(incorrectPoints_IMR90_plus$meanREMethCoef,incorrectPoints_IMR90_plus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(incorrectPoints_IMR90_minus$meanREMethCoef,incorrectPoints_IMR90_minus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(correctPoints_IMR90$meanREMethCoef,correctPoints_IMR90$MethCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if(numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_IMR90_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_IMR90_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_IMR90)[1],cex=1.2,col="red")
}else{}

abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()


####################################
## Without Informativeness filter ##
####################################

meanValuesDIR<-paste(Sys.getenv("HOME"),"/projects/Evaluation_of_single_CpG_sites_as_proxies_of_CpG_island_methylation_states_at_the_genome_scale/results/20120423/",sep="")

meanValuesH1<-read.table(paste(meanValuesDIR,"meanValueH1_v2.txt",sep=""),header=FALSE,sep="\t")
meanValuesIMR90<-read.table(paste(meanValuesDIR,"meanValueIMR90_v2.txt",sep=""),header=FALSE,sep="\t")

colnames(meanValuesH1)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")
colnames(meanValuesIMR90)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")

#H1

meanValuesH1Jitter<-as.data.frame(cbind(MethCoef=jitter(meanValuesH1$MethCoef,factor=25),meanREMethCoef=jitter(meanValuesH1$meanREMethCoef,factor=25)))
incorrectPoints_H1_plus<-subset(meanValuesH1,meanValuesH1$MethCoef-meanValuesH1$meanREMethCoef>0.25)
#x_H1_plus<-subset(meanValuesH1Jitter,meanValuesH1Jitter$MethCoef-meanValuesH1Jitter$meanREMethCoef>0.25)
dim(incorrectPoints_H1_plus)
x_H1_plus<-subset(meanValuesH1Jitter,meanValuesH1Jitter$MethCoef>=0.75)

incorrectPoints_H1_minus<-subset(meanValuesH1,meanValuesH1$MethCoef-meanValuesH1$meanREMethCoef<(-0.25))
#x_H1_minus<-subset(meanValuesH1Jitter,meanValuesH1Jitter$MethCoef-meanValuesH1Jitter$meanREMethCoef<(-0.25))
dim(incorrectPoints_H1_minus)
x_H1_minus<-subset(meanValuesH1Jitter,meanValuesH1Jitter$MethCoef<=0.25)

correctPoints_H1<-subset(meanValuesH1,abs(meanValuesH1$MethCoef-meanValuesH1$meanREMethCoef)<=0.25)
#y_H1<-subset(meanValuesH1Jitter,abs(meanValuesH1Jitter$MethCoef-meanValuesH1Jitter$meanREMethCoef)<=0.25)
dim(correctPoints_H1)
y_H1<-subset(meanValuesH1Jitter,meanValuesH1Jitter$MethCoef>0.25 & meanValuesH1Jitter$MethCoef<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure27_meansH1AllBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure27_meansH1AllBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)

xhist <- hist(meanValuesH1$meanREMethCoef, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(meanValuesH1$MethCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (meanValuesH1$MethCoef,meanValuesH1$meanREMethCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n")


points(x_H1_plus$meanREMethCoef,x_H1_plus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=1)
points(x_H1_minus$meanREMethCoef,x_H1_minus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=1)
points(y_H1$meanREMethCoef,y_H1$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=1)

# points(incorrectPoints_H1_plus$meanREMethCoef,incorrectPoints_H1_plus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(incorrectPoints_H1_minus$meanREMethCoef,incorrectPoints_H1_minus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(correctPoints_H1$meanREMethCoef,correctPoints_H1$MethCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if (numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_H1_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_H1_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_H1)[1],cex=1.2,col="red")
}else{}

abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()

# IMR90

meanValuesIMR90Jitter<-as.data.frame(cbind(MethCoef=jitter(meanValuesIMR90$MethCoef,factor=25),meanREMethCoef=jitter(meanValuesIMR90$meanREMethCoef,factor=25)))
incorrectPoints_IMR90_plus<-subset(meanValuesIMR90,meanValuesIMR90$MethCoef-meanValuesIMR90$meanREMethCoef>0.25)
#x_IMR90_plus<-subset(meanValuesIMR90Jitter,meanValuesIMR90Jitter$MethCoef-meanValuesIMR90Jitter$meanREMethCoef>0.25)
dim(incorrectPoints_IMR90_plus)
x_IMR90_plus<-subset(meanValuesIMR90Jitter,meanValuesIMR90Jitter$MethCoef>=0.75)

incorrectPoints_IMR90_minus<-subset(meanValuesIMR90,meanValuesIMR90$MethCoef-meanValuesIMR90$meanREMethCoef<(-0.25))
#x_IMR90_minus<-subset(meanValuesIMR90Jitter,meanValuesIMR90Jitter$MethCoef-meanValuesIMR90Jitter$meanREMethCoef<(-0.25))
dim(incorrectPoints_IMR90_minus)
x_IMR90_minus<-subset(meanValuesIMR90Jitter,meanValuesIMR90Jitter$MethCoef<=0.25)

correctPoints_IMR90<-subset(meanValuesIMR90,abs(meanValuesIMR90$MethCoef-meanValuesIMR90$meanREMethCoef)<=0.25)
#y_IMR90<-subset(meanValuesIMR90Jitter,abs(meanValuesIMR90Jitter$MethCoef-meanValuesIMR90Jitter$meanREMethCoef)<=0.25)
dim(correctPoints_IMR90)
y_IMR90<-subset(meanValuesIMR90Jitter,meanValuesIMR90Jitter$MethCoef>0.25 & meanValuesIMR90Jitter$MethCoef<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure27_meansIMR90AllBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure27_meansIMR90AllBN.png",sep=""),height=15,width=15,units="cm",res=300)
}

def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)

xhist <- hist(meanValuesIMR90$meanREMethCoef, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(meanValuesIMR90$MethCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (meanValuesIMR90$MethCoef,meanValuesIMR90$meanREMethCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n")


points(x_IMR90_plus$meanREMethCoef,x_IMR90_plus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=1)
points(x_IMR90_minus$meanREMethCoef,x_IMR90_minus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=1)
points(y_IMR90$meanREMethCoef,y_IMR90$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=1)

# points(incorrectPoints_IMR90_plus$meanREMethCoef,incorrectPoints_IMR90_plus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(incorrectPoints_IMR90_minus$meanREMethCoef,incorrectPoints_IMR90_minus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(correctPoints_IMR90$meanREMethCoef,correctPoints_IMR90$MethCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)


if(numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_IMR90_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_IMR90_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_IMR90)[1],cex=1.2,col="red")
}else{}

abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")


par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()