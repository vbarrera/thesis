##################
### Figure 27 ###
##################
#################################
## With Informativeness filter ##
#################################

meanValuesDIR<-paste(Sys.getenv("HOME"),"/projects/Evaluation_of_single_CpG_sites_as_proxies_of_CpG_island_methylation_states_at_the_genome_scale/results/20120802/",sep="")

meanValuesADS<-read.table(paste(meanValuesDIR,"meanValueADS_Filtered.txt",sep=""),header=FALSE,sep="\t")
meanValuesADS_Adipose<-read.table(paste(meanValuesDIR,"meanValueADS_Adipose_Filtered.txt",sep=""),header=FALSE,sep="\t")
meanValuesADS_IPSC<-read.table(paste(meanValuesDIR,"meanValueADS_IPSC_Filtered.txt",sep=""),header=FALSE,sep="\t")


colnames(meanValuesADS)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")
colnames(meanValuesADS_Adipose)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")
colnames(meanValuesADS_IPSC)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")

#ADS
# lmADSFiltered<-lm(meanValuesADS$MethCoef~meanValuesADS$meanREMethCoef)
# summary(lmADSFiltered)

meanValuesADSJitter<-as.data.frame(cbind(MethCoef=jitter(meanValuesADS$MethCoef,factor=25),meanREMethCoef=jitter(meanValuesADS$meanREMethCoef,factor=25)))
incorrectPoints_ADS_plus<-subset(meanValuesADS,meanValuesADS$MethCoef-meanValuesADS$meanREMethCoef>0.25)
#x_ADS_plus<-subset(meanValuesADSJitter,meanValuesADSJitter$MethCoef-meanValuesADSJitter$meanREMethCoef>0.25)
dim(incorrectPoints_ADS_plus)
x_ADS_plus<-subset(meanValuesADSJitter,meanValuesADSJitter$MethCoef>=0.75)

incorrectPoints_ADS_minus<-subset(meanValuesADS,meanValuesADS$MethCoef-meanValuesADS$meanREMethCoef<(-0.25))
#x_ADS_minus<-subset(meanValuesADSJitter,meanValuesADSJitter$MethCoef-meanValuesADSJitter$meanREMethCoef<(-0.25))
dim(incorrectPoints_ADS_minus)
x_ADS_minus<-subset(meanValuesADSJitter,meanValuesADSJitter$MethCoef<=0.25)

correctPoints_ADS<-subset(meanValuesADS,abs(meanValuesADS$MethCoef-meanValuesADS$meanREMethCoef)<=0.25)
#y_ADS<-subset(meanValuesADSJitter,abs(meanValuesADSJitter$MethCoef-meanValuesADSJitter$meanREMethCoef)<=0.25)
dim(correctPoints_ADS)
y_ADS<-subset(meanValuesADSJitter,meanValuesADSJitter$MethCoef>0.25 & meanValuesADSJitter$MethCoef<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure27_meansADSFilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure27_meansADSFilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)

xhist <- hist(meanValuesADS$meanREMethCoef, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(meanValuesADS$MethCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (meanValuesADS$MethCoef,meanValuesADS$meanREMethCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n")


points(x_ADS_plus$meanREMethCoef,x_ADS_plus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_ADS_minus$meanREMethCoef,x_ADS_minus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_ADS$meanREMethCoef,y_ADS$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# points(incorrectPoints_ADS_plus$meanREMethCoef,incorrectPoints_ADS_plus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(incorrectPoints_ADS_minus$meanREMethCoef,incorrectPoints_ADS_minus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(correctPoints_ADS$meanREMethCoef,correctPoints_ADS$MethCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if (numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_ADS_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_ADS_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_ADS)[1],cex=1.2,col="red")
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

# ADS_Adipose

# lmADS_AdiposeFiltered<-lm(meanValuesADS_Adipose$MethCoef~meanValuesADS_Adipose$meanREMethCoef)
# summary(lmADS_AdiposeFiltered)

meanValuesADS_AdiposeJitter<-as.data.frame(cbind(MethCoef=jitter(meanValuesADS_Adipose$MethCoef,factor=25),meanREMethCoef=jitter(meanValuesADS_Adipose$meanREMethCoef,factor=25)))
incorrectPoints_ADS_Adipose_plus<-subset(meanValuesADS_Adipose,meanValuesADS_Adipose$MethCoef-meanValuesADS_Adipose$meanREMethCoef>0.25)
x_ADS_Adipose_plus<-subset(meanValuesADS_AdiposeJitter,meanValuesADS_AdiposeJitter$MethCoef-meanValuesADS_AdiposeJitter$meanREMethCoef>0.25)
dim(incorrectPoints_ADS_Adipose_plus)
x_ADS_Adipose_plus<-subset(meanValuesADS_AdiposeJitter,meanValuesADS_AdiposeJitter$MethCoef>=0.75)


incorrectPoints_ADS_Adipose_minus<-subset(meanValuesADS_Adipose,meanValuesADS_Adipose$MethCoef-meanValuesADS_Adipose$meanREMethCoef<(-0.25))
x_ADS_Adipose_minus<-subset(meanValuesADS_AdiposeJitter,meanValuesADS_AdiposeJitter$MethCoef-meanValuesADS_AdiposeJitter$meanREMethCoef<(-0.25))
dim(incorrectPoints_ADS_Adipose_minus)
x_ADS_Adipose_minus<-subset(meanValuesADS_AdiposeJitter,meanValuesADS_AdiposeJitter$MethCoef<=0.25)


correctPoints_ADS_Adipose<-subset(meanValuesADS_Adipose,abs(meanValuesADS_Adipose$MethCoef-meanValuesADS_Adipose$meanREMethCoef)<=0.25)
y_ADS_Adipose<-subset(meanValuesADS_AdiposeJitter,abs(meanValuesADS_AdiposeJitter$MethCoef-meanValuesADS_AdiposeJitter$meanREMethCoef)<=0.25)
dim(correctPoints_ADS_Adipose)
y_ADS_Adipose<-subset(meanValuesADS_AdiposeJitter,meanValuesADS_AdiposeJitter$MethCoef>0.25 & meanValuesADS_AdiposeJitter$MethCoef<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure27_meansADS_AdiposeFilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure27_meansADS_AdiposeFilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}

def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)

xhist <- hist(meanValuesADS_Adipose$meanREMethCoef, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(meanValuesADS_Adipose$MethCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (meanValuesADS_Adipose$MethCoef,meanValuesADS_Adipose$meanREMethCoef,pch=16,col=rgb(0,0,0,0.2),xlab="",
      ylab="",main="",type="n")


points(x_ADS_Adipose_plus$meanREMethCoef,x_ADS_Adipose_plus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_ADS_Adipose_minus$meanREMethCoef,x_ADS_Adipose_minus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_ADS_Adipose$meanREMethCoef,y_ADS_Adipose$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# points(incorrectPoints_ADS_Adipose_plus$meanREMethCoef,incorrectPoints_ADS_Adipose_plus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(incorrectPoints_ADS_Adipose_minus$meanREMethCoef,incorrectPoints_ADS_Adipose_minus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(correctPoints_ADS_Adipose$meanREMethCoef,correctPoints_ADS_Adipose$MethCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if(numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_ADS_Adipose_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_ADS_Adipose_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_ADS_Adipose)[1],cex=1.2,col="red")
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

# ADS_IPSC

# lmADS_IPSCFiltered<-lm(meanValuesADS_IPSC$MethCoef~meanValuesADS_IPSC$meanREMethCoef)
# summary(lmADS_IPSCFiltered)

meanValuesADS_IPSCJitter<-as.data.frame(cbind(MethCoef=jitter(meanValuesADS_IPSC$MethCoef,factor=25),meanREMethCoef=jitter(meanValuesADS_IPSC$meanREMethCoef,factor=25)))
incorrectPoints_ADS_IPSC_plus<-subset(meanValuesADS_IPSC,meanValuesADS_IPSC$MethCoef-meanValuesADS_IPSC$meanREMethCoef>0.25)
x_ADS_IPSC_plus<-subset(meanValuesADS_IPSCJitter,meanValuesADS_IPSCJitter$MethCoef-meanValuesADS_IPSCJitter$meanREMethCoef>0.25)
dim(incorrectPoints_ADS_IPSC_plus)
x_ADS_IPSC_plus<-subset(meanValuesADS_IPSCJitter,meanValuesADS_IPSCJitter$MethCoef>=0.75)


incorrectPoints_ADS_IPSC_minus<-subset(meanValuesADS_IPSC,meanValuesADS_IPSC$MethCoef-meanValuesADS_IPSC$meanREMethCoef<(-0.25))
x_ADS_IPSC_minus<-subset(meanValuesADS_IPSCJitter,meanValuesADS_IPSCJitter$MethCoef-meanValuesADS_IPSCJitter$meanREMethCoef<(-0.25))
dim(incorrectPoints_ADS_IPSC_minus)
x_ADS_IPSC_minus<-subset(meanValuesADS_IPSCJitter,meanValuesADS_IPSCJitter$MethCoef<=0.25)


correctPoints_ADS_IPSC<-subset(meanValuesADS_IPSC,abs(meanValuesADS_IPSC$MethCoef-meanValuesADS_IPSC$meanREMethCoef)<=0.25)
y_ADS_IPSC<-subset(meanValuesADS_IPSCJitter,abs(meanValuesADS_IPSCJitter$MethCoef-meanValuesADS_IPSCJitter$meanREMethCoef)<=0.25)
dim(correctPoints_ADS_IPSC)
y_ADS_IPSC<-subset(meanValuesADS_IPSCJitter,meanValuesADS_IPSCJitter$MethCoef>0.25 & meanValuesADS_IPSCJitter$MethCoef<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure27_meansADS_IPSCFilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
  png(paste(resultsDIR,"figure27_meansADS_IPSCFilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}

def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)

xhist <- hist(meanValuesADS_IPSC$meanREMethCoef, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(meanValuesADS_IPSC$MethCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (meanValuesADS_IPSC$MethCoef,meanValuesADS_IPSC$meanREMethCoef,pch=16,col=rgb(0,0,0,0.2),xlab="",
      ylab="",main="",type="n")


points(x_ADS_IPSC_plus$meanREMethCoef,x_ADS_IPSC_plus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_ADS_IPSC_minus$meanREMethCoef,x_ADS_IPSC_minus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_ADS_IPSC$meanREMethCoef,y_ADS_IPSC$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# points(incorrectPoints_ADS_IPSC_plus$meanREMethCoef,incorrectPoints_ADS_IPSC_plus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(incorrectPoints_ADS_IPSC_minus$meanREMethCoef,incorrectPoints_ADS_IPSC_minus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(correctPoints_ADS_IPSC$meanREMethCoef,correctPoints_ADS_IPSC$MethCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if(numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_ADS_IPSC_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_ADS_IPSC_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_ADS_IPSC)[1],cex=1.2,col="red")
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

meanValuesDIR<-paste(Sys.getenv("HOME"),"/projects/Evaluation_of_single_CpG_sites_as_proxies_of_CpG_island_methylation_states_at_the_genome_scale/results/20120802/",sep="")

meanValuesADS<-read.table(paste(meanValuesDIR,"meanValueADS_All.txt",sep=""),header=FALSE,sep="\t")
meanValuesADS_Adipose<-read.table(paste(meanValuesDIR,"meanValueADS_Adipose_All.txt",sep=""),header=FALSE,sep="\t")
meanValuesADS_IPSC<-read.table(paste(meanValuesDIR,"meanValueADS_IPSC_All.txt",sep=""),header=FALSE,sep="\t")


colnames(meanValuesADS)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")
colnames(meanValuesADS_Adipose)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")
colnames(meanValuesADS_IPSC)<-c("id_In_Type","MethCoef","nREs","meanREMethCoef")

#ADS
# lmADSAll<-lm(meanValuesADS$MethCoef~meanValuesADS$meanREMethCoef)
# summary(lmADSAll)

meanValuesADSJitter<-as.data.frame(cbind(MethCoef=jitter(meanValuesADS$MethCoef,factor=25),meanREMethCoef=jitter(meanValuesADS$meanREMethCoef,factor=25)))
incorrectPoints_ADS_plus<-subset(meanValuesADS,meanValuesADS$MethCoef-meanValuesADS$meanREMethCoef>0.25)
#x_ADS_plus<-subset(meanValuesADSJitter,meanValuesADSJitter$MethCoef-meanValuesADSJitter$meanREMethCoef>0.25)
dim(incorrectPoints_ADS_plus)
x_ADS_plus<-subset(meanValuesADSJitter,meanValuesADSJitter$MethCoef>=0.75)

incorrectPoints_ADS_minus<-subset(meanValuesADS,meanValuesADS$MethCoef-meanValuesADS$meanREMethCoef<(-0.25))
#x_ADS_minus<-subset(meanValuesADSJitter,meanValuesADSJitter$MethCoef-meanValuesADSJitter$meanREMethCoef<(-0.25))
dim(incorrectPoints_ADS_minus)
x_ADS_minus<-subset(meanValuesADSJitter,meanValuesADSJitter$MethCoef<=0.25)

correctPoints_ADS<-subset(meanValuesADS,abs(meanValuesADS$MethCoef-meanValuesADS$meanREMethCoef)<=0.25)
#y_ADS<-subset(meanValuesADSJitter,abs(meanValuesADSJitter$MethCoef-meanValuesADSJitter$meanREMethCoef)<=0.25)
dim(correctPoints_ADS)
y_ADS<-subset(meanValuesADSJitter,meanValuesADSJitter$MethCoef>0.25 & meanValuesADSJitter$MethCoef<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure27_meansADSAllBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure27_meansADSAllBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)

xhist <- hist(meanValuesADS$meanREMethCoef, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(meanValuesADS$MethCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (meanValuesADS$MethCoef,meanValuesADS$meanREMethCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n")


points(x_ADS_plus$meanREMethCoef,x_ADS_plus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=1)
points(x_ADS_minus$meanREMethCoef,x_ADS_minus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=1)
points(y_ADS$meanREMethCoef,y_ADS$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=1)

# points(incorrectPoints_ADS_plus$meanREMethCoef,incorrectPoints_ADS_plus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(incorrectPoints_ADS_minus$meanREMethCoef,incorrectPoints_ADS_minus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(correctPoints_ADS$meanREMethCoef,correctPoints_ADS$MethCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if (numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_ADS_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_ADS_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_ADS)[1],cex=1.2,col="red")
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

# ADS_Adipose
# lmADS_AdiposeAll<-lm(meanValuesADS_Adipose$MethCoef~meanValuesADS_Adipose$meanREMethCoef)
# summary(lmADS_AdiposeAll)

meanValuesADS_AdiposeJitter<-as.data.frame(cbind(MethCoef=jitter(meanValuesADS_Adipose$MethCoef,factor=25),meanREMethCoef=jitter(meanValuesADS_Adipose$meanREMethCoef,factor=25)))
incorrectPoints_ADS_Adipose_plus<-subset(meanValuesADS_Adipose,meanValuesADS_Adipose$MethCoef-meanValuesADS_Adipose$meanREMethCoef>0.25)
#x_ADS_Adipose_plus<-subset(meanValuesADS_AdiposeJitter,meanValuesADS_AdiposeJitter$MethCoef-meanValuesADS_AdiposeJitter$meanREMethCoef>0.25)
dim(incorrectPoints_ADS_Adipose_plus)
x_ADS_Adipose_plus<-subset(meanValuesADS_AdiposeJitter,meanValuesADS_AdiposeJitter$MethCoef>=0.75)

incorrectPoints_ADS_Adipose_minus<-subset(meanValuesADS_Adipose,meanValuesADS_Adipose$MethCoef-meanValuesADS_Adipose$meanREMethCoef<(-0.25))
#x_ADS_Adipose_minus<-subset(meanValuesADS_AdiposeJitter,meanValuesADS_AdiposeJitter$MethCoef-meanValuesADS_AdiposeJitter$meanREMethCoef<(-0.25))
dim(incorrectPoints_ADS_Adipose_minus)
x_ADS_Adipose_minus<-subset(meanValuesADS_AdiposeJitter,meanValuesADS_AdiposeJitter$MethCoef<=0.25)

correctPoints_ADS_Adipose<-subset(meanValuesADS_Adipose,abs(meanValuesADS_Adipose$MethCoef-meanValuesADS_Adipose$meanREMethCoef)<=0.25)
#y_ADS_Adipose<-subset(meanValuesADS_AdiposeJitter,abs(meanValuesADS_AdiposeJitter$MethCoef-meanValuesADS_AdiposeJitter$meanREMethCoef)<=0.25)
dim(correctPoints_ADS_Adipose)
y_ADS_Adipose<-subset(meanValuesADS_AdiposeJitter,meanValuesADS_AdiposeJitter$MethCoef>0.25 & meanValuesADS_AdiposeJitter$MethCoef<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure27_meansADS_AdiposeAllBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure27_meansADS_AdiposeAllBN.png",sep=""),height=15,width=15,units="cm",res=300)
}

def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)

xhist <- hist(meanValuesADS_Adipose$meanREMethCoef, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(meanValuesADS_Adipose$MethCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (meanValuesADS_Adipose$MethCoef,meanValuesADS_Adipose$meanREMethCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n")


points(x_ADS_Adipose_plus$meanREMethCoef,x_ADS_Adipose_plus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=1)
points(x_ADS_Adipose_minus$meanREMethCoef,x_ADS_Adipose_minus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=1)
points(y_ADS_Adipose$meanREMethCoef,y_ADS_Adipose$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=1)

# points(incorrectPoints_ADS_Adipose_plus$meanREMethCoef,incorrectPoints_ADS_Adipose_plus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(incorrectPoints_ADS_Adipose_minus$meanREMethCoef,incorrectPoints_ADS_Adipose_minus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(correctPoints_ADS_Adipose$meanREMethCoef,correctPoints_ADS_Adipose$MethCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)


if(numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_ADS_Adipose_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_ADS_Adipose_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_ADS_Adipose)[1],cex=1.2,col="red")
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

# ADS_Adipose

# lmADS_IPSCAll<-lm(meanValuesADS_IPSC$MethCoef~meanValuesADS_IPSC$meanREMethCoef)
# summary(lmADS_IPSCAll)

meanValuesADS_IPSCJitter<-as.data.frame(cbind(MethCoef=jitter(meanValuesADS_IPSC$MethCoef,factor=25),meanREMethCoef=jitter(meanValuesADS_IPSC$meanREMethCoef,factor=25)))
incorrectPoints_ADS_IPSC_plus<-subset(meanValuesADS_IPSC,meanValuesADS_IPSC$MethCoef-meanValuesADS_IPSC$meanREMethCoef>0.25)
#x_ADS_IPSC_plus<-subset(meanValuesADS_IPSCJitter,meanValuesADS_IPSCJitter$MethCoef-meanValuesADS_IPSCJitter$meanREMethCoef>0.25)
dim(incorrectPoints_ADS_IPSC_plus)
x_ADS_IPSC_plus<-subset(meanValuesADS_IPSCJitter,meanValuesADS_IPSCJitter$MethCoef>=0.75)

incorrectPoints_ADS_IPSC_minus<-subset(meanValuesADS_IPSC,meanValuesADS_IPSC$MethCoef-meanValuesADS_IPSC$meanREMethCoef<(-0.25))
#x_ADS_IPSC_minus<-subset(meanValuesADS_IPSCJitter,meanValuesADS_IPSCJitter$MethCoef-meanValuesADS_IPSCJitter$meanREMethCoef<(-0.25))
dim(incorrectPoints_ADS_IPSC_minus)
x_ADS_IPSC_minus<-subset(meanValuesADS_IPSCJitter,meanValuesADS_IPSCJitter$MethCoef<=0.25)

correctPoints_ADS_IPSC<-subset(meanValuesADS_IPSC,abs(meanValuesADS_IPSC$MethCoef-meanValuesADS_IPSC$meanREMethCoef)<=0.25)
#y_ADS_IPSC<-subset(meanValuesADS_IPSCJitter,abs(meanValuesADS_IPSCJitter$MethCoef-meanValuesADS_IPSCJitter$meanREMethCoef)<=0.25)
dim(correctPoints_ADS_IPSC)
y_ADS_IPSC<-subset(meanValuesADS_IPSCJitter,meanValuesADS_IPSCJitter$MethCoef>0.25 & meanValuesADS_IPSCJitter$MethCoef<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure27_meansADS_IPSCAllBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
  png(paste(resultsDIR,"figure27_meansADS_IPSCAllBN.png",sep=""),height=15,width=15,units="cm",res=300)
}

def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)

xhist <- hist(meanValuesADS_IPSC$meanREMethCoef, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(meanValuesADS_IPSC$MethCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (meanValuesADS_IPSC$MethCoef,meanValuesADS_IPSC$meanREMethCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n")


points(x_ADS_IPSC_plus$meanREMethCoef,x_ADS_IPSC_plus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=1)
points(x_ADS_IPSC_minus$meanREMethCoef,x_ADS_IPSC_minus$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=1)
points(y_ADS_IPSC$meanREMethCoef,y_ADS_IPSC$MethCoef,pch=16,col=rgb(0,0,0,0.2),cex=1)

# points(incorrectPoints_ADS_IPSC_plus$meanREMethCoef,incorrectPoints_ADS_IPSC_plus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(incorrectPoints_ADS_IPSC_minus$meanREMethCoef,incorrectPoints_ADS_IPSC_minus$MethCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
# points(correctPoints_ADS_IPSC$meanREMethCoef,correctPoints_ADS_IPSC$MethCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)


if(numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_ADS_IPSC_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_ADS_IPSC_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_ADS_IPSC)[1],cex=1.2,col="red")
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