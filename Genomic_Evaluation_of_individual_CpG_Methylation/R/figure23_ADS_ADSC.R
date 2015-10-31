################
### Figure 23 ###
################

ADSCgiFolder<-paste(baseDIR,"results/20120802/",sep="")
ADS_AdiposeCgiFolder<-paste(baseDIR,"results/20120802/",sep="")


datos_ADS<-read.table(paste(ADSCgiFolder,"hg18_reads_ads_CGi_status_sd.txt",sep=""))
datos_ADS_Adipose<-read.table(paste(ADS_AdiposeCgiFolder,"hg18_reads_ads_adipose_CGi_status_sd.txt",sep=""))
colnames(datos_ADS)<-c("id","nCGt","chr","chrStart","chrEnd","CGinf","methCoef","SD","nReads")
colnames(datos_ADS_Adipose)<-c("id","nCGt","chr","chrStart","chrEnd","CGinf","methCoef","SD","nReads")

#################################
## With Informativeness filter ##
#################################
validDatosADS<-subset(datos_ADS,datos_ADS$methCoef != -1.00 & 
  (datos_ADS$CGinf/2)/datos_ADS$nCGt >=0.25)
validDatosADS_Adipose<-subset(datos_ADS_Adipose,datos_ADS_Adipose$methCoef != -1.00 & 
  (datos_ADS_Adipose$CGinf/2)/datos_ADS_Adipose$nCGt >=0.25)
validDatosADS_unmeth<-subset(validDatosADS, 
                            validDatosADS$methCoef<=0.25)
validDatosADS_halfmeth<-subset(validDatosADS, 
                              validDatosADS$methCoef>0.25 & validDatosADS$methCoef<0.75)
validDatosADS_meth<-subset(validDatosADS, 
                          validDatosADS$methCoef>=0.75)
validDatosADS_Adipose_unmeth<-subset(validDatosADS_Adipose, 
                               validDatosADS_Adipose$methCoef<=0.25)
validDatosADS_Adipose_halfmeth<-subset(validDatosADS_Adipose, 
                                 validDatosADS_Adipose$methCoef>0.25 & validDatosADS_Adipose$methCoef<0.75)
validDatosADS_Adipose_meth<-subset(validDatosADS_Adipose, 
                             validDatosADS_Adipose$methCoef>=0.75)

#Plotting

#ADS

if (numbered=='TRUE'){
  png(paste(resultsDIR,"figure23ADSFilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure23ADSFilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)

par(lwd=3)
par(cex.axis=1.2)

xhist <- hist(validDatosADS$methCoef, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(validDatosADS$SD, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)


par(mar=c(5,5,1,1))
plot(validDatosADS$methCoef,validDatosADS$SD,
     xlab="",ylab="",pch=20,xlim=c(0,1),col=rgb(0,0,0,0.2),type="n",ylim=c(0,0.6))

## Plot the points with some noise (jitter)
points(jitter(validDatosADS_unmeth$methCoef),jitter(validDatosADS_unmeth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(jitter(validDatosADS_halfmeth$methCoef),jitter(validDatosADS_halfmeth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(jitter(validDatosADS_meth$methCoef),jitter(validDatosADS_meth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# Plot the points without noise
# points(validDatosADS_unmeth$methCoef,validDatosADS_unmeth$SD,pch=16,col=rgb(0,1,0,0.2),cex=1)
# points(validDatosADS_halfmeth$methCoef,validDatosADS_halfmeth$SD,pch=16,col=rgb(0,0,1,0.2),cex=1)
# points(validDatosADS_meth$methCoef,validDatosADS_meth$SD,pch=16,col=rgb(1,0,0,0.2),cex=1)


if (numbered=='TRUE'){
# Adding the number of points
text(0.1,0.4,dim(validDatosADS_unmeth)[1],cex=1.2,col="red")
text(0.5,0.4,dim(validDatosADS_halfmeth)[1],cex=1.2,col="red")
text(0.9,0.4,dim(validDatosADS_meth)[1],cex=1.2,col="red")
}else{}
abline(v=0.25,lty=2,lwd=1.5)
abline(v=0.75,lty=2,lwd=1.5)
par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        )
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)

dev.off()


#ADS_Adipose 

if (numbered=='TRUE'){
png(paste(resultsDIR,"figure23ADS_AdiposeFilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
  png(paste(resultsDIR,"figure23ADS_AdiposeFilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)

par(lwd=3)
par(cex.axis=1.2)

xhist <- hist(validDatosADS_Adipose$methCoef, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(validDatosADS_Adipose$SD, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot(validDatosADS_Adipose$methCoef,validDatosADS_Adipose$SD,
     xlab="", ylab="",pch=20,xlim=c(0,1),col=rgb(0,0,0,0.2),main="",type="n",ylim=c(0,0.6))

## Plot the points with some noise (jitter)
points(jitter(validDatosADS_Adipose_unmeth$methCoef),jitter(validDatosADS_Adipose_unmeth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(jitter(validDatosADS_Adipose_halfmeth$methCoef),jitter(validDatosADS_Adipose_halfmeth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(jitter(validDatosADS_Adipose_meth$methCoef),jitter(validDatosADS_Adipose_meth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)


# Plot the points without noise
# points(validDatosADS_Adipose_unmeth$methCoef,validDatosADS_Adipose_unmeth$SD,pch=16,col=rgb(0,1,0,0.2),cex=1)
# points(validDatosADS_Adipose_halfmeth$methCoef,validDatosADS_Adipose_halfmeth$SD,pch=16,col=rgb(0,0,1,0.2),cex=1)
# points(validDatosADS_Adipose_meth$methCoef,validDatosADS_Adipose_meth$SD,pch=16,col=rgb(1,0,0,0.2),cex=1)

if (numbered=='TRUE'){
# Adding the number of points
text(0.1,0.4,dim(validDatosADS_Adipose_unmeth)[1],cex=1.2,col="red")
text(0.5,0.4,dim(validDatosADS_Adipose_halfmeth)[1],cex=1.2,col="red")
text(0.9,0.4,dim(validDatosADS_Adipose_meth)[1],cex=1.2,col="red")
}else{}
abline(v=0.25,lty=2,lwd=1.5)
abline(v=0.75,lty=2,lwd=1.5)
par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        )
par(mar=c(5,0,1,1))

barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()

####################################
## Without Informativeness filter ##
####################################
validDatosADS<-subset(datos_ADS,datos_ADS$methCoef != -1.00)
validDatosADS_Adipose<-subset(datos_ADS_Adipose,datos_ADS_Adipose$methCoef != -1.00)

validDatosADS_unmeth<-subset(validDatosADS, 
                            validDatosADS$methCoef<=0.25)
validDatosADS_halfmeth<-subset(validDatosADS, 
                              validDatosADS$methCoef>0.25 & validDatosADS$methCoef<0.75)
validDatosADS_meth<-subset(validDatosADS, 
                          validDatosADS$methCoef>=0.75)
validDatosADS_Adipose_unmeth<-subset(validDatosADS_Adipose, 
                               validDatosADS_Adipose$methCoef<=0.25)
validDatosADS_Adipose_halfmeth<-subset(validDatosADS_Adipose, 
                                 validDatosADS_Adipose$methCoef>0.25 & validDatosADS_Adipose$methCoef<0.75)
validDatosADS_Adipose_meth<-subset(validDatosADS_Adipose, 
                             validDatosADS_Adipose$methCoef>=0.75)

#Plotting

#ADS
if (numbered=='TRUE'){
  png(paste(resultsDIR,"figure23ADSAllBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure23ADSAllBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)

par(lwd=3)
par(cex.axis=1.2)

xhist <- hist(validDatosADS$methCoef, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(validDatosADS$SD, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)


par(mar=c(5,5,1,1))
plot(validDatosADS$methCoef,validDatosADS$SD,
     xlab="",ylab="",pch=20,xlim=c(0,1),col=rgb(0,0,0,0.2),type="n",ylim=c(0,0.6))

## Plot the points with some noise (jitter)
points(jitter(validDatosADS_unmeth$methCoef),jitter(validDatosADS_unmeth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(jitter(validDatosADS_halfmeth$methCoef),jitter(validDatosADS_halfmeth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(jitter(validDatosADS_meth$methCoef),jitter(validDatosADS_meth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# Plot the points without noise
# points(validDatosADS_unmeth$methCoef,validDatosADS_unmeth$SD,pch=16,col=rgb(0,1,0,0.2),cex=1)
# points(validDatosADS_halfmeth$methCoef,validDatosADS_halfmeth$SD,pch=16,col=rgb(0,0,1,0.2),cex=1)
# points(validDatosADS_meth$methCoef,validDatosADS_meth$SD,pch=16,col=rgb(1,0,0,0.2),cex=1)

if(numbered=='TRUE'){
# Adding the number of points
text(0.1,0.4,dim(validDatosADS_unmeth)[1],cex=1.2,col="red")
text(0.5,0.4,dim(validDatosADS_halfmeth)[1],cex=1.2,col="red")
text(0.9,0.4,dim(validDatosADS_meth)[1],cex=1.2,col="red")
}else{}
abline(v=0.25,lty=2,lwd=1.5)
abline(v=0.75,lty=2,lwd=1.5)
par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        )
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)

dev.off()

#ADS_Adipose

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure23ADS_AdiposeAllBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure23ADS_AdiposeAllBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)

par(lwd=3)
par(cex.axis=1.2)

xhist <- hist(validDatosADS_Adipose$methCoef, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(validDatosADS_Adipose$SD, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)


par(mar=c(5,5,1,1))
plot(validDatosADS_Adipose$methCoef,validDatosADS_Adipose$SD,
     xlab="",ylab="",pch=20,xlim=c(0,1),col=rgb(0,0,0,0.2),type="n",ylim=c(0,0.6))

# Plot the points with some noise (jitter)
points(jitter(validDatosADS_Adipose_unmeth$methCoef),jitter(validDatosADS_Adipose_unmeth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(jitter(validDatosADS_Adipose_halfmeth$methCoef),jitter(validDatosADS_Adipose_halfmeth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(jitter(validDatosADS_Adipose_meth$methCoef),jitter(validDatosADS_Adipose_meth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# Plot the points without noise
# points(validDatosADS_Adipose_unmeth$methCoef,validDatosADS_Adipose_unmeth$SD,pch=16,col=rgb(0,1,0,0.2),cex=1)
# points(validDatosADS_Adipose_halfmeth$methCoef,validDatosADS_Adipose_halfmeth$SD,pch=16,col=rgb(0,0,1,0.2),cex=1)
# points(validDatosADS_Adipose_meth$methCoef,validDatosADS_Adipose_meth$SD,pch=16,col=rgb(1,0,0,0.2),cex=1)

if(numbered=='TRUE'){
# Adding the number of points
text(0.1,0.4,dim(validDatosADS_Adipose_unmeth)[1],cex=1.2,col="red")
text(0.5,0.4,dim(validDatosADS_Adipose_halfmeth)[1],cex=1.2,col="red")
text(0.9,0.4,dim(validDatosADS_Adipose_meth)[1],cex=1.2,col="red")
}else{}
abline(v=0.25,lty=2,lwd=1.5)
abline(v=0.75,lty=2,lwd=1.5)
par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        )
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)

dev.off()
