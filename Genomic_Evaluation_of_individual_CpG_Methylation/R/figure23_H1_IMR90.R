################
### Figure 23 ###
################

h1CgiFolder<-paste(baseDIR,"oldOrganization/Cstatus/output/h1/CGi/",sep="")
imr90CgiFolder<-paste(baseDIR,"oldOrganization/Cstatus/output/imr90/CGi/",sep="")


datos_h1<-read.table(paste(h1CgiFolder,"hg18_h1_CGi_status_sd_5READ.txt",sep=""))
datos_imr90<-read.table(paste(imr90CgiFolder,"hg18_imr90_CGi_status_sd_5READ.txt",sep=""))
colnames(datos_h1)<-c("id","nCGt","chr","chrStart","chrEnd","CGinf","methCoef","SD","nReads")
colnames(datos_imr90)<-c("id","nCGt","chr","chrStart","chrEnd","CGinf","methCoef","SD","nReads")

#################################
## With Informativeness filter ##
#################################
validDatosH1<-subset(datos_h1,datos_h1$methCoef != -1.00 & 
  (datos_h1$CGinf/2)/datos_h1$nCGt >=0.25)
validDatosIMR90<-subset(datos_imr90,datos_imr90$methCoef != -1.00 & 
  (datos_imr90$CGinf/2)/datos_imr90$nCGt >=0.25)
validDatosH1_unmeth<-subset(validDatosH1, 
                            validDatosH1$methCoef<=0.25)
validDatosH1_halfmeth<-subset(validDatosH1, 
                              validDatosH1$methCoef>0.25 & validDatosH1$methCoef<0.75)
validDatosH1_meth<-subset(validDatosH1, 
                          validDatosH1$methCoef>=0.75)
validDatosIMR90_unmeth<-subset(validDatosIMR90, 
                               validDatosIMR90$methCoef<=0.25)
validDatosIMR90_halfmeth<-subset(validDatosIMR90, 
                                 validDatosIMR90$methCoef>0.25 & validDatosIMR90$methCoef<0.75)
validDatosIMR90_meth<-subset(validDatosIMR90, 
                             validDatosIMR90$methCoef>=0.75)

#Plotting

#H1

if (numbered=='TRUE'){
  png(paste(resultsDIR,"figure23_H1FilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure23_H1FilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)

par(lwd=3)
par(cex.axis=1.2)

xhist <- hist(validDatosH1$methCoef, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(validDatosH1$SD, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)


par(mar=c(5,5,1,1))
plot(validDatosH1$methCoef,validDatosH1$SD,
     xlab="",ylab="",pch=20,xlim=c(0,1),col=rgb(0,0,0,0.2),type="n",ylim=c(0,0.6))

## Plot the points with some noise (jitter)
points(jitter(validDatosH1_unmeth$methCoef),jitter(validDatosH1_unmeth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(jitter(validDatosH1_halfmeth$methCoef),jitter(validDatosH1_halfmeth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(jitter(validDatosH1_meth$methCoef),jitter(validDatosH1_meth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# Plot the points without noise
# points(validDatosH1_unmeth$methCoef,validDatosH1_unmeth$SD,pch=16,col=rgb(0,1,0,0.2),cex=1)
# points(validDatosH1_halfmeth$methCoef,validDatosH1_halfmeth$SD,pch=16,col=rgb(0,0,1,0.2),cex=1)
# points(validDatosH1_meth$methCoef,validDatosH1_meth$SD,pch=16,col=rgb(1,0,0,0.2),cex=1)


if (numbered=='TRUE'){
# Adding the number of points
text(0.1,0.4,dim(validDatosH1_unmeth)[1],cex=1.2,col="red")
text(0.5,0.4,dim(validDatosH1_halfmeth)[1],cex=1.2,col="red")
text(0.9,0.4,dim(validDatosH1_meth)[1],cex=1.2,col="red")
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


#IMR90 

if (numbered=='TRUE'){
png(paste(resultsDIR,"figure23_IMR90FilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
  png(paste(resultsDIR,"figure23_IMR90FilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)

par(lwd=3)
par(cex.axis=1.2)

xhist <- hist(validDatosIMR90$methCoef, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(validDatosIMR90$SD, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot(validDatosIMR90$methCoef,validDatosIMR90$SD,
     xlab="", ylab="",pch=20,xlim=c(0,1),col=rgb(0,0,0,0.2),main="",type="n",ylim=c(0,0.6))

## Plot the points with some noise (jitter)
points(jitter(validDatosIMR90_unmeth$methCoef),jitter(validDatosIMR90_unmeth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(jitter(validDatosIMR90_halfmeth$methCoef),jitter(validDatosIMR90_halfmeth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(jitter(validDatosIMR90_meth$methCoef),jitter(validDatosIMR90_meth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)


# Plot the points without noise
# points(validDatosIMR90_unmeth$methCoef,validDatosIMR90_unmeth$SD,pch=16,col=rgb(0,1,0,0.2),cex=1)
# points(validDatosIMR90_halfmeth$methCoef,validDatosIMR90_halfmeth$SD,pch=16,col=rgb(0,0,1,0.2),cex=1)
# points(validDatosIMR90_meth$methCoef,validDatosIMR90_meth$SD,pch=16,col=rgb(1,0,0,0.2),cex=1)

if (numbered=='TRUE'){
# Adding the number of points
text(0.1,0.4,dim(validDatosIMR90_unmeth)[1],cex=1.2,col="red")
text(0.5,0.4,dim(validDatosIMR90_halfmeth)[1],cex=1.2,col="red")
text(0.9,0.4,dim(validDatosIMR90_meth)[1],cex=1.2,col="red")
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
validDatosH1<-subset(datos_h1,datos_h1$methCoef != -1.00)
validDatosIMR90<-subset(datos_imr90,datos_imr90$methCoef != -1.00)

validDatosH1_unmeth<-subset(validDatosH1, 
                            validDatosH1$methCoef<=0.25)
validDatosH1_halfmeth<-subset(validDatosH1, 
                              validDatosH1$methCoef>0.25 & validDatosH1$methCoef<0.75)
validDatosH1_meth<-subset(validDatosH1, 
                          validDatosH1$methCoef>=0.75)
validDatosIMR90_unmeth<-subset(validDatosIMR90, 
                               validDatosIMR90$methCoef<=0.25)
validDatosIMR90_halfmeth<-subset(validDatosIMR90, 
                                 validDatosIMR90$methCoef>0.25 & validDatosIMR90$methCoef<0.75)
validDatosIMR90_meth<-subset(validDatosIMR90, 
                             validDatosIMR90$methCoef>=0.75)

#Plotting

#H1
if (numbered=='TRUE'){
  png(paste(resultsDIR,"figure23_H1AllBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure23_H1AllBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)

par(lwd=3)
par(cex.axis=1.2)

xhist <- hist(validDatosH1$methCoef, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(validDatosH1$SD, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)


par(mar=c(5,5,1,1))
plot(validDatosH1$methCoef,validDatosH1$SD,
     xlab="",ylab="",pch=20,xlim=c(0,1),col=rgb(0,0,0,0.2),type="n",ylim=c(0,0.6))

## Plot the points with some noise (jitter)
points(jitter(validDatosH1_unmeth$methCoef),jitter(validDatosH1_unmeth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(jitter(validDatosH1_halfmeth$methCoef),jitter(validDatosH1_halfmeth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(jitter(validDatosH1_meth$methCoef),jitter(validDatosH1_meth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# Plot the points without noise
# points(validDatosH1_unmeth$methCoef,validDatosH1_unmeth$SD,pch=16,col=rgb(0,1,0,0.2),cex=1)
# points(validDatosH1_halfmeth$methCoef,validDatosH1_halfmeth$SD,pch=16,col=rgb(0,0,1,0.2),cex=1)
# points(validDatosH1_meth$methCoef,validDatosH1_meth$SD,pch=16,col=rgb(1,0,0,0.2),cex=1)

if(numbered=='TRUE'){
# Adding the number of points
text(0.1,0.4,dim(validDatosH1_unmeth)[1],cex=1.2,col="red")
text(0.5,0.4,dim(validDatosH1_halfmeth)[1],cex=1.2,col="red")
text(0.9,0.4,dim(validDatosH1_meth)[1],cex=1.2,col="red")
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

#IMR90

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure23_IMR90AllBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure23_IMR90AllBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)

par(lwd=3)
par(cex.axis=1.2)

xhist <- hist(validDatosIMR90$methCoef, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(validDatosIMR90$SD, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)


par(mar=c(5,5,1,1))
plot(validDatosIMR90$methCoef,validDatosIMR90$SD,
     xlab="",ylab="",pch=20,xlim=c(0,1),col=rgb(0,0,0,0.2),type="n",ylim=c(0,0.6))

# Plot the points with some noise (jitter)
points(jitter(validDatosIMR90_unmeth$methCoef),jitter(validDatosIMR90_unmeth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(jitter(validDatosIMR90_halfmeth$methCoef),jitter(validDatosIMR90_halfmeth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(jitter(validDatosIMR90_meth$methCoef),jitter(validDatosIMR90_meth$SD,factor=5),pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# Plot the points without noise
# points(validDatosIMR90_unmeth$methCoef,validDatosIMR90_unmeth$SD,pch=16,col=rgb(0,1,0,0.2),cex=1)
# points(validDatosIMR90_halfmeth$methCoef,validDatosIMR90_halfmeth$SD,pch=16,col=rgb(0,0,1,0.2),cex=1)
# points(validDatosIMR90_meth$methCoef,validDatosIMR90_meth$SD,pch=16,col=rgb(1,0,0,0.2),cex=1)

if(numbered=='TRUE'){
# Adding the number of points
text(0.1,0.4,dim(validDatosIMR90_unmeth)[1],cex=1.2,col="red")
text(0.5,0.4,dim(validDatosIMR90_halfmeth)[1],cex=1.2,col="red")
text(0.9,0.4,dim(validDatosIMR90_meth)[1],cex=1.2,col="red")
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
