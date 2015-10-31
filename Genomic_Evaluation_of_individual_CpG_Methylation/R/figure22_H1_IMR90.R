#################
### figure 22 ###
#################

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

# H1

png(paste(resultsDIR,"figure22_H1Filtered.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosH1$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,4.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)
par(def.par)
dev.off()

png(paste(resultsDIR,"figure22_IMR90Filtered.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosIMR90$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,4.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)

par(def.par)
dev.off

####################################
## Without Informativeness filter ##
####################################
validDatosH1<-subset(datos_h1,datos_h1$methCoef != -1.00)
validDatosIMR90<-subset(datos_imr90,datos_imr90$methCoef != -1.00)


# H1

png(paste(resultsDIR,"figure22_H1All.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosH1$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,5.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)
par(def.par)
dev.off()

png(paste(resultsDIR,"figure22_IMR90All.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosIMR90$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,5.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)

par(def.par)
dev.off()
