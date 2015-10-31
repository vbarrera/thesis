#################
### Figure 24 ###
#################
h1CgiFolder<-paste(baseDIR,"oldOrganization/Cstatus/output/h1/CGi/",sep="")
imr90CgiFolder<-paste(baseDIR,"oldOrganization/Cstatus/output/imr90/CGi/",sep="")


datos_h1<-read.table(paste(h1CgiFolder,"hg18_h1_CGi_status_sd_5READ.txt",sep=""))
datos_imr90<-read.table(paste(imr90CgiFolder,"hg18_imr90_CGi_status_sd_5READ.txt",sep=""))
colnames(datos_h1)<-c("id","nCGt","chr","chrStart","chrEnd","CGinf","methCoef","SD","nReads")
colnames(datos_imr90)<-c("id","nCGt","chr","chrStart","chrEnd","CGinf","methCoef","SD","nReads")

ADSCgiFolder<-paste(baseDIR,"results/20120802/",sep="")
ADS_AdiposeCgiFolder<-paste(baseDIR,"results/20120802/",sep="")


datos_ADS<-read.table(paste(ADSCgiFolder,"hg18_reads_ads_CGi_status_sd.txt",sep=""))
datos_ADS_Adipose<-read.table(paste(ADS_AdiposeCgiFolder,"hg18_reads_ads_adipose_CGi_status_sd.txt",sep=""))
colnames(datos_ADS)<-c("id","nCGt","chr","chrStart","chrEnd","CGinf","methCoef","SD","nReads")
colnames(datos_ADS_Adipose)<-c("id","nCGt","chr","chrStart","chrEnd","CGinf","methCoef","SD","nReads")

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

png(paste(resultsDIR,"figure24Filtered.png",sep=""),height=12,width=12,units="cm",res=300)
def.par <- par(no.readonly = TRUE)

par(mfrow=c(2,2))
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosADS$SD), xlim=c(0,1), 
     ylim=c(0,20),type="n",xlab="",main="",ylab="")
lines(density(validDatosH1$SD),col="darkblue",lwd=0.8)
lines(density(validDatosIMR90$SD),col="red",lwd=0.8)
lines(density(validDatosADS$SD),col="green",lwd=0.8)
lines(density(validDatosADS_Adipose$SD),col="purple",lwd=0.8)

plot(density(validDatosADS_unmeth$SD), xlim=c(0,1), 
     ylim=c(0,20),type="n",xlab="",main="",ylab="")
lines(density(validDatosH1_unmeth$SD),col="darkblue",lwd=0.8)
lines(density(validDatosIMR90_unmeth$SD),col="red",lwd=0.8)
lines(density(validDatosADS_unmeth$SD,width=0.05),col="green",lwd=0.8)
lines(density(validDatosADS_Adipose_unmeth$SD,width=0.05),col="purple",lwd=0.8)

plot(density(validDatosADS_meth$SD), xlim=c(0,1), 
     ylim=c(0,20),type="n",xlab="",main="",ylab="")
lines(density(validDatosH1_meth$SD),col="darkblue",lwd=0.8)
lines(density(validDatosIMR90_meth$SD),col="red",lwd=0.8)
lines(density(validDatosADS_meth$SD),col="green",lwd=0.8)
lines(density(validDatosADS_Adipose_meth$SD),col="purple",lwd=0.8)

plot(density(validDatosADS_halfmeth$SD), xlim=c(0,1), 
     ylim=c(0,20),type="n",xlab="",main="",ylab="")
lines(density(validDatosH1_halfmeth$SD),col="darkblue",lwd=0.8)
lines(density(validDatosIMR90_halfmeth$SD),col="red",lwd=0.8)
lines(density(validDatosADS_halfmeth$SD),col="green",lwd=0.8)
lines(density(validDatosADS_Adipose_halfmeth$SD),col="purple",lwd=0.8)

par(mfrow=c(1,1))
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

png(paste(resultsDIR,"figure24All.png",sep=""),height=12,width=12,units="cm",res=300)
def.par <- par(no.readonly = TRUE)

par(mfrow=c(2,2))
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosADS$SD), xlim=c(0,1), 
     ylim=c(0,20),type="n",xlab="",main="",ylab="")
lines(density(validDatosH1$SD),col="darkblue",lwd=0.8)
lines(density(validDatosIMR90$SD),col="red",lwd=0.8)
lines(density(validDatosADS$SD),col="green",lwd=0.8)
lines(density(validDatosADS_Adipose$SD),col="purple",lwd=0.8)

plot(density(validDatosADS_unmeth$SD), xlim=c(0,1), 
     ylim=c(0,20),type="n",xlab="",main="",ylab="")
lines(density(validDatosH1_unmeth$SD),col="darkblue",lwd=0.8)
lines(density(validDatosIMR90_unmeth$SD),col="red",lwd=0.8)
lines(density(validDatosADS_unmeth$SD,width=0.05),col="green",lwd=0.8)
lines(density(validDatosADS_Adipose_unmeth$SD,width=0.05),col="purple",lwd=0.8)

plot(density(validDatosADS_meth$SD), xlim=c(0,1), 
     ylim=c(0,20),type="n",xlab="",main="",ylab="")
lines(density(validDatosH1_meth$SD),col="darkblue",lwd=0.8)
lines(density(validDatosIMR90_meth$SD),col="red",lwd=0.8)
lines(density(validDatosADS_meth$SD),col="green",lwd=0.8)
lines(density(validDatosADS_Adipose_meth$SD),col="purple",lwd=0.8)

plot(density(validDatosADS_halfmeth$SD), xlim=c(0,1), 
     ylim=c(0,20),type="n",xlab="",main="",ylab="")
lines(density(validDatosH1_halfmeth$SD),col="darkblue",lwd=0.8)
lines(density(validDatosIMR90_halfmeth$SD),col="red",lwd=0.8)
lines(density(validDatosADS_halfmeth$SD),col="green",lwd=0.8)
lines(density(validDatosADS_Adipose_halfmeth$SD),col="purple",lwd=0.8)

par(mfrow=c(1,1))
par(def.par)
dev.off()
