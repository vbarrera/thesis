#################
### figure 22 ###
#################

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

# ADS

png(paste(resultsDIR,"figure22_ADSFiltered.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosADS$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,9.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)
par(def.par)
dev.off()

png(paste(resultsDIR,"figure22_ADS_AdiposeFiltered.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosADS_Adipose$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,9.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)

par(def.par)
dev.off()

####################################
## Without Informativeness filter ##
####################################
validDatosADS<-subset(datos_ADS,datos_ADS$methCoef != -1.00)
validDatosADS_Adipose<-subset(datos_ADS_Adipose,datos_ADS_Adipose$methCoef != -1.00)


# ADS

png(paste(resultsDIR,"figure22_ADSAll.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosADS$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,9.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)
par(def.par)
dev.off()

png(paste(resultsDIR,"figure22_ADS_AdiposeAll.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosADS_Adipose$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,9.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)

par(def.par)
dev.off()
#################
### figure 22 ###
#################

reads_ff_ipsc_19.11CgiFolder<-paste(baseDIR,"results/20120724/",sep="")
reads_ff_ipsc_19.7CgiFolder<-paste(baseDIR,"results/20120724/",sep="")


datos_reads_ff_ipsc_19.11<-read.table(paste(reads_ff_ipsc_19.11CgiFolder,"hg18_reads_ff_ipsc_19.11_CGi_status_sd.txt",sep=""))
datos_reads_ff_ipsc_19.7<-read.table(paste(reads_ff_ipsc_19.7CgiFolder,"hg18_reads_ff_ipsc_19.7_CGi_status_sd.txt",sep=""))
colnames(datos_reads_ff_ipsc_19.11)<-c("id","nCGt","chr","chrStart","chrEnd","CGinf","methCoef","SD","nReads")
colnames(datos_reads_ff_ipsc_19.7)<-c("id","nCGt","chr","chrStart","chrEnd","CGinf","methCoef","SD","nReads")


#################################
## With Informativeness filter ##
#################################
validDatosreads_ff_ipsc_19.11<-subset(datos_reads_ff_ipsc_19.11,datos_reads_ff_ipsc_19.11$methCoef != -1.00 & 
  (datos_reads_ff_ipsc_19.11$CGinf/2)/datos_reads_ff_ipsc_19.11$nCGt >=0.25)
validDatosreads_ff_ipsc_19.7<-subset(datos_reads_ff_ipsc_19.7,datos_reads_ff_ipsc_19.7$methCoef != -1.00 & 
  (datos_reads_ff_ipsc_19.7$CGinf/2)/datos_reads_ff_ipsc_19.7$nCGt >=0.25)

# reads_ff_ipsc_19.11

png(paste(resultsDIR,"figure22_reads_ff_ipsc_19.11Filtered.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosreads_ff_ipsc_19.11$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,9.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)
par(def.par)
dev.off()

png(paste(resultsDIR,"figure22_reads_ff_ipsc_19.7Filtered.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosreads_ff_ipsc_19.7$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,9.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)

par(def.par)
dev.off()

####################################
## Without Informativeness filter ##
####################################
validDatosreads_ff_ipsc_19.11<-subset(datos_reads_ff_ipsc_19.11,datos_reads_ff_ipsc_19.11$methCoef != -1.00)
validDatosreads_ff_ipsc_19.7<-subset(datos_reads_ff_ipsc_19.7,datos_reads_ff_ipsc_19.7$methCoef != -1.00)


# reads_ff_ipsc_19.11

png(paste(resultsDIR,"figure22_reads_ff_ipsc_19.11All.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosreads_ff_ipsc_19.11$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,9.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)
par(def.par)
dev.off()

png(paste(resultsDIR,"figure22_reads_ff_ipsc_19.7All.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosreads_ff_ipsc_19.7$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,9.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)

par(def.par)
dev.off()
#################
### figure 22 ###
#################

reads_ff_ipsc_6.9CgiFolder<-paste(baseDIR,"results/20120724/",sep="")
reads_h9CgiFolder<-paste(baseDIR,"results/20120724/",sep="")


datos_reads_ff_ipsc_6.9<-read.table(paste(reads_ff_ipsc_6.9CgiFolder,"hg18_reads_ff_ipsc_6.9_CGi_status_sd.txt",sep=""))
datos_reads_h9<-read.table(paste(reads_h9CgiFolder,"hg18_reads_h9_CGi_status_sd.txt",sep=""))
colnames(datos_reads_ff_ipsc_6.9)<-c("id","nCGt","chr","chrStart","chrEnd","CGinf","methCoef","SD","nReads")
colnames(datos_reads_h9)<-c("id","nCGt","chr","chrStart","chrEnd","CGinf","methCoef","SD","nReads")


#################################
## With Informativeness filter ##
#################################
validDatosreads_ff_ipsc_6.9<-subset(datos_reads_ff_ipsc_6.9,datos_reads_ff_ipsc_6.9$methCoef != -1.00 & 
  (datos_reads_ff_ipsc_6.9$CGinf/2)/datos_reads_ff_ipsc_6.9$nCGt >=0.25)
validDatosreads_h9<-subset(datos_reads_h9,datos_reads_h9$methCoef != -1.00 & 
  (datos_reads_h9$CGinf/2)/datos_reads_h9$nCGt >=0.25)

# reads_ff_ipsc_6.9

png(paste(resultsDIR,"figure22_reads_ff_ipsc_6.9Filtered.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosreads_ff_ipsc_6.9$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,9.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)
par(def.par)
dev.off()

png(paste(resultsDIR,"figure22_reads_h9Filtered.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosreads_h9$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,9.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)

par(def.par)
dev.off()

####################################
## Without Informativeness filter ##
####################################
validDatosreads_ff_ipsc_6.9<-subset(datos_reads_ff_ipsc_6.9,datos_reads_ff_ipsc_6.9$methCoef != -1.00)
validDatosreads_h9<-subset(datos_reads_h9,datos_reads_h9$methCoef != -1.00)


# reads_ff_ipsc_6.9

png(paste(resultsDIR,"figure22_reads_ff_ipsc_6.9All.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosreads_ff_ipsc_6.9$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,9.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)
par(def.par)
dev.off()

png(paste(resultsDIR,"figure22_reads_h9All.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosreads_h9$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,9.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)

par(def.par)
dev.off()
#################
### figure 22 ###
#################

reads_imr90_ipscCgiFolder<-paste(baseDIR,"results/20120724/",sep="")
reads_imr90_ipscCgiFolder<-paste(baseDIR,"results/20120724/",sep="")


datos_reads_imr90_ipsc<-read.table(paste(reads_imr90_ipscCgiFolder,"hg18_reads_imr90_ipsc_CGi_status_sd.txt",sep=""))
datos_reads_imr90_ipsc<-read.table(paste(reads_imr90_ipscCgiFolder,"hg18_reads_imr90_ipsc_CGi_status_sd.txt",sep=""))
colnames(datos_reads_imr90_ipsc)<-c("id","nCGt","chr","chrStart","chrEnd","CGinf","methCoef","SD","nReads")
colnames(datos_reads_imr90_ipsc)<-c("id","nCGt","chr","chrStart","chrEnd","CGinf","methCoef","SD","nReads")


#################################
## With Informativeness filter ##
#################################
validDatosreads_imr90_ipsc<-subset(datos_reads_imr90_ipsc,datos_reads_imr90_ipsc$methCoef != -1.00 & 
  (datos_reads_imr90_ipsc$CGinf/2)/datos_reads_imr90_ipsc$nCGt >=0.25)
validDatosreads_imr90_ipsc<-subset(datos_reads_imr90_ipsc,datos_reads_imr90_ipsc$methCoef != -1.00 & 
  (datos_reads_imr90_ipsc$CGinf/2)/datos_reads_imr90_ipsc$nCGt >=0.25)

# reads_imr90_ipsc

png(paste(resultsDIR,"figure22_reads_imr90_ipscFiltered.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosreads_imr90_ipsc$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,9.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)
par(def.par)
dev.off()

png(paste(resultsDIR,"figure22_reads_imr90_ipscFiltered.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosreads_imr90_ipsc$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,9.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)

par(def.par)
dev.off()

####################################
## Without Informativeness filter ##
####################################
validDatosreads_imr90_ipsc<-subset(datos_reads_imr90_ipsc,datos_reads_imr90_ipsc$methCoef != -1.00)
validDatosreads_imr90_ipsc<-subset(datos_reads_imr90_ipsc,datos_reads_imr90_ipsc$methCoef != -1.00)


# reads_imr90_ipsc

png(paste(resultsDIR,"figure22_reads_imr90_ipscAll.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosreads_imr90_ipsc$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,9.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)
par(def.par)
dev.off()

png(paste(resultsDIR,"figure22_reads_imr90_ipscAll.png",sep=""),height=8,width=8,units="cm",res=300)
def.par <- par(no.readonly = TRUE)
par(lwd=1.5)
par(cex.axis=0.8)
plot(density(validDatosreads_imr90_ipsc$methCoef),
     main="",xlab="",ylab="",xlim=c(0,1),ylim=c(0.0,9.0),lwd=1)
abline(v=0.25,lty=2,lwd=1)
abline(v=0.75,lty=2,lwd=1)

par(def.par)
dev.off()
