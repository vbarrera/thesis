#################
### Figure 30 ###
#################

#################################
## With Informativeness filter ##
#################################

selected_reDinucleotide_Element_ADS<-subset(reDinucleotide_Element_ADS,
                                           reDinucleotide_Element_ADS$meanMethCoef_ADS!=-1 & 
                                             reDinucleotide_Element_ADS$E_methCoef!=-1 & ((reDinucleotide_Element_ADS$E_posInf/2)/reDinucleotide_Element_ADS$E_nCG)>=0.25)
selected_reDinucleotide_Element_ADS_Adipose<-subset(reDinucleotide_Element_ADS_Adipose,
                                              reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose!=-1 & 
                                                reDinucleotide_Element_ADS_Adipose$E_methCoef!=-1 & ((reDinucleotide_Element_ADS_Adipose$E_posInf/2)/reDinucleotide_Element_ADS_Adipose$E_nCG)>=0.25)
##########################


library(Epi)

#Methylated ROC


groupsADS<-as.numeric(selected_reDinucleotide_Element_ADS$E_methCoef>=0.75)
groupsADS_Adipose<-as.numeric(selected_reDinucleotide_Element_ADS_Adipose$E_methCoef>=0.75)

predictionsADS<-selected_reDinucleotide_Element_ADS$meanMethCoef_ADS
predictionsADS_Adipose<-selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose

png(paste(resultsDIR,"figure30_MethADSFilteredData.png",sep=""))
q<-ROC(test=predictionsADS,stat=groupsADS,plot=c("ROC"),PV=TRUE,MX=TRUE,
       MI=TRUE,AUC=TRUE)
dev.off()

q<-ROC(test=predictionsADS,stat=groupsADS,plot=c("ROC"),PV=FALSE,MX=FALSE,
       MI=FALSE,AUC=FALSE)

png(paste(resultsDIR,"figure30_MethADSFiltered.png",sep=""),height=8,width=8,units="cm",res=600)
par(lwd=1.5)
par(cex.axis=0.8)
plot(1-q$res$spec,q$res$sens,type="n",lwd=1,xlab="",ylab="")
lines(1-q$res$spec,q$res$sens,lwd=1)
abline(0.0,1,lty=1,lwd=0.5,col="grey")
dev.off()

png(paste(resultsDIR,"figure30_MethADS_AdiposeFilteredData.png",sep=""))
q<-ROC(test=predictionsADS_Adipose,stat=groupsADS_Adipose,plot=c("ROC"),PV=TRUE,MX=TRUE,
       MI=TRUE,AUC=TRUE)
dev.off()

q<-ROC(test=predictionsADS_Adipose,stat=groupsADS_Adipose,plot=c("ROC"),PV=FALSE,MX=FALSE,
       MI=FALSE,AUC=FALSE)

png(paste(resultsDIR,"figure30_MethADS_AdiposeFiltered.png",sep=""),height=8,width=8,units="cm",res=600)
par(lwd=1.5)
par(cex.axis=0.8)
plot(1-q$res$spec,q$res$sens,type="n",lwd=1,xlab="",ylab="")
lines(1-q$res$spec,q$res$sens,lwd=1)
abline(0.0,1,lty=1,lwd=0.5,col="grey")
dev.off()

#UnMethylated ROC

groupsADS<-as.numeric(selected_reDinucleotide_Element_ADS$E_methCoef>=0.25)
groupsADS_Adipose<-as.numeric(selected_reDinucleotide_Element_ADS_Adipose$E_methCoef>=0.25)

predictionsADS<-selected_reDinucleotide_Element_ADS$meanMethCoef_ADS
predictionsADS_Adipose<-selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose

png(paste(resultsDIR,"figure30_UnmethADSFilteredData.png",sep=""))
q<-ROC(test=predictionsADS,stat=groupsADS,plot=c("ROC"),PV=TRUE,MX=TRUE,
       MI=TRUE,AUC=TRUE)
dev.off()

q<-ROC(test=predictionsADS,stat=groupsADS,plot=c("ROC"),PV=FALSE,MX=FALSE,
       MI=FALSE,AUC=FALSE)

png(paste(resultsDIR,"figure30_UnmethADSFiltered.png",sep=""),height=8,width=8,units="cm",res=600)
par(lwd=1.5)
par(cex.axis=0.8)
plot(1-q$res$spec,q$res$sens,type="n",lwd=1,xlab="",ylab="")
lines(1-q$res$spec,q$res$sens,lwd=1)
abline(0.0,1,lty=1,lwd=0.5,col="grey")
dev.off()

png(paste(resultsDIR,"figure30_UnMethADS_AdiposeFilteredData.png",sep=""))
q<-ROC(test=predictionsADS_Adipose,stat=groupsADS_Adipose,plot=c("ROC"),PV=TRUE,MX=TRUE,
       MI=TRUE,AUC=TRUE)
dev.off()

q<-ROC(test=predictionsADS_Adipose,stat=groupsADS_Adipose,plot=c("ROC"),PV=FALSE,MX=FALSE,
       MI=FALSE,AUC=FALSE)

png(paste(resultsDIR,"figure30_UnmethADS_AdiposeFiltered.png",sep=""),height=8,width=8,units="cm",res=600)
par(lwd=1.5)
par(cex.axis=0.8)
plot(1-q$res$spec,q$res$sens,type="n",lwd=1,xlab="",ylab="")
lines(1-q$res$spec,q$res$sens,lwd=1)
abline(0.0,1,lty=1,lwd=0.5,col="grey")
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

#######################

#Methylated ROC


groupsADS<-as.numeric(selected_reDinucleotide_Element_ADS$E_methCoef>=0.75)
groupsADS_Adipose<-as.numeric(selected_reDinucleotide_Element_ADS_Adipose$E_methCoef>=0.75)

predictionsADS<-selected_reDinucleotide_Element_ADS$meanMethCoef_ADS
predictionsADS_Adipose<-selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose

png(paste(resultsDIR,"figure30_MethADSAllData.png",sep=""))
q<-ROC(test=predictionsADS,stat=groupsADS,plot=c("ROC"),PV=TRUE,MX=TRUE,
       MI=TRUE,AUC=TRUE)
dev.off()

q<-ROC(test=predictionsADS,stat=groupsADS,plot=c("ROC"),PV=FALSE,MX=FALSE,
       MI=FALSE,AUC=FALSE)

png(paste(resultsDIR,"figure30_MethADSAll.png",sep=""),height=8,width=8,units="cm",res=600)
par(lwd=1.5)
par(cex.axis=0.8)
plot(1-q$res$spec,q$res$sens,type="n",lwd=1,xlab="",ylab="")
lines(1-q$res$spec,q$res$sens,lwd=1)
abline(0.0,1,lty=1,lwd=0.5,col="grey")
dev.off()

png(paste(resultsDIR,"figure30_MethADS_AdiposeAllData.png",sep=""))
q<-ROC(test=predictionsADS_Adipose,stat=groupsADS_Adipose,plot=c("ROC"),PV=TRUE,MX=TRUE,
       MI=TRUE,AUC=TRUE)
dev.off()

q<-ROC(test=predictionsADS_Adipose,stat=groupsADS_Adipose,plot=c("ROC"),PV=FALSE,MX=FALSE,
       MI=FALSE,AUC=FALSE)

png(paste(resultsDIR,"figure30_MethADS_AdiposeAll.png",sep=""),height=8,width=8,units="cm",res=600)
par(lwd=1.5)
par(cex.axis=0.8)
plot(1-q$res$spec,q$res$sens,type="n",lwd=1,xlab="",ylab="")
lines(1-q$res$spec,q$res$sens,lwd=1)
abline(0.0,1,lty=1,lwd=0.5,col="grey")
dev.off()

#UnMethylated ROC

groupsADS<-as.numeric(selected_reDinucleotide_Element_ADS$E_methCoef>=0.25)
groupsADS_Adipose<-as.numeric(selected_reDinucleotide_Element_ADS_Adipose$E_methCoef>=0.25)

predictionsADS<-selected_reDinucleotide_Element_ADS$meanMethCoef_ADS
predictionsADS_Adipose<-selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose

png(paste(resultsDIR,"figure30_UnmethADSAllData.png",sep=""))
q<-ROC(test=predictionsADS,stat=groupsADS,plot=c("ROC"),PV=TRUE,MX=TRUE,
       MI=TRUE,AUC=TRUE)
dev.off()

q<-ROC(test=predictionsADS,stat=groupsADS,plot=c("ROC"),PV=FALSE,MX=FALSE,
       MI=FALSE,AUC=FALSE)

png(paste(resultsDIR,"figure30_UnmethADSAll.png",sep=""),height=8,width=8,units="cm",res=600)
par(lwd=1.5)
par(cex.axis=0.8)
plot(1-q$res$spec,q$res$sens,type="n",lwd=1,xlab="",ylab="")
lines(1-q$res$spec,q$res$sens,lwd=1)
abline(0.0,1,lty=1,lwd=0.5,col="grey")
dev.off()

png(paste(resultsDIR,"figure30_UnMethADS_AdiposeAllData.png",sep=""))
q<-ROC(test=predictionsADS_Adipose,stat=groupsADS_Adipose,plot=c("ROC"),PV=TRUE,MX=TRUE,
       MI=TRUE,AUC=TRUE)
dev.off()

q<-ROC(test=predictionsADS_Adipose,stat=groupsADS_Adipose,plot=c("ROC"),PV=FALSE,MX=FALSE,
       MI=FALSE,AUC=FALSE)

png(paste(resultsDIR,"figure30_UnmethADS_AdiposeAll.png",sep=""),height=8,width=8,units="cm",res=600)
par(lwd=1.5)
par(cex.axis=0.8)
plot(1-q$res$spec,q$res$sens,type="n",lwd=1,xlab="",ylab="")
lines(1-q$res$spec,q$res$sens,lwd=1)
abline(0.0,1,lty=1,lwd=0.5,col="grey")
dev.off()

