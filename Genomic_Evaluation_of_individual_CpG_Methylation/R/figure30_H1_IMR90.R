#################
### Figure 30 ###
#################

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


library(Epi)

#Methylated ROC


groupsH1<-as.numeric(selected_reDinucleotide_Element_H1$E_methCoef>=0.75)
groupsIMR90<-as.numeric(selected_reDinucleotide_Element_IMR90$E_methCoef>=0.75)

predictionsH1<-selected_reDinucleotide_Element_H1$meanMethCoef_H1
predictionsIMR90<-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90

png(paste(resultsDIR,"figure30_MethH1FilteredData.png",sep=""))
q<-ROC(test=predictionsH1,stat=groupsH1,plot=c("ROC"),PV=TRUE,MX=TRUE,
       MI=TRUE,AUC=TRUE)
dev.off()

q<-ROC(test=predictionsH1,stat=groupsH1,plot=c("ROC"),PV=FALSE,MX=FALSE,
       MI=FALSE,AUC=FALSE)

png(paste(resultsDIR,"figure30_MethH1Filtered.png",sep=""),height=8,width=8,units="cm",res=600)
par(lwd=1.5)
par(cex.axis=0.8)
plot(1-q$res$spec,q$res$sens,type="n",lwd=1,xlab="",ylab="")
lines(1-q$res$spec,q$res$sens,lwd=1)
abline(0.0,1,lty=1,lwd=0.5,col="grey")
dev.off()

png(paste(resultsDIR,"figure30_MethIMR90FilteredData.png",sep=""))
q<-ROC(test=predictionsIMR90,stat=groupsIMR90,plot=c("ROC"),PV=TRUE,MX=TRUE,
       MI=TRUE,AUC=TRUE)
dev.off()

q<-ROC(test=predictionsIMR90,stat=groupsIMR90,plot=c("ROC"),PV=FALSE,MX=FALSE,
       MI=FALSE,AUC=FALSE)

png(paste(resultsDIR,"figure30_MethIMR90Filtered.png",sep=""),height=8,width=8,units="cm",res=600)
par(lwd=1.5)
par(cex.axis=0.8)
plot(1-q$res$spec,q$res$sens,type="n",lwd=1,xlab="",ylab="")
lines(1-q$res$spec,q$res$sens,lwd=1)
abline(0.0,1,lty=1,lwd=0.5,col="grey")
dev.off()

#UnMethylated ROC

groupsH1<-as.numeric(selected_reDinucleotide_Element_H1$E_methCoef>=0.25)
groupsIMR90<-as.numeric(selected_reDinucleotide_Element_IMR90$E_methCoef>=0.25)

predictionsH1<-selected_reDinucleotide_Element_H1$meanMethCoef_H1
predictionsIMR90<-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90

png(paste(resultsDIR,"figure30_UnmethH1FilteredData.png",sep=""))
q<-ROC(test=predictionsH1,stat=groupsH1,plot=c("ROC"),PV=TRUE,MX=TRUE,
       MI=TRUE,AUC=TRUE)
dev.off()

q<-ROC(test=predictionsH1,stat=groupsH1,plot=c("ROC"),PV=FALSE,MX=FALSE,
       MI=FALSE,AUC=FALSE)

png(paste(resultsDIR,"figure30_UnmethH1Filtered.png",sep=""),height=8,width=8,units="cm",res=600)
par(lwd=1.5)
par(cex.axis=0.8)
plot(1-q$res$spec,q$res$sens,type="n",lwd=1,xlab="",ylab="")
lines(1-q$res$spec,q$res$sens,lwd=1)
abline(0.0,1,lty=1,lwd=0.5,col="grey")
dev.off()

png(paste(resultsDIR,"figure30_UnMethIMR90FilteredData.png",sep=""))
q<-ROC(test=predictionsIMR90,stat=groupsIMR90,plot=c("ROC"),PV=TRUE,MX=TRUE,
       MI=TRUE,AUC=TRUE)
dev.off()

q<-ROC(test=predictionsIMR90,stat=groupsIMR90,plot=c("ROC"),PV=FALSE,MX=FALSE,
       MI=FALSE,AUC=FALSE)

png(paste(resultsDIR,"figure30_UnmethIMR90Filtered.png",sep=""),height=8,width=8,units="cm",res=600)
par(lwd=1.5)
par(cex.axis=0.8)
plot(1-q$res$spec,q$res$sens,type="n",lwd=1,xlab="",ylab="")
lines(1-q$res$spec,q$res$sens,lwd=1)
abline(0.0,1,lty=1,lwd=0.5,col="grey")
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

#Methylated ROC


groupsH1<-as.numeric(selected_reDinucleotide_Element_H1$E_methCoef>=0.75)
groupsIMR90<-as.numeric(selected_reDinucleotide_Element_IMR90$E_methCoef>=0.75)

predictionsH1<-selected_reDinucleotide_Element_H1$meanMethCoef_H1
predictionsIMR90<-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90

png(paste(resultsDIR,"figure30_MethH1AllData.png",sep=""))
q<-ROC(test=predictionsH1,stat=groupsH1,plot=c("ROC"),PV=TRUE,MX=TRUE,
       MI=TRUE,AUC=TRUE)
dev.off()

q<-ROC(test=predictionsH1,stat=groupsH1,plot=c("ROC"),PV=FALSE,MX=FALSE,
       MI=FALSE,AUC=FALSE)

png(paste(resultsDIR,"figure30_MethH1All.png",sep=""),height=8,width=8,units="cm",res=600)
par(lwd=1.5)
par(cex.axis=0.8)
plot(1-q$res$spec,q$res$sens,type="n",lwd=1,xlab="",ylab="")
lines(1-q$res$spec,q$res$sens,lwd=1)
abline(0.0,1,lty=1,lwd=0.5,col="grey")
dev.off()

png(paste(resultsDIR,"figure30_MethIMR90AllData.png",sep=""))
q<-ROC(test=predictionsIMR90,stat=groupsIMR90,plot=c("ROC"),PV=TRUE,MX=TRUE,
       MI=TRUE,AUC=TRUE)
dev.off()

q<-ROC(test=predictionsIMR90,stat=groupsIMR90,plot=c("ROC"),PV=FALSE,MX=FALSE,
       MI=FALSE,AUC=FALSE)

png(paste(resultsDIR,"figure30_MethIMR90All.png",sep=""),height=8,width=8,units="cm",res=600)
par(lwd=1.5)
par(cex.axis=0.8)
plot(1-q$res$spec,q$res$sens,type="n",lwd=1,xlab="",ylab="")
lines(1-q$res$spec,q$res$sens,lwd=1)
abline(0.0,1,lty=1,lwd=0.5,col="grey")
dev.off()

#UnMethylated ROC

groupsH1<-as.numeric(selected_reDinucleotide_Element_H1$E_methCoef>=0.25)
groupsIMR90<-as.numeric(selected_reDinucleotide_Element_IMR90$E_methCoef>=0.25)

predictionsH1<-selected_reDinucleotide_Element_H1$meanMethCoef_H1
predictionsIMR90<-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90

png(paste(resultsDIR,"figure30_UnmethH1AllData.png",sep=""))
q<-ROC(test=predictionsH1,stat=groupsH1,plot=c("ROC"),PV=TRUE,MX=TRUE,
       MI=TRUE,AUC=TRUE)
dev.off()

q<-ROC(test=predictionsH1,stat=groupsH1,plot=c("ROC"),PV=FALSE,MX=FALSE,
       MI=FALSE,AUC=FALSE)

png(paste(resultsDIR,"figure30_UnmethH1All.png",sep=""),height=8,width=8,units="cm",res=600)
par(lwd=1.5)
par(cex.axis=0.8)
plot(1-q$res$spec,q$res$sens,type="n",lwd=1,xlab="",ylab="")
lines(1-q$res$spec,q$res$sens,lwd=1)
abline(0.0,1,lty=1,lwd=0.5,col="grey")
dev.off()

png(paste(resultsDIR,"figure30_UnMethIMR90AllData.png",sep=""))
q<-ROC(test=predictionsIMR90,stat=groupsIMR90,plot=c("ROC"),PV=TRUE,MX=TRUE,
       MI=TRUE,AUC=TRUE)
dev.off()

q<-ROC(test=predictionsIMR90,stat=groupsIMR90,plot=c("ROC"),PV=FALSE,MX=FALSE,
       MI=FALSE,AUC=FALSE)

png(paste(resultsDIR,"figure30_UnmethIMR90All.png",sep=""),height=8,width=8,units="cm",res=600)
par(lwd=1.5)
par(cex.axis=0.8)
plot(1-q$res$spec,q$res$sens,type="n",lwd=1,xlab="",ylab="")
lines(1-q$res$spec,q$res$sens,lwd=1)
abline(0.0,1,lty=1,lwd=0.5,col="grey")
dev.off()

