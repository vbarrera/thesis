################
### Figure 25 ###
################
require(RMySQL)
con<-dbConnect(MySQL(),dbname="MethData_Lister_hg18")

reNucleotide_Element<-dbGetQuery(con,"SELECT R.nucleotide,R.id_In_RE,R.position,MPA.methCoef,
MPA.name,MPA.nReads,E.id_In_Type,E.nCG,E.chrom,E.chromStart,E.chromEnd,MEA.methCoef,MEA.posInf,MEA.Std_Dev 
FROM (((R_POS R JOIN METH_POS_ASSIGNMENT MPA ON R.nucleotide=MPA.nucleotide 
AND R.RE_name=MPA.RE_name AND R.id_In_RE=MPA.id_In_RE) 
JOIN CORRESPONDENCE C ON R.nucleotide=C.nucleotide 
AND R.RE_name=C.RE_name AND R.id_In_RE=C.id_In_RE) 
JOIN ELEMENT E ON C.type=E.type AND C.id_In_Type=E.id_In_Type) 
JOIN METH_ELEM_ASSIGNMENT MEA ON MEA.type=E.type 
AND MEA.id_In_Type=E.id_In_Type 
WHERE MPA.name=MEA.name AND R.RE_name='HpaII' AND E.type='CpGisland' 
AND MPA.RE_name='HpaII' AND MEA.type='CpGisland' AND C.type='CpGisland'")

reDinucleotide_Element<-reshape(reNucleotide_Element,idvar=c("id_In_RE","name"),
                                direction="wide",timevar="nucleotide")

reDinucleotide_Element<-subset(reDinucleotide_Element,select=c(id_In_RE,name,position.C,nReads.C,methCoef.C,
                                                               methCoef.G,id_In_Type.C,chrom.C,chromStart.C,chromEnd.C,nCG.C,methCoef.1.C,posInf.C,Std_Dev.C)
)

colnames(reDinucleotide_Element)<-c("id_In_RE","CLine","RE_position","RE_nReads",
                                    "C_methCoef","G_methCoef","id_In_Type","E_chrom","E_chromStart","E_chromEnd","E_nCG",
                                    "E_methCoef","E_posInf","E_Std_Dev")

methCoefMean<-function(x){
  if (as.numeric(x[5])==-1 & as.numeric(x[6])!=-1){
    x[6]
  }
  else if (as.numeric(x[6])==-1 & as.numeric(x[5])!=-1){
    x[5]
  }
  else{
    (as.numeric(x[5])+as.numeric(x[6]))/2
  }
  
}

reDinucleotide_Element_H1<-subset(reDinucleotide_Element,reDinucleotide_Element$CLine=="H1")
reDinucleotide_Element_IMR90<-subset(reDinucleotide_Element,reDinucleotide_Element$CLine=="IMR90")

meanMethCoef_H1<-apply(reDinucleotide_Element_H1,1,FUN=methCoefMean)
meanMethCoef_IMR90<-apply(reDinucleotide_Element_IMR90,1,FUN=methCoefMean)

reDinucleotide_Element_H1<-cbind(reDinucleotide_Element_H1,meanMethCoef_H1=as.numeric(as.vector(meanMethCoef_H1)))
reDinucleotide_Element_IMR90<-cbind(reDinucleotide_Element_IMR90,meanMethCoef_IMR90=as.numeric(as.vector(meanMethCoef_IMR90)))

reDinucleotide_Element_ADS<-subset(reDinucleotide_Element,reDinucleotide_Element$CLine=="ADS")
reDinucleotide_Element_ADS_Adipose<-subset(reDinucleotide_Element,reDinucleotide_Element$CLine=="ADS_Adipose")

meanMethCoef_ADS<-apply(reDinucleotide_Element_ADS,1,FUN=methCoefMean)
meanMethCoef_ADS_Adipose<-apply(reDinucleotide_Element_ADS_Adipose,1,FUN=methCoefMean)

reDinucleotide_Element_ADS<-cbind(reDinucleotide_Element_ADS,meanMethCoef_ADS=as.numeric(as.vector(meanMethCoef_ADS)))
reDinucleotide_Element_ADS_Adipose<-cbind(reDinucleotide_Element_ADS_Adipose,meanMethCoef_ADS_Adipose=as.numeric(as.vector(meanMethCoef_ADS_Adipose)))

#################################
## With Informativeness filter ##
#################################
selected_reDinucleotide_Element_H1<-subset(reDinucleotide_Element_H1,
                                           reDinucleotide_Element_H1$meanMethCoef_H1!=-1 & 
                                             reDinucleotide_Element_H1$E_methCoef!=-1 & ((reDinucleotide_Element_H1$E_posInf/2)/reDinucleotide_Element_H1$E_nCG)>=0.25)
selected_reDinucleotide_Element_IMR90<-subset(reDinucleotide_Element_IMR90,
                                              reDinucleotide_Element_IMR90$meanMethCoef_IMR90!=-1 & 
                                                reDinucleotide_Element_IMR90$E_methCoef!=-1 & ((reDinucleotide_Element_IMR90$E_posInf/2)/reDinucleotide_Element_IMR90$E_nCG)>=0.25)
selected_reDinucleotide_Element_ADS<-subset(reDinucleotide_Element_ADS,
                                           reDinucleotide_Element_ADS$meanMethCoef_ADS!=-1 & 
                                             reDinucleotide_Element_ADS$E_methCoef!=-1 & ((reDinucleotide_Element_ADS$E_posInf/2)/reDinucleotide_Element_ADS$E_nCG)>=0.25)
selected_reDinucleotide_Element_ADS_Adipose<-subset(reDinucleotide_Element_ADS_Adipose,
                                              reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose!=-1 & 
                                                reDinucleotide_Element_ADS_Adipose$E_methCoef!=-1 & ((reDinucleotide_Element_ADS_Adipose$E_posInf/2)/reDinucleotide_Element_ADS_Adipose$E_nCG)>=0.25)
png(paste(resultsDIR,"figure25Filtered.png",sep=""),height=12,width=12,units="cm",res=300)
par(lwd=1.5)
par(cex.axis=0.8)

diffCpGiHpaIIH1<-as.data.frame(cbind(selected_reDinucleotide_Element_H1$E_methCoef-selected_reDinucleotide_Element_H1$meanMethCoef_H1,
                                     selected_reDinucleotide_Element_H1$E_methCoef))

diffCpGiHpaIIIMR90<-as.data.frame(cbind(selected_reDinucleotide_Element_IMR90$E_methCoef-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90,
                                        selected_reDinucleotide_Element_IMR90$E_methCoef))

diffCpGiHpaIIADS<-as.data.frame(cbind(selected_reDinucleotide_Element_ADS$E_methCoef-selected_reDinucleotide_Element_ADS$meanMethCoef_ADS,
                                     selected_reDinucleotide_Element_ADS$E_methCoef))

diffCpGiHpaIIADS_Adipose<-as.data.frame(cbind(selected_reDinucleotide_Element_ADS_Adipose$E_methCoef-selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose,
                                        selected_reDinucleotide_Element_ADS_Adipose$E_methCoef))

par(mfrow=c(2,2))

plot(density(diffCpGiHpaIIADS[,1]),type="n",xlim=c(-1,1),ylim=c(0,20),main="",
     xlab="",ylab="")
lines(density(diffCpGiHpaIIH1[,1]),col="darkblue",lwd=0.8)
lines(density(diffCpGiHpaIIIMR90[,1]),col="red",lwd=0.8)
lines(density(diffCpGiHpaIIADS[,1],width=0.05),col="green",lwd=0.8)
lines(density(diffCpGiHpaIIADS_Adipose[,1],width=0.05),col="purple",lwd=0.8)

plot(density(diffCpGiHpaIIADS[,1]),type="n",xlim=c(-1,1),ylim=c(0,20),main="",
     xlab="",ylab="")
lines(density(subset(diffCpGiHpaIIH1,diffCpGiHpaIIH1$V2>=0.75)[,1]),col="darkblue",lwd=0.8)
lines(density(subset(diffCpGiHpaIIIMR90,diffCpGiHpaIIIMR90$V2>=0.75)[,1]),col="red",lwd=0.8)
lines(density(subset(diffCpGiHpaIIADS,diffCpGiHpaIIADS$V2>=0.75)[,1]),col="green",lwd=0.8)
lines(density(subset(diffCpGiHpaIIADS_Adipose,diffCpGiHpaIIADS_Adipose$V2>=0.75)[,1]),col="purple",lwd=0.8)

plot(density(diffCpGiHpaIIADS[,1]),type="n",xlim=c(-1,1),ylim=c(0,20),main="",
     xlab="",ylab="")
lines(density(subset(diffCpGiHpaIIH1,diffCpGiHpaIIH1$V2<=0.25)[,1],width=0.05),col="darkblue",lwd=0.8)
lines(density(subset(diffCpGiHpaIIIMR90,diffCpGiHpaIIIMR90$V2<=0.25)[,1],width=0.05),col="red",lwd=0.8)
lines(density(subset(diffCpGiHpaIIADS,diffCpGiHpaIIADS$V2<=0.25)[,1],width=0.05),col="green",lwd=0.8)
lines(density(subset(diffCpGiHpaIIADS_Adipose,diffCpGiHpaIIADS_Adipose$V2<=0.25)[,1],width=0.05),col="purple",lwd=0.8)


plot(density(diffCpGiHpaIIADS[,1]),type="n",xlim=c(-1,1),ylim=c(0,20),main="",
     xlab="",ylab="")
lines(density(subset(diffCpGiHpaIIH1,diffCpGiHpaIIH1$V2>0.25 & diffCpGiHpaIIH1$V2<0.75)[,1]),col="darkblue",lwd=0.8)
lines(density(subset(diffCpGiHpaIIIMR90,diffCpGiHpaIIIMR90$V2>0.25 & diffCpGiHpaIIIMR90$V2<0.75)[,1]),col="red",lwd=0.8)
lines(density(subset(diffCpGiHpaIIADS,diffCpGiHpaIIADS$V2>0.25 & diffCpGiHpaIIADS$V2<0.75)[,1]),col="green",lwd=0.8)
lines(density(subset(diffCpGiHpaIIADS_Adipose,diffCpGiHpaIIADS_Adipose$V2>0.25 & diffCpGiHpaIIADS_Adipose$V2<0.75)[,1]),col="purple",lwd=0.8)

par(mfrow=c(1,1))

dev.off()

dataMatrix<-matrix(ncol=2,nrow=4)
colnames(dataMatrix)<-c("ADS","ADS_Adipose")
rownames(dataMatrix)<-c("All","Methylated","Unmethylated","Intermediate")
dataMatrix[1,1]<-length(diffCpGiHpaIIADS[,1])
dataMatrix[1,2]<-length(diffCpGiHpaIIADS_Adipose[,1])
dataMatrix[2,1]<-length(subset(diffCpGiHpaIIADS,diffCpGiHpaIIADS$V2>=0.75)[,1])
dataMatrix[2,2]<-length(subset(diffCpGiHpaIIADS_Adipose,diffCpGiHpaIIADS_Adipose$V2>=0.75)[,1])
dataMatrix[3,1]<-length(subset(diffCpGiHpaIIADS,diffCpGiHpaIIADS$V2<=0.25)[,1])
dataMatrix[3,2]<-length(subset(diffCpGiHpaIIADS_Adipose,diffCpGiHpaIIADS_Adipose$V2<=0.25)[,1])
dataMatrix[4,1]<-length(subset(diffCpGiHpaIIADS,diffCpGiHpaIIADS$V2>0.25 & diffCpGiHpaIIADS$V2<0.75)[,1])
dataMatrix[4,2]<-length(subset(diffCpGiHpaIIADS_Adipose,diffCpGiHpaIIADS_Adipose$V2>0.25 & diffCpGiHpaIIADS_Adipose$V2<0.75)[,1])

write.table(dataMatrix,file=paste(resultsDIR,"figure25FilteredValues.txt",sep=""),sep="\t")

####################################
## Without Informativeness filter ##
####################################
selected_reDinucleotide_Element_H1<-subset(reDinucleotide_Element_H1,
                                           reDinucleotide_Element_H1$meanMethCoef_H1!=-1 & 
                                             reDinucleotide_Element_H1$E_methCoef!=-1)
selected_reDinucleotide_Element_IMR90<-subset(reDinucleotide_Element_IMR90,
                                              reDinucleotide_Element_IMR90$meanMethCoef_IMR90!=-1 & 
                                                reDinucleotide_Element_IMR90$E_methCoef!=-1)

selected_reDinucleotide_Element_ADS<-subset(reDinucleotide_Element_ADS,
                                           reDinucleotide_Element_ADS$meanMethCoef_ADS!=-1 & 
                                             reDinucleotide_Element_ADS$E_methCoef!=-1)
selected_reDinucleotide_Element_ADS_Adipose<-subset(reDinucleotide_Element_ADS_Adipose,
                                              reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose!=-1 & 
                                                reDinucleotide_Element_ADS_Adipose$E_methCoef!=-1)

#######################

png(paste(resultsDIR,"figure25All.png",sep=""),height=12,width=12,units="cm",res=300)
par(lwd=1.5)
par(cex.axis=0.8)

diffCpGiHpaIIH1<-as.data.frame(cbind(selected_reDinucleotide_Element_H1$E_methCoef-selected_reDinucleotide_Element_H1$meanMethCoef_H1,
                                     selected_reDinucleotide_Element_H1$E_methCoef))

diffCpGiHpaIIIMR90<-as.data.frame(cbind(selected_reDinucleotide_Element_IMR90$E_methCoef-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90,
                                        selected_reDinucleotide_Element_IMR90$E_methCoef))

diffCpGiHpaIIADS<-as.data.frame(cbind(selected_reDinucleotide_Element_ADS$E_methCoef-selected_reDinucleotide_Element_ADS$meanMethCoef_ADS,
                                     selected_reDinucleotide_Element_ADS$E_methCoef))

diffCpGiHpaIIADS_Adipose<-as.data.frame(cbind(selected_reDinucleotide_Element_ADS_Adipose$E_methCoef-selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose,
                                        selected_reDinucleotide_Element_ADS_Adipose$E_methCoef))

par(mfrow=c(2,2))

plot(density(diffCpGiHpaIIADS[,1]),type="n",xlim=c(-1,1),ylim=c(0,20),main="",
     xlab="",ylab="")
lines(density(diffCpGiHpaIIH1[,1],width=0.05),col="darkblue",lwd=0.8)
lines(density(diffCpGiHpaIIIMR90[,1],width=0.05),col="red",lwd=0.8)
lines(density(diffCpGiHpaIIADS[,1],width=0.05),col="green",lwd=0.8)
lines(density(diffCpGiHpaIIADS_Adipose[,1],width=0.05),col="purple",lwd=0.8)


plot(density(diffCpGiHpaIIADS[,1]),type="n",xlim=c(-1,1),ylim=c(0,20),main="",
     xlab="",ylab="")
lines(density(subset(diffCpGiHpaIIH1,diffCpGiHpaIIH1$V2>=0.75)[,1]),col="darkblue",lwd=0.8)
lines(density(subset(diffCpGiHpaIIIMR90,diffCpGiHpaIIIMR90$V2>=0.75)[,1]),col="red",lwd=0.8)
lines(density(subset(diffCpGiHpaIIADS,diffCpGiHpaIIADS$V2>=0.75)[,1]),col="green",lwd=0.8)
lines(density(subset(diffCpGiHpaIIADS_Adipose,diffCpGiHpaIIADS_Adipose$V2>=0.75)[,1]),col="purple",lwd=0.8)



plot(density(diffCpGiHpaIIADS[,1]),type="n",xlim=c(-1,1),ylim=c(0,20),main="",
     xlab="",ylab="")
lines(density(subset(diffCpGiHpaIIH1,diffCpGiHpaIIH1$V2<=0.25)[,1],width=0.05),col="darkblue",lwd=0.8)
lines(density(subset(diffCpGiHpaIIIMR90,diffCpGiHpaIIIMR90$V2<=0.25)[,1],width=0.05),col="red",lwd=0.8)
lines(density(subset(diffCpGiHpaIIADS,diffCpGiHpaIIADS$V2<=0.25)[,1],width=0.05),col="green",lwd=0.8)
lines(density(subset(diffCpGiHpaIIADS_Adipose,diffCpGiHpaIIADS_Adipose$V2<=0.25)[,1],width=0.05),col="purple",lwd=0.8)


plot(density(diffCpGiHpaIIADS[,1]),type="n",xlim=c(-1,1),ylim=c(0,20),main="",
     xlab="",ylab="")
lines(density(subset(diffCpGiHpaIIH1,diffCpGiHpaIIH1$V2>0.25 & diffCpGiHpaIIH1$V2<0.75)[,1]),col="darkblue",lwd=0.8)
lines(density(subset(diffCpGiHpaIIIMR90,diffCpGiHpaIIIMR90$V2>0.25 & diffCpGiHpaIIIMR90$V2<0.75)[,1]),col="red",lwd=0.8)
lines(density(subset(diffCpGiHpaIIADS,diffCpGiHpaIIADS$V2>0.25 & diffCpGiHpaIIADS$V2<0.75)[,1]),col="green",lwd=0.8)
lines(density(subset(diffCpGiHpaIIADS_Adipose,diffCpGiHpaIIADS_Adipose$V2>0.25 & diffCpGiHpaIIADS_Adipose$V2<0.75)[,1]),col="purple",lwd=0.8)

par(mfrow=c(1,1))

dev.off()

dataMatrix<-matrix(ncol=2,nrow=4)
colnames(dataMatrix)<-c("ADS","ADS_Adipose")
rownames(dataMatrix)<-c("All","Methylated","Unmethylated","Intermediate")
dataMatrix[1,1]<-length(diffCpGiHpaIIADS[,1])
dataMatrix[1,2]<-length(diffCpGiHpaIIADS_Adipose[,1])
dataMatrix[2,1]<-length(subset(diffCpGiHpaIIADS,diffCpGiHpaIIADS$V2>=0.75)[,1])
dataMatrix[2,2]<-length(subset(diffCpGiHpaIIADS_Adipose,diffCpGiHpaIIADS_Adipose$V2>=0.75)[,1])
dataMatrix[3,1]<-length(subset(diffCpGiHpaIIADS,diffCpGiHpaIIADS$V2<=0.25)[,1])
dataMatrix[3,2]<-length(subset(diffCpGiHpaIIADS_Adipose,diffCpGiHpaIIADS_Adipose$V2<=0.25)[,1])
dataMatrix[4,1]<-length(subset(diffCpGiHpaIIADS,diffCpGiHpaIIADS$V2>0.25 & diffCpGiHpaIIADS$V2<0.75)[,1])
dataMatrix[4,2]<-length(subset(diffCpGiHpaIIADS_Adipose,diffCpGiHpaIIADS_Adipose$V2>0.25 & diffCpGiHpaIIADS_Adipose$V2<0.75)[,1])

write.table(dataMatrix,file=paste(resultsDIR,"figure25ValuesAll.txt",sep=""),sep="\t")