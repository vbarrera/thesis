##################
### Figure 27 ###
##################

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

reDinucleotide_Element_ADS<-subset(reDinucleotide_Element,reDinucleotide_Element$CLine=="ADS")
reDinucleotide_Element_ADS_Adipose<-subset(reDinucleotide_Element,reDinucleotide_Element$CLine=="ADS_Adipose")
reDinucleotide_Element_ADS_IPSC<-subset(reDinucleotide_Element,reDinucleotide_Element$CLine=="ADS_IPSC")

meanMethCoef_ADS<-apply(reDinucleotide_Element_ADS,1,FUN=methCoefMean)
meanMethCoef_ADS_Adipose<-apply(reDinucleotide_Element_ADS_Adipose,1,FUN=methCoefMean)
meanMethCoef_ADS_IPSC<-apply(reDinucleotide_Element_ADS_IPSC,1,FUN=methCoefMean)

reDinucleotide_Element_ADS<-cbind(reDinucleotide_Element_ADS,meanMethCoef_ADS=as.numeric(as.vector(meanMethCoef_ADS)))
reDinucleotide_Element_ADS_Adipose<-cbind(reDinucleotide_Element_ADS_Adipose,meanMethCoef_ADS_Adipose=as.numeric(as.vector(meanMethCoef_ADS_Adipose)))
reDinucleotide_Element_ADS_IPSC<-cbind(reDinucleotide_Element_ADS_IPSC,meanMethCoef_ADS_IPSC=as.numeric(as.vector(meanMethCoef_ADS_IPSC)))

#################################
## With Informativeness filter ##
#################################

selected_reDinucleotide_Element_ADS<-subset(reDinucleotide_Element_ADS,
                                           reDinucleotide_Element_ADS$meanMethCoef_ADS!=-1 & 
                                             reDinucleotide_Element_ADS$E_methCoef!=-1 & ((reDinucleotide_Element_ADS$E_posInf/2)/reDinucleotide_Element_ADS$E_nCG)>=0.25)
selected_reDinucleotide_Element_ADS_Adipose<-subset(reDinucleotide_Element_ADS_Adipose,
                                              reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose!=-1 & 
                                                reDinucleotide_Element_ADS_Adipose$E_methCoef!=-1 & ((reDinucleotide_Element_ADS_Adipose$E_posInf/2)/reDinucleotide_Element_ADS_Adipose$E_nCG)>=0.25)
selected_reDinucleotide_Element_ADS_IPSC<-subset(reDinucleotide_Element_ADS_IPSC,
                                                    reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC!=-1 & 
                                                      reDinucleotide_Element_ADS_IPSC$E_methCoef!=-1 & ((reDinucleotide_Element_ADS_IPSC$E_posInf/2)/reDinucleotide_Element_ADS_IPSC$E_nCG)>=0.25)


# write.table(cbind(selected_reDinucleotide_Element_ADS$id_In_RE,selected_reDinucleotide_Element_ADS$meanMethCoef_ADS,selected_reDinucleotide_Element_ADS$id_In_Type,selected_reDinucleotide_Element_ADS$E_methCoef),file=paste(resultsDIR,"reducedDataADS_Filtered.txt",sep=""),sep="\t",col.names=FALSE,row.names=F)
# write.table(cbind(selected_reDinucleotide_Element_ADS_Adipose$id_In_RE,selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose,selected_reDinucleotide_Element_ADS_Adipose$id_In_Type,selected_reDinucleotide_Element_ADS_Adipose$E_methCoef),file=paste(resultsDIR,"reducedDataADS_Adipose_Filtered.txt",sep=""),sep="\t",col.names=FALSE,row.names=F)
# write.table(cbind(selected_reDinucleotide_Element_ADS_IPSC$id_In_RE,selected_reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC,selected_reDinucleotide_Element_ADS_IPSC$id_In_Type,selected_reDinucleotide_Element_ADS_IPSC$E_methCoef),file=paste(resultsDIR,"reducedDataADS_IPSC_Filtered.txt",sep=""),sep="\t",col.names=FALSE,row.names=F)

##########################
# Distance to extreme
t<-cbind(x=selected_reDinucleotide_Element_ADS$RE_position-selected_reDinucleotide_Element_ADS$E_chromStart,y=selected_reDinucleotide_Element_ADS$E_chromEnd-selected_reDinucleotide_Element_ADS$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_ADS<-cbind(selected_reDinucleotide_Element_ADS,dist2Ext=m)

t<-cbind(x=selected_reDinucleotide_Element_ADS_Adipose$RE_position-selected_reDinucleotide_Element_ADS_Adipose$E_chromStart,y=selected_reDinucleotide_Element_ADS_Adipose$E_chromEnd-selected_reDinucleotide_Element_ADS_Adipose$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_ADS_Adipose<-cbind(selected_reDinucleotide_Element_ADS_Adipose,dist2Ext=m)

t<-cbind(x=selected_reDinucleotide_Element_ADS_IPSC$RE_position-selected_reDinucleotide_Element_ADS_IPSC$E_chromStart,y=selected_reDinucleotide_Element_ADS_IPSC$E_chromEnd-selected_reDinucleotide_Element_ADS_IPSC$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_ADS_IPSC<-cbind(selected_reDinucleotide_Element_ADS_IPSC,dist2Ext=m)


#Obtain O/E ratio of the element

require(RMySQL)
con<-dbConnect(MySQL(),dbname="SolexaDB")
query<-"select id_cpgi_pk,obsExp FROM CpGisland"
rs1<-dbSendQuery(con,query)
OEdata<-fetch(rs1,n=-1) 
colnames(OEdata)=c("id_In_Type","E_OE")

selected_reDinucleotide_Element_ADS<-merge(selected_reDinucleotide_Element_ADS,OEdata,group="id_In_Type")
selected_reDinucleotide_Element_ADS_Adipose<-merge(selected_reDinucleotide_Element_ADS_Adipose,OEdata,group="id_In_Type")
selected_reDinucleotide_Element_ADS_IPSC<-merge(selected_reDinucleotide_Element_ADS_IPSC,OEdata,group="id_In_Type")

selected_reDinucleotide_Element_ADS<-selected_reDinucleotide_Element_ADS[,c(2:7,1,8:17)]
selected_reDinucleotide_Element_ADS_Adipose<-selected_reDinucleotide_Element_ADS_Adipose[,c(2:7,1,8:17)]
selected_reDinucleotide_Element_ADS_IPSC<-selected_reDinucleotide_Element_ADS_IPSC[,c(2:7,1,8:17)]

#Obtain length of the element

selected_reDinucleotide_Element_ADS<-cbind(selected_reDinucleotide_Element_ADS,
                                          cgiLength=selected_reDinucleotide_Element_ADS$E_chromEnd-
                                            selected_reDinucleotide_Element_ADS$E_chromStart)

selected_reDinucleotide_Element_ADS_Adipose<-cbind(selected_reDinucleotide_Element_ADS_Adipose,
                                             cgiLength=selected_reDinucleotide_Element_ADS_Adipose$E_chromEnd-
                                               selected_reDinucleotide_Element_ADS_Adipose$E_chromStart)
selected_reDinucleotide_Element_ADS_IPSC<-cbind(selected_reDinucleotide_Element_ADS_IPSC,
                                                   cgiLength=selected_reDinucleotide_Element_ADS_IPSC$E_chromEnd-
                                                     selected_reDinucleotide_Element_ADS_IPSC$E_chromStart)

############
#ADS

#lmADSFiltered<-lm(selected_reDinucleotide_Element_ADS$E_methCoef~selected_reDinucleotide_Element_ADS$meanMethCoef_ADS)
#summary(lmADSFiltered)
selected_reDinucleotide_Element_ADSJitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_ADS$E_methCoef,factor=25),meanMethCoef_ADS=jitter(selected_reDinucleotide_Element_ADS$meanMethCoef_ADS,factor=25)))

  
incorrectPoints_ADS_plus<-subset(selected_reDinucleotide_Element_ADS,selected_reDinucleotide_Element_ADS$E_methCoef-selected_reDinucleotide_Element_ADS$meanMethCoef_ADS>0.25)
#x_ADS_plus<-subset(selected_reDinucleotide_Element_ADSJitter,selected_reDinucleotide_Element_ADSJitter$E_methCoef-selected_reDinucleotide_Element_ADSJitter$meanMethCoef_ADS>0.25)
dim(incorrectPoints_ADS_plus)
x_ADS_plus<-subset(selected_reDinucleotide_Element_ADSJitter,selected_reDinucleotide_Element_ADSJitter$E_methCoef>=0.75)

incorrectPoints_ADS_minus<-subset(selected_reDinucleotide_Element_ADS,selected_reDinucleotide_Element_ADS$E_methCoef-selected_reDinucleotide_Element_ADS$meanMethCoef_ADS<(-0.25))
#x_ADS_minus<-subset(selected_reDinucleotide_Element_ADSJitter,selected_reDinucleotide_Element_ADSJitter$E_methCoef-selected_reDinucleotide_Element_ADSJitter$meanMethCoef_ADS<(-0.25))
dim(incorrectPoints_ADS_minus)
x_ADS_minus<-subset(selected_reDinucleotide_Element_ADSJitter,selected_reDinucleotide_Element_ADSJitter$E_methCoef<=0.25)

correctPoints_ADS<-subset(selected_reDinucleotide_Element_ADS,abs(selected_reDinucleotide_Element_ADS$E_methCoef-selected_reDinucleotide_Element_ADS$meanMethCoef_ADS)<=0.25)
# y_ADS<-subset(selected_reDinucleotide_Element_ADSJitter,abs(selected_reDinucleotide_Element_ADSJitter$E_methCoef-selected_reDinucleotide_Element_ADSJitter$meanMethCoef_ADS)<=0.25)
dim(correctPoints_ADS)
y_ADS<-subset(selected_reDinucleotide_Element_ADSJitter,selected_reDinucleotide_Element_ADSJitter$E_methCoef>0.25 & selected_reDinucleotide_Element_ADSJitter<0.75)

if (numbered=='TRUE'){
  png(paste(resultsDIR,"figure27_ADSFilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
} else {
  png(paste(resultsDIR,"figure27_ADSFilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_ADS$meanMethCoef_ADS, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_ADS$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_ADS$meanMethCoef_ADS,selected_reDinucleotide_Element_ADS$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n",xlim=c(0,1),ylim=c(0,1))

## Plot the points with some noise (jitter)
points(x_ADS_plus$meanMethCoef_ADS,x_ADS_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_ADS_minus$meanMethCoef_ADS,x_ADS_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_ADS$meanMethCoef_ADS,y_ADS$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# Plot the points without noise
#points(incorrectPoints_ADS_plus$meanMethCoef_ADS,incorrectPoints_ADS_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_ADS_minus$meanMethCoef_ADS,incorrectPoints_ADS_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_ADS$meanMethCoef_ADS,correctPoints_ADS$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

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

# lmADS_AdiposeFiltered<-lm(selected_reDinucleotide_Element_ADS_Adipose$E_methCoef~selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose)
# summary(lmADS_AdiposeFiltered)

selected_reDinucleotide_Element_ADS_AdiposeJitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_ADS_Adipose$E_methCoef,factor=25),meanMethCoef_ADS_Adipose=jitter(selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose,factor=25)))

incorrectPoints_ADS_Adipose_plus<-subset(selected_reDinucleotide_Element_ADS_Adipose,selected_reDinucleotide_Element_ADS_Adipose$E_methCoef-selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose>0.25)
#x_ADS_Adipose_plus<-subset(selected_reDinucleotide_Element_ADS_AdiposeJitter,selected_reDinucleotide_Element_ADS_AdiposeJitter$E_methCoef-selected_reDinucleotide_Element_ADS_AdiposeJitter$meanMethCoef_ADS_Adipose>0.25)
dim(incorrectPoints_ADS_Adipose_plus)
x_ADS_Adipose_plus<-subset(selected_reDinucleotide_Element_ADS_AdiposeJitter,selected_reDinucleotide_Element_ADS_AdiposeJitter$E_methCoef>=0.75)

incorrectPoints_ADS_Adipose_minus<-subset(selected_reDinucleotide_Element_ADS_Adipose,selected_reDinucleotide_Element_ADS_Adipose$E_methCoef-selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose<(-0.25))
#x_ADS_Adipose_minus<-subset(selected_reDinucleotide_Element_ADS_AdiposeJitter,selected_reDinucleotide_Element_ADS_AdiposeJitter$E_methCoef-selected_reDinucleotide_Element_ADS_AdiposeJitter$meanMethCoef_ADS_Adipose<(-0.25))
dim(incorrectPoints_ADS_Adipose_minus)
x_ADS_Adipose_minus<-subset(selected_reDinucleotide_Element_ADS_AdiposeJitter,selected_reDinucleotide_Element_ADS_AdiposeJitter$E_methCoef<=0.25)

correctPoints_ADS_Adipose<-subset(selected_reDinucleotide_Element_ADS_Adipose,abs(selected_reDinucleotide_Element_ADS_Adipose$E_methCoef-selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose)<=0.25)
#y_ADS_Adipose<-subset(selected_reDinucleotide_Element_ADS_AdiposeJitter,abs(selected_reDinucleotide_Element_ADS_AdiposeJitter$E_methCoef-selected_reDinucleotide_Element_ADS_AdiposeJitter$meanMethCoef_ADS_Adipose)<=0.25)
dim(correctPoints_ADS_Adipose)
y_ADS_Adipose<-subset(selected_reDinucleotide_Element_ADS_AdiposeJitter,selected_reDinucleotide_Element_ADS_AdiposeJitter$E_methCoef>0.25 & selected_reDinucleotide_Element_ADS_AdiposeJitter<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure27_ADS_AdiposeFilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
  png(paste(resultsDIR,"figure27_ADS_AdiposeFilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_ADS_Adipose$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose,selected_reDinucleotide_Element_ADS_Adipose$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n",xlim=c(0,1),ylim=c(0,1))

## Plot the points with some noise (jitter)
points(x_ADS_Adipose_plus$meanMethCoef_ADS_Adipose,x_ADS_Adipose_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_ADS_Adipose_minus$meanMethCoef_ADS_Adipose,x_ADS_Adipose_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_ADS_Adipose$meanMethCoef_ADS_Adipose,y_ADS_Adipose$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# Plot the point without noise
#points(incorrectPoints_ADS_Adipose_plus$meanMethCoef_ADS_Adipose,incorrectPoints_ADS_Adipose_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_ADS_Adipose_minus$meanMethCoef_ADS_Adipose,incorrectPoints_ADS_Adipose_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_ADS_Adipose$meanMethCoef_ADS_Adipose,correctPoints_ADS_Adipose$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if(numbered=='TRUE'){
# Adding the number of points to the plot
text(0.2,0.8,dim(incorrectPoints_ADS_Adipose_plus)[1],cex=1.2,col="red")
text(0.8,0.2,dim(incorrectPoints_ADS_Adipose_minus)[1],cex=1.2,col="red")
text(0.5,0.5,dim(correctPoints_ADS_Adipose)[1],cex=1.2,col="red")
}else{}
abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

groupCPD_ADS_Adipose<-subset(correctPoints_ADS_Adipose,correctPoints_ADS_Adipose$E_methCoef<=0.25)
groupCPM_ADS_Adipose<-subset(correctPoints_ADS_Adipose,correctPoints_ADS_Adipose$E_methCoef>0.25 & correctPoints_ADS_Adipose$E_methCoef<0.75)
groupCPU_ADS_Adipose<-subset(correctPoints_ADS_Adipose,correctPoints_ADS_Adipose$E_methCoef>=0.75)

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
# lmADS_IPSCFiltered<-lm(selected_reDinucleotide_Element_ADS_IPSC$E_methCoef~selected_reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC)
# summary(lmADS_IPSCFiltered)

selected_reDinucleotide_Element_ADS_IPSCJitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_ADS_IPSC$E_methCoef,factor=25),meanMethCoef_ADS_IPSC=jitter(selected_reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC,factor=25)))

incorrectPoints_ADS_IPSC_plus<-subset(selected_reDinucleotide_Element_ADS_IPSC,selected_reDinucleotide_Element_ADS_IPSC$E_methCoef-selected_reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC>0.25)
#x_ADS_IPSC_plus<-subset(selected_reDinucleotide_Element_ADS_IPSCJitter,selected_reDinucleotide_Element_ADS_IPSCJitter$E_methCoef-selected_reDinucleotide_Element_ADS_IPSCJitter$meanMethCoef_ADS_IPSC>0.25)
dim(incorrectPoints_ADS_IPSC_plus)
x_ADS_IPSC_plus<-subset(selected_reDinucleotide_Element_ADS_IPSCJitter,selected_reDinucleotide_Element_ADS_IPSCJitter$E_methCoef>=0.75)

incorrectPoints_ADS_IPSC_minus<-subset(selected_reDinucleotide_Element_ADS_IPSC,selected_reDinucleotide_Element_ADS_IPSC$E_methCoef-selected_reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC<(-0.25))
#x_ADS_IPSC_minus<-subset(selected_reDinucleotide_Element_ADS_IPSCJitter,selected_reDinucleotide_Element_ADS_IPSCJitter$E_methCoef-selected_reDinucleotide_Element_ADS_IPSCJitter$meanMethCoef_ADS_IPSC<(-0.25))
dim(incorrectPoints_ADS_IPSC_minus)
x_ADS_IPSC_minus<-subset(selected_reDinucleotide_Element_ADS_IPSCJitter,selected_reDinucleotide_Element_ADS_IPSCJitter$E_methCoef<=0.25)

correctPoints_ADS_IPSC<-subset(selected_reDinucleotide_Element_ADS_IPSC,abs(selected_reDinucleotide_Element_ADS_IPSC$E_methCoef-selected_reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC)<=0.25)
#y_ADS_IPSC<-subset(selected_reDinucleotide_Element_ADS_IPSCJitter,abs(selected_reDinucleotide_Element_ADS_IPSCJitter$E_methCoef-selected_reDinucleotide_Element_ADS_IPSCJitter$meanMethCoef_ADS_IPSC)<=0.25)
dim(correctPoints_ADS_IPSC)
y_ADS_IPSC<-subset(selected_reDinucleotide_Element_ADS_IPSCJitter,selected_reDinucleotide_Element_ADS_IPSCJitter$E_methCoef>0.25 & selected_reDinucleotide_Element_ADS_IPSCJitter<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure28_ADS_IPSCFilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
  png(paste(resultsDIR,"figure28_ADS_IPSCFilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_ADS_IPSC$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC,selected_reDinucleotide_Element_ADS_IPSC$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n",xlim=c(0,1),ylim=c(0,1))

## Plot the points with some noise (jitter)
points(x_ADS_IPSC_plus$meanMethCoef_ADS_IPSC,x_ADS_IPSC_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_ADS_IPSC_minus$meanMethCoef_ADS_IPSC,x_ADS_IPSC_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_ADS_IPSC$meanMethCoef_ADS_IPSC,y_ADS_IPSC$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# Plot the point without noise
#points(incorrectPoints_ADS_IPSC_plus$meanMethCoef_ADS_IPSC,incorrectPoints_ADS_IPSC_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_ADS_IPSC_minus$meanMethCoef_ADS_IPSC,incorrectPoints_ADS_IPSC_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_ADS_IPSC$meanMethCoef_ADS_IPSC,correctPoints_ADS_IPSC$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if(numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_ADS_IPSC_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_ADS_IPSC_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_ADS_IPSC)[1],cex=1.2,col="red")
}else{}
abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

groupCPD_ADS_IPSC<-subset(correctPoints_ADS_IPSC,correctPoints_ADS_IPSC$E_methCoef<=0.25)
groupCPM_ADS_IPSC<-subset(correctPoints_ADS_IPSC,correctPoints_ADS_IPSC$E_methCoef>0.25 & correctPoints_ADS_IPSC$E_methCoef<0.75)
groupCPU_ADS_IPSC<-subset(correctPoints_ADS_IPSC,correctPoints_ADS_IPSC$E_methCoef>=0.75)

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

selected_reDinucleotide_Element_ADS<-subset(reDinucleotide_Element_ADS,
                                           reDinucleotide_Element_ADS$meanMethCoef_ADS!=-1 & 
                                             reDinucleotide_Element_ADS$E_methCoef!=-1)
selected_reDinucleotide_Element_ADS_Adipose<-subset(reDinucleotide_Element_ADS_Adipose,
                                              reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose!=-1 & 
                                                reDinucleotide_Element_ADS_Adipose$E_methCoef!=-1)
selected_reDinucleotide_Element_ADS_IPSC<-subset(reDinucleotide_Element_ADS_IPSC,
                                                    reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC!=-1 & 
                                                      reDinucleotide_Element_ADS_IPSC$E_methCoef!=-1)

# write.table(cbind(selected_reDinucleotide_Element_ADS$id_In_RE,selected_reDinucleotide_Element_ADS$meanMethCoef_ADS,selected_reDinucleotide_Element_ADS$id_In_Type,selected_reDinucleotide_Element_ADS$E_methCoef),file=paste(resultsDIR,"reducedDataADS_All.txt",sep=""),sep="\t",col.names=FALSE,row.names=F)
# write.table(cbind(selected_reDinucleotide_Element_ADS_Adipose$id_In_RE,selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose,selected_reDinucleotide_Element_ADS_Adipose$id_In_Type,selected_reDinucleotide_Element_ADS_Adipose$E_methCoef),file=paste(resultsDIR,"reducedDataADS_Adipose_All.txt",sep=""),sep="\t",col.names=FALSE,row.names=F)
# write.table(cbind(selected_reDinucleotide_Element_ADS_IPSC$id_In_RE,selected_reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC,selected_reDinucleotide_Element_ADS_IPSC$id_In_Type,selected_reDinucleotide_Element_ADS_IPSC$E_methCoef),file=paste(resultsDIR,"reducedDataADS_IPSC_All.txt",sep=""),sep="\t",col.names=FALSE,row.names=F)



#######################
# Distance to extreme
t<-cbind(x=selected_reDinucleotide_Element_ADS$RE_position-selected_reDinucleotide_Element_ADS$E_chromStart,y=selected_reDinucleotide_Element_ADS$E_chromEnd-selected_reDinucleotide_Element_ADS$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_ADS<-cbind(selected_reDinucleotide_Element_ADS,dist2Ext=m)

t<-cbind(x=selected_reDinucleotide_Element_ADS_Adipose$RE_position-selected_reDinucleotide_Element_ADS_Adipose$E_chromStart,y=selected_reDinucleotide_Element_ADS_Adipose$E_chromEnd-selected_reDinucleotide_Element_ADS_Adipose$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_ADS_Adipose<-cbind(selected_reDinucleotide_Element_ADS_Adipose,dist2Ext=m)


t<-cbind(x=selected_reDinucleotide_Element_ADS_IPSC$RE_position-selected_reDinucleotide_Element_ADS_IPSC$E_chromStart,y=selected_reDinucleotide_Element_ADS_IPSC$E_chromEnd-selected_reDinucleotide_Element_ADS_IPSC$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_ADS_IPSC<-cbind(selected_reDinucleotide_Element_ADS_IPSC,dist2Ext=m)

#Obtain O/E ratio of the element

require(RMySQL)
con<-dbConnect(MySQL(),dbname="SolexaDB")
query<-"select id_cpgi_pk,obsExp FROM CpGisland"
rs1<-dbSendQuery(con,query)
OEdata<-fetch(rs1,n=-1) 
colnames(OEdata)=c("id_In_Type","E_OE")

selected_reDinucleotide_Element_ADS<-merge(selected_reDinucleotide_Element_ADS,OEdata,group="id_In_Type")
selected_reDinucleotide_Element_ADS_Adipose<-merge(selected_reDinucleotide_Element_ADS_Adipose,OEdata,group="id_In_Type")
selected_reDinucleotide_Element_ADS_IPSC<-merge(selected_reDinucleotide_Element_ADS_IPSC,OEdata,group="id_In_Type")


selected_reDinucleotide_Element_ADS<-selected_reDinucleotide_Element_ADS[,c(2:7,1,8:17)]
selected_reDinucleotide_Element_ADS_Adipose<-selected_reDinucleotide_Element_ADS_Adipose[,c(2:7,1,8:17)]
selected_reDinucleotide_Element_ADS_IPSC<-selected_reDinucleotide_Element_ADS_IPSC[,c(2:7,1,8:17)]

#Obtain length of the element

selected_reDinucleotide_Element_ADS<-cbind(selected_reDinucleotide_Element_ADS,
                                          cgiLength=selected_reDinucleotide_Element_ADS$E_chromEnd-
                                            selected_reDinucleotide_Element_ADS$E_chromStart)

selected_reDinucleotide_Element_ADS_Adipose<-cbind(selected_reDinucleotide_Element_ADS_Adipose,
                                             cgiLength=selected_reDinucleotide_Element_ADS_Adipose$E_chromEnd-
                                               selected_reDinucleotide_Element_ADS_Adipose$E_chromStart)
selected_reDinucleotide_Element_ADS_IPSC<-cbind(selected_reDinucleotide_Element_ADS_IPSC,
                                                   cgiLength=selected_reDinucleotide_Element_ADS_IPSC$E_chromEnd-
                                                     selected_reDinucleotide_Element_ADS_IPSC$E_chromStart)

############
#ADS
# lmADSAll<-lm(selected_reDinucleotide_Element_ADS$E_methCoef~selected_reDinucleotide_Element_ADS$meanMethCoef_ADS)
# summary(lmADSAll)
selected_reDinucleotide_Element_ADSJitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_ADS$E_methCoef,factor=25),meanMethCoef_ADS=jitter(selected_reDinucleotide_Element_ADS$meanMethCoef_ADS,factor=25)))

incorrectPoints_ADS_plus<-subset(selected_reDinucleotide_Element_ADS,selected_reDinucleotide_Element_ADS$E_methCoef-selected_reDinucleotide_Element_ADS$meanMethCoef_ADS>0.25)
#x_ADS_plus<-subset(selected_reDinucleotide_Element_ADSJitter,selected_reDinucleotide_Element_ADSJitter$E_methCoef-selected_reDinucleotide_Element_ADSJitter$meanMethCoef_ADS>0.25)
dim(incorrectPoints_ADS_plus)
x_ADS_plus<-subset(selected_reDinucleotide_Element_ADSJitter,selected_reDinucleotide_Element_ADSJitter$E_methCoef>=0.75)

incorrectPoints_ADS_minus<-subset(selected_reDinucleotide_Element_ADS,selected_reDinucleotide_Element_ADS$E_methCoef-selected_reDinucleotide_Element_ADS$meanMethCoef_ADS<(-0.25))
#x_ADS_minus<-subset(selected_reDinucleotide_Element_ADSJitter,selected_reDinucleotide_Element_ADSJitter$E_methCoef-selected_reDinucleotide_Element_ADSJitter$meanMethCoef_ADS<(-0.25))
dim(incorrectPoints_ADS_minus)
x_ADS_minus<-subset(selected_reDinucleotide_Element_ADSJitter,selected_reDinucleotide_Element_ADSJitter$E_methCoef<=0.25)

correctPoints_ADS<-subset(selected_reDinucleotide_Element_ADS,abs(selected_reDinucleotide_Element_ADS$E_methCoef-selected_reDinucleotide_Element_ADS$meanMethCoef_ADS)<=0.25)
# y_ADS<-subset(selected_reDinucleotide_Element_ADSJitter,abs(selected_reDinucleotide_Element_ADSJitter$E_methCoef-selected_reDinucleotide_Element_ADSJitter$meanMethCoef_ADS)<=0.25)
dim(correctPoints_ADS)
y_ADS<-subset(selected_reDinucleotide_Element_ADSJitter,selected_reDinucleotide_Element_ADSJitter$E_methCoef>0.25 & selected_reDinucleotide_Element_ADSJitter<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure27_ADSAllBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure27_ADSAllBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_ADS$meanMethCoef_ADS, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_ADS$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_ADS$meanMethCoef_ADS,selected_reDinucleotide_Element_ADS$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n",xlim=c(0,1),ylim=c(0,1))
points(x_ADS_plus$meanMethCoef_ADS,x_ADS_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_ADS_minus$meanMethCoef_ADS,x_ADS_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_ADS$meanMethCoef_ADS,y_ADS$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

#points(incorrectPoints_ADS_plus$meanMethCoef_ADS,incorrectPoints_ADS_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_ADS_minus$meanMethCoef_ADS,incorrectPoints_ADS_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_ADS$meanMethCoef_ADS,correctPoints_ADS$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if (numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_ADS_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_ADS_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_ADS)[1],cex=1.2,col="red")
}else{}

abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

groupCPD_ADS<-subset(correctPoints_ADS,correctPoints_ADS$E_methCoef<=0.25)
groupCPM_ADS<-subset(correctPoints_ADS,correctPoints_ADS$E_methCoef>0.25 & correctPoints_ADS$E_methCoef<0.75)
groupCPU_ADS<-subset(correctPoints_ADS,correctPoints_ADS$E_methCoef>=0.75)



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

# lmADS_AdiposeAll<-lm(selected_reDinucleotide_Element_ADS_Adipose$E_methCoef~selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose)
# summary(lmADS_AdiposeAll)

selected_reDinucleotide_Element_ADS_AdiposeJitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_ADS_Adipose$E_methCoef,factor=25),meanMethCoef_ADS_Adipose=jitter(selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose,factor=25)))

incorrectPoints_ADS_Adipose_plus<-subset(selected_reDinucleotide_Element_ADS_Adipose,selected_reDinucleotide_Element_ADS_Adipose$E_methCoef-selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose>0.25)
#x_ADS_Adipose_plus<-subset(selected_reDinucleotide_Element_ADS_AdiposeJitter,selected_reDinucleotide_Element_ADS_AdiposeJitter$E_methCoef-selected_reDinucleotide_Element_ADS_AdiposeJitter$meanMethCoef_ADS_Adipose>0.25)
dim(incorrectPoints_ADS_Adipose_plus)
x_ADS_Adipose_plus<-subset(selected_reDinucleotide_Element_ADS_AdiposeJitter,selected_reDinucleotide_Element_ADS_AdiposeJitter$E_methCoef>=0.75)

incorrectPoints_ADS_Adipose_minus<-subset(selected_reDinucleotide_Element_ADS_Adipose,selected_reDinucleotide_Element_ADS_Adipose$E_methCoef-selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose<(-0.25))
#x_ADS_Adipose_minus<-subset(selected_reDinucleotide_Element_ADS_AdiposeJitter,selected_reDinucleotide_Element_ADS_AdiposeJitter$E_methCoef-selected_reDinucleotide_Element_ADS_AdiposeJitter$meanMethCoef_ADS_Adipose<(-0.25))
dim(incorrectPoints_ADS_Adipose_minus)
x_ADS_Adipose_minus<-subset(selected_reDinucleotide_Element_ADS_AdiposeJitter,selected_reDinucleotide_Element_ADS_AdiposeJitter$E_methCoef<=0.25)

correctPoints_ADS_Adipose<-subset(selected_reDinucleotide_Element_ADS_Adipose,abs(selected_reDinucleotide_Element_ADS_Adipose$E_methCoef-selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose)<=0.25)
# y_ADS_Adipose<-subset(selected_reDinucleotide_Element_ADS_AdiposeJitter,abs(selected_reDinucleotide_Element_ADS_AdiposeJitter$E_methCoef-selected_reDinucleotide_Element_ADS_AdiposeJitter$meanMethCoef_ADS_Adipose)<=0.25)
dim(correctPoints_ADS_Adipose)
y_ADS_Adipose<-subset(selected_reDinucleotide_Element_ADS_AdiposeJitter,selected_reDinucleotide_Element_ADS_AdiposeJitter$E_methCoef>0.25 & selected_reDinucleotide_Element_ADS_AdiposeJitter<0.75)

if (numbered=='TRUE'){
  png(paste(resultsDIR,"figure27_AllADS_AdiposeBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure27_AllADS_AdiposeBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_ADS_Adipose$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose,selected_reDinucleotide_Element_ADS_Adipose$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n",xlim=c(0,1),ylim=c(0,1))
points(x_ADS_Adipose_plus$meanMethCoef_ADS_Adipose,x_ADS_Adipose_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_ADS_Adipose_minus$meanMethCoef_ADS_Adipose,x_ADS_Adipose_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_ADS_Adipose$meanMethCoef_ADS_Adipose,y_ADS_Adipose$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

#points(incorrectPoints_ADS_Adipose_plus$meanMethCoef_ADS_Adipose,incorrectPoints_ADS_Adipose_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_ADS_Adipose_minus$meanMethCoef_ADS_Adipose,incorrectPoints_ADS_Adipose_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_ADS_Adipose$meanMethCoef_ADS_Adipose,correctPoints_ADS_Adipose$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if(numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_ADS_Adipose_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_ADS_Adipose_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_ADS_Adipose)[1],cex=1.2,col="red")
}else{}
abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

groupCPD_ADS_Adipose<-subset(correctPoints_ADS_Adipose,correctPoints_ADS_Adipose$E_methCoef<=0.25)
groupCPM_ADS_Adipose<-subset(correctPoints_ADS_Adipose,correctPoints_ADS_Adipose$E_methCoef>0.25 & correctPoints_ADS_Adipose$E_methCoef<0.75)
groupCPU_ADS_Adipose<-subset(correctPoints_ADS_Adipose,correctPoints_ADS_Adipose$E_methCoef>=0.75)

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
# lmADS_IPSCAll<-lm(selected_reDinucleotide_Element_ADS_IPSC$E_methCoef~selected_reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC)
# summary(lmADS_IPSCAll)

selected_reDinucleotide_Element_ADS_IPSCJitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_ADS_IPSC$E_methCoef,factor=25),meanMethCoef_ADS_IPSC=jitter(selected_reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC,factor=25)))

incorrectPoints_ADS_IPSC_plus<-subset(selected_reDinucleotide_Element_ADS_IPSC,selected_reDinucleotide_Element_ADS_IPSC$E_methCoef-selected_reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC>0.25)
#x_ADS_IPSC_plus<-subset(selected_reDinucleotide_Element_ADS_IPSCJitter,selected_reDinucleotide_Element_ADS_IPSCJitter$E_methCoef-selected_reDinucleotide_Element_ADS_IPSCJitter$meanMethCoef_ADS_IPSC>0.25)
dim(incorrectPoints_ADS_IPSC_plus)
x_ADS_IPSC_plus<-subset(selected_reDinucleotide_Element_ADS_IPSCJitter,selected_reDinucleotide_Element_ADS_IPSCJitter$E_methCoef>=0.75)

incorrectPoints_ADS_IPSC_minus<-subset(selected_reDinucleotide_Element_ADS_IPSC,selected_reDinucleotide_Element_ADS_IPSC$E_methCoef-selected_reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC<(-0.25))
#x_ADS_IPSC_minus<-subset(selected_reDinucleotide_Element_ADS_IPSCJitter,selected_reDinucleotide_Element_ADS_IPSCJitter$E_methCoef-selected_reDinucleotide_Element_ADS_IPSCJitter$meanMethCoef_ADS_IPSC<(-0.25))
dim(incorrectPoints_ADS_IPSC_minus)
x_ADS_IPSC_minus<-subset(selected_reDinucleotide_Element_ADS_IPSCJitter,selected_reDinucleotide_Element_ADS_IPSCJitter$E_methCoef<=0.25)

correctPoints_ADS_IPSC<-subset(selected_reDinucleotide_Element_ADS_IPSC,abs(selected_reDinucleotide_Element_ADS_IPSC$E_methCoef-selected_reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC)<=0.25)
# y_ADS_IPSC<-subset(selected_reDinucleotide_Element_ADS_IPSCJitter,abs(selected_reDinucleotide_Element_ADS_IPSCJitter$E_methCoef-selected_reDinucleotide_Element_ADS_IPSCJitter$meanMethCoef_ADS_IPSC)<=0.25)
dim(correctPoints_ADS_IPSC)
y_ADS_IPSC<-subset(selected_reDinucleotide_Element_ADS_IPSCJitter,selected_reDinucleotide_Element_ADS_IPSCJitter$E_methCoef>0.25 & selected_reDinucleotide_Element_ADS_IPSCJitter<0.75)

if (numbered=='TRUE'){
  png(paste(resultsDIR,"figure28_AllADS_IPSCBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
  png(paste(resultsDIR,"figure28_AllADS_IPSCBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_ADS_IPSC$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC,selected_reDinucleotide_Element_ADS_IPSC$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n",xlim=c(0,1),ylim=c(0,1))
points(x_ADS_IPSC_plus$meanMethCoef_ADS_IPSC,x_ADS_IPSC_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_ADS_IPSC_minus$meanMethCoef_ADS_IPSC,x_ADS_IPSC_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_ADS_IPSC$meanMethCoef_ADS_IPSC,y_ADS_IPSC$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

#points(incorrectPoints_ADS_IPSC_plus$meanMethCoef_ADS_IPSC,incorrectPoints_ADS_IPSC_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_ADS_IPSC_minus$meanMethCoef_ADS_IPSC,incorrectPoints_ADS_IPSC_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_ADS_IPSC$meanMethCoef_ADS_IPSC,correctPoints_ADS_IPSC$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if(numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_ADS_IPSC_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_ADS_IPSC_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_ADS_IPSC)[1],cex=1.2,col="red")
}else{}
abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

groupCPD_ADS_IPSC<-subset(correctPoints_ADS_IPSC,correctPoints_ADS_IPSC$E_methCoef<=0.25)
groupCPM_ADS_IPSC<-subset(correctPoints_ADS_IPSC,correctPoints_ADS_IPSC$E_methCoef>0.25 & correctPoints_ADS_IPSC$E_methCoef<0.75)
groupCPU_ADS_IPSC<-subset(correctPoints_ADS_IPSC,correctPoints_ADS_IPSC$E_methCoef>=0.75)

par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()


##################
### Figure 28 ###
##################


reDinucleotide_Element_FF_IPSC_19.11<-subset(reDinucleotide_Element,reDinucleotide_Element$CLine=="FF_IPSC_19.11")
reDinucleotide_Element_FF_IPSC_19.7<-subset(reDinucleotide_Element,reDinucleotide_Element$CLine=="FF_IPSC_19.7")
reDinucleotide_Element_FF_IPSC_6.9<-subset(reDinucleotide_Element,reDinucleotide_Element$CLine=="FF_IPSC_6.9")

meanMethCoef_FF_IPSC_19.11<-apply(reDinucleotide_Element_FF_IPSC_19.11,1,FUN=methCoefMean)
meanMethCoef_FF_IPSC_19.7<-apply(reDinucleotide_Element_FF_IPSC_19.7,1,FUN=methCoefMean)
meanMethCoef_FF_IPSC_6.9<-apply(reDinucleotide_Element_FF_IPSC_6.9,1,FUN=methCoefMean)

reDinucleotide_Element_FF_IPSC_19.11<-cbind(reDinucleotide_Element_FF_IPSC_19.11,meanMethCoef_FF_IPSC_19.11=as.numeric(as.vector(meanMethCoef_FF_IPSC_19.11)))
reDinucleotide_Element_FF_IPSC_19.7<-cbind(reDinucleotide_Element_FF_IPSC_19.7,meanMethCoef_FF_IPSC_19.7=as.numeric(as.vector(meanMethCoef_FF_IPSC_19.7)))
reDinucleotide_Element_FF_IPSC_6.9<-cbind(reDinucleotide_Element_FF_IPSC_6.9,meanMethCoef_FF_IPSC_6.9=as.numeric(as.vector(meanMethCoef_FF_IPSC_6.9)))

#################################
## With Informativeness filter ##
#################################

selected_reDinucleotide_Element_FF_IPSC_19.11<-subset(reDinucleotide_Element_FF_IPSC_19.11,
                                           reDinucleotide_Element_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11!=-1 & 
                                             reDinucleotide_Element_FF_IPSC_19.11$E_methCoef!=-1 & ((reDinucleotide_Element_FF_IPSC_19.11$E_posInf/2)/reDinucleotide_Element_FF_IPSC_19.11$E_nCG)>=0.25)
selected_reDinucleotide_Element_FF_IPSC_19.7<-subset(reDinucleotide_Element_FF_IPSC_19.7,
                                              reDinucleotide_Element_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7!=-1 & 
                                                reDinucleotide_Element_FF_IPSC_19.7$E_methCoef!=-1 & ((reDinucleotide_Element_FF_IPSC_19.7$E_posInf/2)/reDinucleotide_Element_FF_IPSC_19.7$E_nCG)>=0.25)
selected_reDinucleotide_Element_FF_IPSC_6.9<-subset(reDinucleotide_Element_FF_IPSC_6.9,
                                                    reDinucleotide_Element_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9!=-1 & 
                                                      reDinucleotide_Element_FF_IPSC_6.9$E_methCoef!=-1 & ((reDinucleotide_Element_FF_IPSC_6.9$E_posInf/2)/reDinucleotide_Element_FF_IPSC_6.9$E_nCG)>=0.25)


# write.table(cbind(selected_reDinucleotide_Element_FF_IPSC_19.11$id_In_RE,selected_reDinucleotide_Element_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11,selected_reDinucleotide_Element_FF_IPSC_19.11$id_In_Type,selected_reDinucleotide_Element_FF_IPSC_19.11$E_methCoef),file="reducedDataFF_IPSC_19.11_Filtered.txt",sep="\t",col.names=FALSE,row.names=F)
# write.table(cbind(selected_reDinucleotide_Element_FF_IPSC_19.7$id_In_RE,selected_reDinucleotide_Element_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7,selected_reDinucleotide_Element_FF_IPSC_19.7$id_In_Type,selected_reDinucleotide_Element_FF_IPSC_19.7$E_methCoef),file="reducedDataFF_IPSC_19.7_Filtered.txt",sep="\t",col.names=FALSE,row.names=F)
# write.table(cbind(selected_reDinucleotide_Element_FF_IPSC_6.9$id_In_RE,selected_reDinucleotide_Element_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9,selected_reDinucleotide_Element_FF_IPSC_6.9$id_In_Type,selected_reDinucleotide_Element_FF_IPSC_6.9$E_methCoef),file="reducedDataFF_IPSC_6.9_Filtered.txt",sep="\t",col.names=FALSE,row.names=F)

##########################
# Distance to extreme
t<-cbind(x=selected_reDinucleotide_Element_FF_IPSC_19.11$RE_position-selected_reDinucleotide_Element_FF_IPSC_19.11$E_chromStart,y=selected_reDinucleotide_Element_FF_IPSC_19.11$E_chromEnd-selected_reDinucleotide_Element_FF_IPSC_19.11$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_FF_IPSC_19.11<-cbind(selected_reDinucleotide_Element_FF_IPSC_19.11,dist2Ext=m)

t<-cbind(x=selected_reDinucleotide_Element_FF_IPSC_19.7$RE_position-selected_reDinucleotide_Element_FF_IPSC_19.7$E_chromStart,y=selected_reDinucleotide_Element_FF_IPSC_19.7$E_chromEnd-selected_reDinucleotide_Element_FF_IPSC_19.7$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_FF_IPSC_19.7<-cbind(selected_reDinucleotide_Element_FF_IPSC_19.7,dist2Ext=m)

t<-cbind(x=selected_reDinucleotide_Element_FF_IPSC_6.9$RE_position-selected_reDinucleotide_Element_FF_IPSC_6.9$E_chromStart,y=selected_reDinucleotide_Element_FF_IPSC_6.9$E_chromEnd-selected_reDinucleotide_Element_FF_IPSC_6.9$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_FF_IPSC_6.9<-cbind(selected_reDinucleotide_Element_FF_IPSC_6.9,dist2Ext=m)


#Obtain O/E ratio of the element

require(RMySQL)
con<-dbConnect(MySQL(),dbname="SolexaDB")
query<-"select id_cpgi_pk,obsExp FROM CpGisland"
rs1<-dbSendQuery(con,query)
OEdata<-fetch(rs1,n=-1) 
colnames(OEdata)=c("id_In_Type","E_OE")

selected_reDinucleotide_Element_FF_IPSC_19.11<-merge(selected_reDinucleotide_Element_FF_IPSC_19.11,OEdata,group="id_In_Type")
selected_reDinucleotide_Element_FF_IPSC_19.7<-merge(selected_reDinucleotide_Element_FF_IPSC_19.7,OEdata,group="id_In_Type")
selected_reDinucleotide_Element_FF_IPSC_6.9<-merge(selected_reDinucleotide_Element_FF_IPSC_6.9,OEdata,group="id_In_Type")

selected_reDinucleotide_Element_FF_IPSC_19.11<-selected_reDinucleotide_Element_FF_IPSC_19.11[,c(2:7,1,8:17)]
selected_reDinucleotide_Element_FF_IPSC_19.7<-selected_reDinucleotide_Element_FF_IPSC_19.7[,c(2:7,1,8:17)]
selected_reDinucleotide_Element_FF_IPSC_6.9<-selected_reDinucleotide_Element_FF_IPSC_6.9[,c(2:7,1,8:17)]

#Obtain length of the element

selected_reDinucleotide_Element_FF_IPSC_19.11<-cbind(selected_reDinucleotide_Element_FF_IPSC_19.11,
                                          cgiLength=selected_reDinucleotide_Element_FF_IPSC_19.11$E_chromEnd-
                                            selected_reDinucleotide_Element_FF_IPSC_19.11$E_chromStart)

selected_reDinucleotide_Element_FF_IPSC_19.7<-cbind(selected_reDinucleotide_Element_FF_IPSC_19.7,
                                             cgiLength=selected_reDinucleotide_Element_FF_IPSC_19.7$E_chromEnd-
                                               selected_reDinucleotide_Element_FF_IPSC_19.7$E_chromStart)
selected_reDinucleotide_Element_FF_IPSC_6.9<-cbind(selected_reDinucleotide_Element_FF_IPSC_6.9,
                                                   cgiLength=selected_reDinucleotide_Element_FF_IPSC_6.9$E_chromEnd-
                                                     selected_reDinucleotide_Element_FF_IPSC_6.9$E_chromStart)

############
#FF_IPSC_19.11

# lmFF_IPSC_19.11Filtered<-lm(selected_reDinucleotide_Element_FF_IPSC_19.11$E_methCoef~selected_reDinucleotide_Element_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11)
# summary(lmFF_IPSC_19.11Filtered)
selected_reDinucleotide_Element_FF_IPSC_19.11Jitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_FF_IPSC_19.11$E_methCoef,factor=25),meanMethCoef_FF_IPSC_19.11=jitter(selected_reDinucleotide_Element_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11,factor=25)))

  
incorrectPoints_FF_IPSC_19.11_plus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.11,selected_reDinucleotide_Element_FF_IPSC_19.11$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11>0.25)
#x_FF_IPSC_19.11_plus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.11Jitter,selected_reDinucleotide_Element_FF_IPSC_19.11Jitter$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.11Jitter$meanMethCoef_FF_IPSC_19.11>0.25)
dim(incorrectPoints_FF_IPSC_19.11_plus)
x_FF_IPSC_19.11_plus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.11Jitter,selected_reDinucleotide_Element_FF_IPSC_19.11Jitter$E_methCoef>=0.75)

incorrectPoints_FF_IPSC_19.11_minus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.11,selected_reDinucleotide_Element_FF_IPSC_19.11$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11<(-0.25))
#x_FF_IPSC_19.11_minus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.11Jitter,selected_reDinucleotide_Element_FF_IPSC_19.11Jitter$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.11Jitter$meanMethCoef_FF_IPSC_19.11<(-0.25))
dim(incorrectPoints_FF_IPSC_19.11_minus)
x_FF_IPSC_19.11_minus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.11Jitter,selected_reDinucleotide_Element_FF_IPSC_19.11Jitter$E_methCoef<=0.25)

correctPoints_FF_IPSC_19.11<-subset(selected_reDinucleotide_Element_FF_IPSC_19.11,abs(selected_reDinucleotide_Element_FF_IPSC_19.11$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11)<=0.25)
# y_FF_IPSC_19.11<-subset(selected_reDinucleotide_Element_FF_IPSC_19.11Jitter,abs(selected_reDinucleotide_Element_FF_IPSC_19.11Jitter$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.11Jitter$meanMethCoef_FF_IPSC_19.11)<=0.25)
dim(correctPoints_FF_IPSC_19.11)
y_FF_IPSC_19.11<-subset(selected_reDinucleotide_Element_FF_IPSC_19.11Jitter,selected_reDinucleotide_Element_FF_IPSC_19.11Jitter$E_methCoef>0.25 & selected_reDinucleotide_Element_FF_IPSC_19.11Jitter<0.75)

if (numbered=='TRUE'){
  png(paste(resultsDIR,"figure28_FF_IPSC_19.11FilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
} else {
  png(paste(resultsDIR,"figure28_FF_IPSC_19.11FilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_FF_IPSC_19.11$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11,selected_reDinucleotide_Element_FF_IPSC_19.11$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n",xlim=c(0,1),ylim=c(0,1))

## Plot the points with some noise (jitter)
points(x_FF_IPSC_19.11_plus$meanMethCoef_FF_IPSC_19.11,x_FF_IPSC_19.11_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_FF_IPSC_19.11_minus$meanMethCoef_FF_IPSC_19.11,x_FF_IPSC_19.11_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11,y_FF_IPSC_19.11$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# Plot the points without noise
#points(incorrectPoints_FF_IPSC_19.11_plus$meanMethCoef_FF_IPSC_19.11,incorrectPoints_FF_IPSC_19.11_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_FF_IPSC_19.11_minus$meanMethCoef_FF_IPSC_19.11,incorrectPoints_FF_IPSC_19.11_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11,correctPoints_FF_IPSC_19.11$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if (numbered=='TRUE'){
# Adding the number of points to the plot
text(0.2,0.8,dim(incorrectPoints_FF_IPSC_19.11_plus)[1],cex=1.2,col="red")
text(0.8,0.2,dim(incorrectPoints_FF_IPSC_19.11_minus)[1],cex=1.2,col="red")
text(0.5,0.5,dim(correctPoints_FF_IPSC_19.11)[1],cex=1.2,col="red")
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


# FF_IPSC_19.7

# lmFF_IPSC_19.7Filtered<-lm(selected_reDinucleotide_Element_FF_IPSC_19.7$E_methCoef~selected_reDinucleotide_Element_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7)
# summary(lmFF_IPSC_19.7Filtered)

selected_reDinucleotide_Element_FF_IPSC_19.7Jitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_FF_IPSC_19.7$E_methCoef,factor=25),meanMethCoef_FF_IPSC_19.7=jitter(selected_reDinucleotide_Element_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7,factor=25)))

incorrectPoints_FF_IPSC_19.7_plus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.7,selected_reDinucleotide_Element_FF_IPSC_19.7$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7>0.25)
#x_FF_IPSC_19.7_plus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.7Jitter,selected_reDinucleotide_Element_FF_IPSC_19.7Jitter$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.7Jitter$meanMethCoef_FF_IPSC_19.7>0.25)
dim(incorrectPoints_FF_IPSC_19.7_plus)
x_FF_IPSC_19.7_plus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.7Jitter,selected_reDinucleotide_Element_FF_IPSC_19.7Jitter$E_methCoef>=0.75)

incorrectPoints_FF_IPSC_19.7_minus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.7,selected_reDinucleotide_Element_FF_IPSC_19.7$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7<(-0.25))
#x_FF_IPSC_19.7_minus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.7Jitter,selected_reDinucleotide_Element_FF_IPSC_19.7Jitter$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.7Jitter$meanMethCoef_FF_IPSC_19.7<(-0.25))
dim(incorrectPoints_FF_IPSC_19.7_minus)
x_FF_IPSC_19.7_minus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.7Jitter,selected_reDinucleotide_Element_FF_IPSC_19.7Jitter$E_methCoef<=0.25)

correctPoints_FF_IPSC_19.7<-subset(selected_reDinucleotide_Element_FF_IPSC_19.7,abs(selected_reDinucleotide_Element_FF_IPSC_19.7$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7)<=0.25)
#y_FF_IPSC_19.7<-subset(selected_reDinucleotide_Element_FF_IPSC_19.7Jitter,abs(selected_reDinucleotide_Element_FF_IPSC_19.7Jitter$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.7Jitter$meanMethCoef_FF_IPSC_19.7)<=0.25)
dim(correctPoints_FF_IPSC_19.7)
y_FF_IPSC_19.7<-subset(selected_reDinucleotide_Element_FF_IPSC_19.7Jitter,selected_reDinucleotide_Element_FF_IPSC_19.7Jitter$E_methCoef>0.25 & selected_reDinucleotide_Element_FF_IPSC_19.7Jitter<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure28_FF_IPSC_19.7FilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
  png(paste(resultsDIR,"figure28_FF_IPSC_19.7FilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_FF_IPSC_19.7$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7,selected_reDinucleotide_Element_FF_IPSC_19.7$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n",xlim=c(0,1),ylim=c(0,1))

## Plot the points with some noise (jitter)
points(x_FF_IPSC_19.7_plus$meanMethCoef_FF_IPSC_19.7,x_FF_IPSC_19.7_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_FF_IPSC_19.7_minus$meanMethCoef_FF_IPSC_19.7,x_FF_IPSC_19.7_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7,y_FF_IPSC_19.7$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# Plot the point without noise
#points(incorrectPoints_FF_IPSC_19.7_plus$meanMethCoef_FF_IPSC_19.7,incorrectPoints_FF_IPSC_19.7_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_FF_IPSC_19.7_minus$meanMethCoef_FF_IPSC_19.7,incorrectPoints_FF_IPSC_19.7_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7,correctPoints_FF_IPSC_19.7$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if(numbered=='TRUE'){
# Adding the number of points to the plot
text(0.2,0.8,dim(incorrectPoints_FF_IPSC_19.7_plus)[1],cex=1.2,col="red")
text(0.8,0.2,dim(incorrectPoints_FF_IPSC_19.7_minus)[1],cex=1.2,col="red")
text(0.5,0.5,dim(correctPoints_FF_IPSC_19.7)[1],cex=1.2,col="red")
}else{}
abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

groupCPD_FF_IPSC_19.7<-subset(correctPoints_FF_IPSC_19.7,correctPoints_FF_IPSC_19.7$E_methCoef<=0.25)
groupCPM_FF_IPSC_19.7<-subset(correctPoints_FF_IPSC_19.7,correctPoints_FF_IPSC_19.7$E_methCoef>0.25 & correctPoints_FF_IPSC_19.7$E_methCoef<0.75)
groupCPU_FF_IPSC_19.7<-subset(correctPoints_FF_IPSC_19.7,correctPoints_FF_IPSC_19.7$E_methCoef>=0.75)

par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()

# FF_IPSC_6.9
# lmFF_IPSC_6.9Filtered<-lm(selected_reDinucleotide_Element_FF_IPSC_6.9$E_methCoef~selected_reDinucleotide_Element_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9)
# summary(lmFF_IPSC_6.9Filtered)

selected_reDinucleotide_Element_FF_IPSC_6.9Jitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_FF_IPSC_6.9$E_methCoef,factor=25),meanMethCoef_FF_IPSC_6.9=jitter(selected_reDinucleotide_Element_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9,factor=25)))

incorrectPoints_FF_IPSC_6.9_plus<-subset(selected_reDinucleotide_Element_FF_IPSC_6.9,selected_reDinucleotide_Element_FF_IPSC_6.9$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9>0.25)
#x_FF_IPSC_6.9_plus<-subset(selected_reDinucleotide_Element_FF_IPSC_6.9Jitter,selected_reDinucleotide_Element_FF_IPSC_6.9Jitter$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_6.9Jitter$meanMethCoef_FF_IPSC_6.9>0.25)
dim(incorrectPoints_FF_IPSC_6.9_plus)
x_FF_IPSC_6.9_plus<-subset(selected_reDinucleotide_Element_FF_IPSC_6.9Jitter,selected_reDinucleotide_Element_FF_IPSC_6.9Jitter$E_methCoef>=0.75)

incorrectPoints_FF_IPSC_6.9_minus<-subset(selected_reDinucleotide_Element_FF_IPSC_6.9,selected_reDinucleotide_Element_FF_IPSC_6.9$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9<(-0.25))
#x_FF_IPSC_6.9_minus<-subset(selected_reDinucleotide_Element_FF_IPSC_6.9Jitter,selected_reDinucleotide_Element_FF_IPSC_6.9Jitter$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_6.9Jitter$meanMethCoef_FF_IPSC_6.9<(-0.25))
dim(incorrectPoints_FF_IPSC_6.9_minus)
x_FF_IPSC_6.9_minus<-subset(selected_reDinucleotide_Element_FF_IPSC_6.9Jitter,selected_reDinucleotide_Element_FF_IPSC_6.9Jitter$E_methCoef<=0.25)

correctPoints_FF_IPSC_6.9<-subset(selected_reDinucleotide_Element_FF_IPSC_6.9,abs(selected_reDinucleotide_Element_FF_IPSC_6.9$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9)<=0.25)
#y_FF_IPSC_6.9<-subset(selected_reDinucleotide_Element_FF_IPSC_6.9Jitter,abs(selected_reDinucleotide_Element_FF_IPSC_6.9Jitter$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_6.9Jitter$meanMethCoef_FF_IPSC_6.9)<=0.25)
dim(correctPoints_FF_IPSC_6.9)
y_FF_IPSC_6.9<-subset(selected_reDinucleotide_Element_FF_IPSC_6.9Jitter,selected_reDinucleotide_Element_FF_IPSC_6.9Jitter$E_methCoef>0.25 & selected_reDinucleotide_Element_FF_IPSC_6.9Jitter<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure28_FF_IPSC_6.9FilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
  png(paste(resultsDIR,"figure28_FF_IPSC_6.9FilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_FF_IPSC_6.9$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9,selected_reDinucleotide_Element_FF_IPSC_6.9$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n",xlim=c(0,1),ylim=c(0,1))

## Plot the points with some noise (jitter)
points(x_FF_IPSC_6.9_plus$meanMethCoef_FF_IPSC_6.9,x_FF_IPSC_6.9_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_FF_IPSC_6.9_minus$meanMethCoef_FF_IPSC_6.9,x_FF_IPSC_6.9_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9,y_FF_IPSC_6.9$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# Plot the point without noise
#points(incorrectPoints_FF_IPSC_6.9_plus$meanMethCoef_FF_IPSC_6.9,incorrectPoints_FF_IPSC_6.9_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_FF_IPSC_6.9_minus$meanMethCoef_FF_IPSC_6.9,incorrectPoints_FF_IPSC_6.9_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9,correctPoints_FF_IPSC_6.9$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if(numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_FF_IPSC_6.9_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_FF_IPSC_6.9_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_FF_IPSC_6.9)[1],cex=1.2,col="red")
}else{}
abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

groupCPD_FF_IPSC_6.9<-subset(correctPoints_FF_IPSC_6.9,correctPoints_FF_IPSC_6.9$E_methCoef<=0.25)
groupCPM_FF_IPSC_6.9<-subset(correctPoints_FF_IPSC_6.9,correctPoints_FF_IPSC_6.9$E_methCoef>0.25 & correctPoints_FF_IPSC_6.9$E_methCoef<0.75)
groupCPU_FF_IPSC_6.9<-subset(correctPoints_FF_IPSC_6.9,correctPoints_FF_IPSC_6.9$E_methCoef>=0.75)

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

selected_reDinucleotide_Element_FF_IPSC_19.11<-subset(reDinucleotide_Element_FF_IPSC_19.11,
                                           reDinucleotide_Element_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11!=-1 & 
                                             reDinucleotide_Element_FF_IPSC_19.11$E_methCoef!=-1)
selected_reDinucleotide_Element_FF_IPSC_19.7<-subset(reDinucleotide_Element_FF_IPSC_19.7,
                                              reDinucleotide_Element_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7!=-1 & 
                                                reDinucleotide_Element_FF_IPSC_19.7$E_methCoef!=-1)
selected_reDinucleotide_Element_FF_IPSC_6.9<-subset(reDinucleotide_Element_FF_IPSC_6.9,
                                                    reDinucleotide_Element_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9!=-1 & 
                                                      reDinucleotide_Element_FF_IPSC_6.9$E_methCoef!=-1)

# write.table(cbind(selected_reDinucleotide_Element_FF_IPSC_19.11$id_In_RE,selected_reDinucleotide_Element_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11,selected_reDinucleotide_Element_FF_IPSC_19.11$id_In_Type,selected_reDinucleotide_Element_FF_IPSC_19.11$E_methCoef),file="reducedDataFF_IPSC_19.11_All.txt",sep="\t",col.names=FALSE,row.names=F)
# write.table(cbind(selected_reDinucleotide_Element_FF_IPSC_19.7$id_In_RE,selected_reDinucleotide_Element_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7,selected_reDinucleotide_Element_FF_IPSC_19.7$id_In_Type,selected_reDinucleotide_Element_FF_IPSC_19.7$E_methCoef),file="reducedDataFF_IPSC_19.7_All.txt",sep="\t",col.names=FALSE,row.names=F)
# write.table(cbind(selected_reDinucleotide_Element_FF_IPSC_6.9$id_In_RE,selected_reDinucleotide_Element_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9,selected_reDinucleotide_Element_FF_IPSC_6.9$id_In_Type,selected_reDinucleotide_Element_FF_IPSC_6.9$E_methCoef),file="reducedDataFF_IPSC_6.9_All.txt",sep="\t",col.names=FALSE,row.names=F)



#######################
# Distance to extreme
t<-cbind(x=selected_reDinucleotide_Element_FF_IPSC_19.11$RE_position-selected_reDinucleotide_Element_FF_IPSC_19.11$E_chromStart,y=selected_reDinucleotide_Element_FF_IPSC_19.11$E_chromEnd-selected_reDinucleotide_Element_FF_IPSC_19.11$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_FF_IPSC_19.11<-cbind(selected_reDinucleotide_Element_FF_IPSC_19.11,dist2Ext=m)

t<-cbind(x=selected_reDinucleotide_Element_FF_IPSC_19.7$RE_position-selected_reDinucleotide_Element_FF_IPSC_19.7$E_chromStart,y=selected_reDinucleotide_Element_FF_IPSC_19.7$E_chromEnd-selected_reDinucleotide_Element_FF_IPSC_19.7$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_FF_IPSC_19.7<-cbind(selected_reDinucleotide_Element_FF_IPSC_19.7,dist2Ext=m)


t<-cbind(x=selected_reDinucleotide_Element_FF_IPSC_6.9$RE_position-selected_reDinucleotide_Element_FF_IPSC_6.9$E_chromStart,y=selected_reDinucleotide_Element_FF_IPSC_6.9$E_chromEnd-selected_reDinucleotide_Element_FF_IPSC_6.9$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_FF_IPSC_6.9<-cbind(selected_reDinucleotide_Element_FF_IPSC_6.9,dist2Ext=m)

#Obtain O/E ratio of the element

require(RMySQL)
con<-dbConnect(MySQL(),dbname="SolexaDB")
query<-"select id_cpgi_pk,obsExp FROM CpGisland"
rs1<-dbSendQuery(con,query)
OEdata<-fetch(rs1,n=-1) 
colnames(OEdata)=c("id_In_Type","E_OE")

selected_reDinucleotide_Element_FF_IPSC_19.11<-merge(selected_reDinucleotide_Element_FF_IPSC_19.11,OEdata,group="id_In_Type")
selected_reDinucleotide_Element_FF_IPSC_19.7<-merge(selected_reDinucleotide_Element_FF_IPSC_19.7,OEdata,group="id_In_Type")
selected_reDinucleotide_Element_FF_IPSC_6.9<-merge(selected_reDinucleotide_Element_FF_IPSC_6.9,OEdata,group="id_In_Type")


selected_reDinucleotide_Element_FF_IPSC_19.11<-selected_reDinucleotide_Element_FF_IPSC_19.11[,c(2:7,1,8:17)]
selected_reDinucleotide_Element_FF_IPSC_19.7<-selected_reDinucleotide_Element_FF_IPSC_19.7[,c(2:7,1,8:17)]
selected_reDinucleotide_Element_FF_IPSC_6.9<-selected_reDinucleotide_Element_FF_IPSC_6.9[,c(2:7,1,8:17)]

#Obtain length of the element

selected_reDinucleotide_Element_FF_IPSC_19.11<-cbind(selected_reDinucleotide_Element_FF_IPSC_19.11,
                                          cgiLength=selected_reDinucleotide_Element_FF_IPSC_19.11$E_chromEnd-
                                            selected_reDinucleotide_Element_FF_IPSC_19.11$E_chromStart)

selected_reDinucleotide_Element_FF_IPSC_19.7<-cbind(selected_reDinucleotide_Element_FF_IPSC_19.7,
                                             cgiLength=selected_reDinucleotide_Element_FF_IPSC_19.7$E_chromEnd-
                                               selected_reDinucleotide_Element_FF_IPSC_19.7$E_chromStart)
selected_reDinucleotide_Element_FF_IPSC_6.9<-cbind(selected_reDinucleotide_Element_FF_IPSC_6.9,
                                                   cgiLength=selected_reDinucleotide_Element_FF_IPSC_6.9$E_chromEnd-
                                                     selected_reDinucleotide_Element_FF_IPSC_6.9$E_chromStart)

############
#FF_IPSC_19.11
# lmFF_IPSC_19.11All<-lm(selected_reDinucleotide_Element_FF_IPSC_19.11$E_methCoef~selected_reDinucleotide_Element_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11)
# summary(lmFF_IPSC_19.11All)
selected_reDinucleotide_Element_FF_IPSC_19.11Jitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_FF_IPSC_19.11$E_methCoef,factor=25),meanMethCoef_FF_IPSC_19.11=jitter(selected_reDinucleotide_Element_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11,factor=25)))

incorrectPoints_FF_IPSC_19.11_plus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.11,selected_reDinucleotide_Element_FF_IPSC_19.11$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11>0.25)
#x_FF_IPSC_19.11_plus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.11Jitter,selected_reDinucleotide_Element_FF_IPSC_19.11Jitter$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.11Jitter$meanMethCoef_FF_IPSC_19.11>0.25)
dim(incorrectPoints_FF_IPSC_19.11_plus)
x_FF_IPSC_19.11_plus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.11Jitter,selected_reDinucleotide_Element_FF_IPSC_19.11Jitter$E_methCoef>=0.75)

incorrectPoints_FF_IPSC_19.11_minus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.11,selected_reDinucleotide_Element_FF_IPSC_19.11$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11<(-0.25))
#x_FF_IPSC_19.11_minus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.11Jitter,selected_reDinucleotide_Element_FF_IPSC_19.11Jitter$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.11Jitter$meanMethCoef_FF_IPSC_19.11<(-0.25))
dim(incorrectPoints_FF_IPSC_19.11_minus)
x_FF_IPSC_19.11_minus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.11Jitter,selected_reDinucleotide_Element_FF_IPSC_19.11Jitter$E_methCoef<=0.25)

correctPoints_FF_IPSC_19.11<-subset(selected_reDinucleotide_Element_FF_IPSC_19.11,abs(selected_reDinucleotide_Element_FF_IPSC_19.11$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11)<=0.25)
# y_FF_IPSC_19.11<-subset(selected_reDinucleotide_Element_FF_IPSC_19.11Jitter,abs(selected_reDinucleotide_Element_FF_IPSC_19.11Jitter$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.11Jitter$meanMethCoef_FF_IPSC_19.11)<=0.25)
dim(correctPoints_FF_IPSC_19.11)
y_FF_IPSC_19.11<-subset(selected_reDinucleotide_Element_FF_IPSC_19.11Jitter,selected_reDinucleotide_Element_FF_IPSC_19.11Jitter$E_methCoef>0.25 & selected_reDinucleotide_Element_FF_IPSC_19.11Jitter<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure28_FF_IPSC_19.11AllBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure28_FF_IPSC_19.11AllBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_FF_IPSC_19.11$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11,selected_reDinucleotide_Element_FF_IPSC_19.11$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n",xlim=c(0,1),ylim=c(0,1))
points(x_FF_IPSC_19.11_plus$meanMethCoef_FF_IPSC_19.11,x_FF_IPSC_19.11_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_FF_IPSC_19.11_minus$meanMethCoef_FF_IPSC_19.11,x_FF_IPSC_19.11_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11,y_FF_IPSC_19.11$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

#points(incorrectPoints_FF_IPSC_19.11_plus$meanMethCoef_FF_IPSC_19.11,incorrectPoints_FF_IPSC_19.11_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_FF_IPSC_19.11_minus$meanMethCoef_FF_IPSC_19.11,incorrectPoints_FF_IPSC_19.11_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_FF_IPSC_19.11$meanMethCoef_FF_IPSC_19.11,correctPoints_FF_IPSC_19.11$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if (numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_FF_IPSC_19.11_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_FF_IPSC_19.11_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_FF_IPSC_19.11)[1],cex=1.2,col="red")
}else{}

abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

groupCPD_FF_IPSC_19.11<-subset(correctPoints_FF_IPSC_19.11,correctPoints_FF_IPSC_19.11$E_methCoef<=0.25)
groupCPM_FF_IPSC_19.11<-subset(correctPoints_FF_IPSC_19.11,correctPoints_FF_IPSC_19.11$E_methCoef>0.25 & correctPoints_FF_IPSC_19.11$E_methCoef<0.75)
groupCPU_FF_IPSC_19.11<-subset(correctPoints_FF_IPSC_19.11,correctPoints_FF_IPSC_19.11$E_methCoef>=0.75)



par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()

# FF_IPSC_19.7

# lmFF_IPSC_19.7All<-lm(selected_reDinucleotide_Element_FF_IPSC_19.7$E_methCoef~selected_reDinucleotide_Element_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7)
# summary(lmFF_IPSC_19.7All)

selected_reDinucleotide_Element_FF_IPSC_19.7Jitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_FF_IPSC_19.7$E_methCoef,factor=25),meanMethCoef_FF_IPSC_19.7=jitter(selected_reDinucleotide_Element_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7,factor=25)))

incorrectPoints_FF_IPSC_19.7_plus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.7,selected_reDinucleotide_Element_FF_IPSC_19.7$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7>0.25)
#x_FF_IPSC_19.7_plus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.7Jitter,selected_reDinucleotide_Element_FF_IPSC_19.7Jitter$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.7Jitter$meanMethCoef_FF_IPSC_19.7>0.25)
dim(incorrectPoints_FF_IPSC_19.7_plus)
x_FF_IPSC_19.7_plus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.7Jitter,selected_reDinucleotide_Element_FF_IPSC_19.7Jitter$E_methCoef>=0.75)

incorrectPoints_FF_IPSC_19.7_minus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.7,selected_reDinucleotide_Element_FF_IPSC_19.7$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7<(-0.25))
#x_FF_IPSC_19.7_minus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.7Jitter,selected_reDinucleotide_Element_FF_IPSC_19.7Jitter$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.7Jitter$meanMethCoef_FF_IPSC_19.7<(-0.25))
dim(incorrectPoints_FF_IPSC_19.7_minus)
x_FF_IPSC_19.7_minus<-subset(selected_reDinucleotide_Element_FF_IPSC_19.7Jitter,selected_reDinucleotide_Element_FF_IPSC_19.7Jitter$E_methCoef<=0.25)

correctPoints_FF_IPSC_19.7<-subset(selected_reDinucleotide_Element_FF_IPSC_19.7,abs(selected_reDinucleotide_Element_FF_IPSC_19.7$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7)<=0.25)
# y_FF_IPSC_19.7<-subset(selected_reDinucleotide_Element_FF_IPSC_19.7Jitter,abs(selected_reDinucleotide_Element_FF_IPSC_19.7Jitter$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_19.7Jitter$meanMethCoef_FF_IPSC_19.7)<=0.25)
dim(correctPoints_FF_IPSC_19.7)
y_FF_IPSC_19.7<-subset(selected_reDinucleotide_Element_FF_IPSC_19.7Jitter,selected_reDinucleotide_Element_FF_IPSC_19.7Jitter$E_methCoef>0.25 & selected_reDinucleotide_Element_FF_IPSC_19.7Jitter<0.75)

if (numbered=='TRUE'){
  png(paste(resultsDIR,"figure28_AllFF_IPSC_19.7BNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure28_AllFF_IPSC_19.7BN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_FF_IPSC_19.7$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7,selected_reDinucleotide_Element_FF_IPSC_19.7$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n",xlim=c(0,1),ylim=c(0,1))
points(x_FF_IPSC_19.7_plus$meanMethCoef_FF_IPSC_19.7,x_FF_IPSC_19.7_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_FF_IPSC_19.7_minus$meanMethCoef_FF_IPSC_19.7,x_FF_IPSC_19.7_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7,y_FF_IPSC_19.7$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

#points(incorrectPoints_FF_IPSC_19.7_plus$meanMethCoef_FF_IPSC_19.7,incorrectPoints_FF_IPSC_19.7_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_FF_IPSC_19.7_minus$meanMethCoef_FF_IPSC_19.7,incorrectPoints_FF_IPSC_19.7_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_FF_IPSC_19.7$meanMethCoef_FF_IPSC_19.7,correctPoints_FF_IPSC_19.7$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if(numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_FF_IPSC_19.7_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_FF_IPSC_19.7_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_FF_IPSC_19.7)[1],cex=1.2,col="red")
}else{}
abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

groupCPD_FF_IPSC_19.7<-subset(correctPoints_FF_IPSC_19.7,correctPoints_FF_IPSC_19.7$E_methCoef<=0.25)
groupCPM_FF_IPSC_19.7<-subset(correctPoints_FF_IPSC_19.7,correctPoints_FF_IPSC_19.7$E_methCoef>0.25 & correctPoints_FF_IPSC_19.7$E_methCoef<0.75)
groupCPU_FF_IPSC_19.7<-subset(correctPoints_FF_IPSC_19.7,correctPoints_FF_IPSC_19.7$E_methCoef>=0.75)

par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()


# FF_IPSC_6.9
# lmFF_IPSC_6.9All<-lm(selected_reDinucleotide_Element_FF_IPSC_6.9$E_methCoef~selected_reDinucleotide_Element_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9)
# summary(lmFF_IPSC_6.9All)

selected_reDinucleotide_Element_FF_IPSC_6.9Jitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_FF_IPSC_6.9$E_methCoef,factor=25),meanMethCoef_FF_IPSC_6.9=jitter(selected_reDinucleotide_Element_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9,factor=25)))

incorrectPoints_FF_IPSC_6.9_plus<-subset(selected_reDinucleotide_Element_FF_IPSC_6.9,selected_reDinucleotide_Element_FF_IPSC_6.9$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9>0.25)
#x_FF_IPSC_6.9_plus<-subset(selected_reDinucleotide_Element_FF_IPSC_6.9Jitter,selected_reDinucleotide_Element_FF_IPSC_6.9Jitter$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_6.9Jitter$meanMethCoef_FF_IPSC_6.9>0.25)
dim(incorrectPoints_FF_IPSC_6.9_plus)
x_FF_IPSC_6.9_plus<-subset(selected_reDinucleotide_Element_FF_IPSC_6.9Jitter,selected_reDinucleotide_Element_FF_IPSC_6.9Jitter$E_methCoef>=0.75)

incorrectPoints_FF_IPSC_6.9_minus<-subset(selected_reDinucleotide_Element_FF_IPSC_6.9,selected_reDinucleotide_Element_FF_IPSC_6.9$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9<(-0.25))
#x_FF_IPSC_6.9_minus<-subset(selected_reDinucleotide_Element_FF_IPSC_6.9Jitter,selected_reDinucleotide_Element_FF_IPSC_6.9Jitter$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_6.9Jitter$meanMethCoef_FF_IPSC_6.9<(-0.25))
dim(incorrectPoints_FF_IPSC_6.9_minus)
x_FF_IPSC_6.9_minus<-subset(selected_reDinucleotide_Element_FF_IPSC_6.9Jitter,selected_reDinucleotide_Element_FF_IPSC_6.9Jitter$E_methCoef<=0.25)

correctPoints_FF_IPSC_6.9<-subset(selected_reDinucleotide_Element_FF_IPSC_6.9,abs(selected_reDinucleotide_Element_FF_IPSC_6.9$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9)<=0.25)
# y_FF_IPSC_6.9<-subset(selected_reDinucleotide_Element_FF_IPSC_6.9Jitter,abs(selected_reDinucleotide_Element_FF_IPSC_6.9Jitter$E_methCoef-selected_reDinucleotide_Element_FF_IPSC_6.9Jitter$meanMethCoef_FF_IPSC_6.9)<=0.25)
dim(correctPoints_FF_IPSC_6.9)
y_FF_IPSC_6.9<-subset(selected_reDinucleotide_Element_FF_IPSC_6.9Jitter,selected_reDinucleotide_Element_FF_IPSC_6.9Jitter$E_methCoef>0.25 & selected_reDinucleotide_Element_FF_IPSC_6.9Jitter<0.75)

if (numbered=='TRUE'){
  png(paste(resultsDIR,"figure28_AllFF_IPSC_6.9BNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
  png(paste(resultsDIR,"figure28_AllFF_IPSC_6.9BN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_FF_IPSC_6.9$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9,selected_reDinucleotide_Element_FF_IPSC_6.9$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n",xlim=c(0,1),ylim=c(0,1))
points(x_FF_IPSC_6.9_plus$meanMethCoef_FF_IPSC_6.9,x_FF_IPSC_6.9_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_FF_IPSC_6.9_minus$meanMethCoef_FF_IPSC_6.9,x_FF_IPSC_6.9_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9,y_FF_IPSC_6.9$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

#points(incorrectPoints_FF_IPSC_6.9_plus$meanMethCoef_FF_IPSC_6.9,incorrectPoints_FF_IPSC_6.9_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_FF_IPSC_6.9_minus$meanMethCoef_FF_IPSC_6.9,incorrectPoints_FF_IPSC_6.9_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_FF_IPSC_6.9$meanMethCoef_FF_IPSC_6.9,correctPoints_FF_IPSC_6.9$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if(numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_FF_IPSC_6.9_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_FF_IPSC_6.9_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_FF_IPSC_6.9)[1],cex=1.2,col="red")
}else{}
abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

groupCPD_FF_IPSC_6.9<-subset(correctPoints_FF_IPSC_6.9,correctPoints_FF_IPSC_6.9$E_methCoef<=0.25)
groupCPM_FF_IPSC_6.9<-subset(correctPoints_FF_IPSC_6.9,correctPoints_FF_IPSC_6.9$E_methCoef>0.25 & correctPoints_FF_IPSC_6.9$E_methCoef<0.75)
groupCPU_FF_IPSC_6.9<-subset(correctPoints_FF_IPSC_6.9,correctPoints_FF_IPSC_6.9$E_methCoef>=0.75)

par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()


##################
### Figure 28 ###
##################

reDinucleotide_Element_IMR90_IPSC<-subset(reDinucleotide_Element,reDinucleotide_Element$CLine=="IMR90_IPSC")
reDinucleotide_Element_H9<-subset(reDinucleotide_Element,reDinucleotide_Element$CLine=="H9")
reDinucleotide_Element_H9<-subset(reDinucleotide_Element,reDinucleotide_Element$CLine=="H9")

meanMethCoef_IMR90_IPSC<-apply(reDinucleotide_Element_IMR90_IPSC,1,FUN=methCoefMean)
meanMethCoef_H9<-apply(reDinucleotide_Element_H9,1,FUN=methCoefMean)
meanMethCoef_H9<-apply(reDinucleotide_Element_H9,1,FUN=methCoefMean)

reDinucleotide_Element_IMR90_IPSC<-cbind(reDinucleotide_Element_IMR90_IPSC,meanMethCoef_IMR90_IPSC=as.numeric(as.vector(meanMethCoef_IMR90_IPSC)))
reDinucleotide_Element_H9<-cbind(reDinucleotide_Element_H9,meanMethCoef_H9=as.numeric(as.vector(meanMethCoef_H9)))
reDinucleotide_Element_H9<-cbind(reDinucleotide_Element_H9,meanMethCoef_H9=as.numeric(as.vector(meanMethCoef_H9)))

#################################
## With Informativeness filter ##
#################################

selected_reDinucleotide_Element_IMR90_IPSC<-subset(reDinucleotide_Element_IMR90_IPSC,
                                           reDinucleotide_Element_IMR90_IPSC$meanMethCoef_IMR90_IPSC!=-1 & 
                                             reDinucleotide_Element_IMR90_IPSC$E_methCoef!=-1 & ((reDinucleotide_Element_IMR90_IPSC$E_posInf/2)/reDinucleotide_Element_IMR90_IPSC$E_nCG)>=0.25)
selected_reDinucleotide_Element_H9<-subset(reDinucleotide_Element_H9,
                                              reDinucleotide_Element_H9$meanMethCoef_H9!=-1 & 
                                                reDinucleotide_Element_H9$E_methCoef!=-1 & ((reDinucleotide_Element_H9$E_posInf/2)/reDinucleotide_Element_H9$E_nCG)>=0.25)
selected_reDinucleotide_Element_H9<-subset(reDinucleotide_Element_H9,
                                                    reDinucleotide_Element_H9$meanMethCoef_H9!=-1 & 
                                                      reDinucleotide_Element_H9$E_methCoef!=-1 & ((reDinucleotide_Element_H9$E_posInf/2)/reDinucleotide_Element_H9$E_nCG)>=0.25)


# write.table(cbind(selected_reDinucleotide_Element_IMR90_IPSC$id_In_RE,selected_reDinucleotide_Element_IMR90_IPSC$meanMethCoef_IMR90_IPSC,selected_reDinucleotide_Element_IMR90_IPSC$id_In_Type,selected_reDinucleotide_Element_IMR90_IPSC$E_methCoef),file="reducedDataIMR90_IPSC_Filtered.txt",sep="\t",col.names=FALSE,row.names=F)
# write.table(cbind(selected_reDinucleotide_Element_H9$id_In_RE,selected_reDinucleotide_Element_H9$meanMethCoef_H9,selected_reDinucleotide_Element_H9$id_In_Type,selected_reDinucleotide_Element_H9$E_methCoef),file="reducedDataH9_Filtered.txt",sep="\t",col.names=FALSE,row.names=F)
# write.table(cbind(selected_reDinucleotide_Element_H9$id_In_RE,selected_reDinucleotide_Element_H9$meanMethCoef_H9,selected_reDinucleotide_Element_H9$id_In_Type,selected_reDinucleotide_Element_H9$E_methCoef),file="reducedDataH9_Filtered.txt",sep="\t",col.names=FALSE,row.names=F)

##########################
# Distance to extreme
t<-cbind(x=selected_reDinucleotide_Element_IMR90_IPSC$RE_position-selected_reDinucleotide_Element_IMR90_IPSC$E_chromStart,y=selected_reDinucleotide_Element_IMR90_IPSC$E_chromEnd-selected_reDinucleotide_Element_IMR90_IPSC$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_IMR90_IPSC<-cbind(selected_reDinucleotide_Element_IMR90_IPSC,dist2Ext=m)

t<-cbind(x=selected_reDinucleotide_Element_H9$RE_position-selected_reDinucleotide_Element_H9$E_chromStart,y=selected_reDinucleotide_Element_H9$E_chromEnd-selected_reDinucleotide_Element_H9$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_H9<-cbind(selected_reDinucleotide_Element_H9,dist2Ext=m)

t<-cbind(x=selected_reDinucleotide_Element_H9$RE_position-selected_reDinucleotide_Element_H9$E_chromStart,y=selected_reDinucleotide_Element_H9$E_chromEnd-selected_reDinucleotide_Element_H9$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_H9<-cbind(selected_reDinucleotide_Element_H9,dist2Ext=m)


#Obtain O/E ratio of the element

require(RMySQL)
con<-dbConnect(MySQL(),dbname="SolexaDB")
query<-"select id_cpgi_pk,obsExp FROM CpGisland"
rs1<-dbSendQuery(con,query)
OEdata<-fetch(rs1,n=-1) 
colnames(OEdata)=c("id_In_Type","E_OE")

selected_reDinucleotide_Element_IMR90_IPSC<-merge(selected_reDinucleotide_Element_IMR90_IPSC,OEdata,group="id_In_Type")
selected_reDinucleotide_Element_H9<-merge(selected_reDinucleotide_Element_H9,OEdata,group="id_In_Type")
selected_reDinucleotide_Element_H9<-merge(selected_reDinucleotide_Element_H9,OEdata,group="id_In_Type")

selected_reDinucleotide_Element_IMR90_IPSC<-selected_reDinucleotide_Element_IMR90_IPSC[,c(2:7,1,8:17)]
selected_reDinucleotide_Element_H9<-selected_reDinucleotide_Element_H9[,c(2:7,1,8:17)]
selected_reDinucleotide_Element_H9<-selected_reDinucleotide_Element_H9[,c(2:7,1,8:17)]

#Obtain length of the element

selected_reDinucleotide_Element_IMR90_IPSC<-cbind(selected_reDinucleotide_Element_IMR90_IPSC,
                                          cgiLength=selected_reDinucleotide_Element_IMR90_IPSC$E_chromEnd-
                                            selected_reDinucleotide_Element_IMR90_IPSC$E_chromStart)

selected_reDinucleotide_Element_H9<-cbind(selected_reDinucleotide_Element_H9,
                                             cgiLength=selected_reDinucleotide_Element_H9$E_chromEnd-
                                               selected_reDinucleotide_Element_H9$E_chromStart)
selected_reDinucleotide_Element_H9<-cbind(selected_reDinucleotide_Element_H9,
                                                   cgiLength=selected_reDinucleotide_Element_H9$E_chromEnd-
                                                     selected_reDinucleotide_Element_H9$E_chromStart)

############
#IMR90_IPSC

#lmIMR90_IPSCFiltered<-lm(selected_reDinucleotide_Element_IMR90_IPSC$E_methCoef~selected_reDinucleotide_Element_IMR90_IPSC$meanMethCoef_IMR90_IPSC)
#summary(lmIMR90_IPSCFiltered)
selected_reDinucleotide_Element_IMR90_IPSCJitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_IMR90_IPSC$E_methCoef,factor=25),meanMethCoef_IMR90_IPSC=jitter(selected_reDinucleotide_Element_IMR90_IPSC$meanMethCoef_IMR90_IPSC,factor=25)))

  
incorrectPoints_IMR90_IPSC_plus<-subset(selected_reDinucleotide_Element_IMR90_IPSC,selected_reDinucleotide_Element_IMR90_IPSC$E_methCoef-selected_reDinucleotide_Element_IMR90_IPSC$meanMethCoef_IMR90_IPSC>0.25)
#x_IMR90_IPSC_plus<-subset(selected_reDinucleotide_Element_IMR90_IPSCJitter,selected_reDinucleotide_Element_IMR90_IPSCJitter$E_methCoef-selected_reDinucleotide_Element_IMR90_IPSCJitter$meanMethCoef_IMR90_IPSC>0.25)
dim(incorrectPoints_IMR90_IPSC_plus)
x_IMR90_IPSC_plus<-subset(selected_reDinucleotide_Element_IMR90_IPSCJitter,selected_reDinucleotide_Element_IMR90_IPSCJitter$E_methCoef>=0.75)

incorrectPoints_IMR90_IPSC_minus<-subset(selected_reDinucleotide_Element_IMR90_IPSC,selected_reDinucleotide_Element_IMR90_IPSC$E_methCoef-selected_reDinucleotide_Element_IMR90_IPSC$meanMethCoef_IMR90_IPSC<(-0.25))
#x_IMR90_IPSC_minus<-subset(selected_reDinucleotide_Element_IMR90_IPSCJitter,selected_reDinucleotide_Element_IMR90_IPSCJitter$E_methCoef-selected_reDinucleotide_Element_IMR90_IPSCJitter$meanMethCoef_IMR90_IPSC<(-0.25))
dim(incorrectPoints_IMR90_IPSC_minus)
x_IMR90_IPSC_minus<-subset(selected_reDinucleotide_Element_IMR90_IPSCJitter,selected_reDinucleotide_Element_IMR90_IPSCJitter$E_methCoef<=0.25)

correctPoints_IMR90_IPSC<-subset(selected_reDinucleotide_Element_IMR90_IPSC,abs(selected_reDinucleotide_Element_IMR90_IPSC$E_methCoef-selected_reDinucleotide_Element_IMR90_IPSC$meanMethCoef_IMR90_IPSC)<=0.25)
# y_IMR90_IPSC<-subset(selected_reDinucleotide_Element_IMR90_IPSCJitter,abs(selected_reDinucleotide_Element_IMR90_IPSCJitter$E_methCoef-selected_reDinucleotide_Element_IMR90_IPSCJitter$meanMethCoef_IMR90_IPSC)<=0.25)
dim(correctPoints_IMR90_IPSC)
y_IMR90_IPSC<-subset(selected_reDinucleotide_Element_IMR90_IPSCJitter,selected_reDinucleotide_Element_IMR90_IPSCJitter$E_methCoef>0.25 & selected_reDinucleotide_Element_IMR90_IPSCJitter<0.75)

if (numbered=='TRUE'){
  png(paste(resultsDIR,"figure28_IMR90_IPSCFilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
} else {
  png(paste(resultsDIR,"figure28_IMR90_IPSCFilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_IMR90_IPSC$meanMethCoef_IMR90_IPSC, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_IMR90_IPSC$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_IMR90_IPSC$meanMethCoef_IMR90_IPSC,selected_reDinucleotide_Element_IMR90_IPSC$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n",xlim=c(0,1),ylim=c(0,1))

## Plot the points with some noise (jitter)
points(x_IMR90_IPSC_plus$meanMethCoef_IMR90_IPSC,x_IMR90_IPSC_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_IMR90_IPSC_minus$meanMethCoef_IMR90_IPSC,x_IMR90_IPSC_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_IMR90_IPSC$meanMethCoef_IMR90_IPSC,y_IMR90_IPSC$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# Plot the points without noise
#points(incorrectPoints_IMR90_IPSC_plus$meanMethCoef_IMR90_IPSC,incorrectPoints_IMR90_IPSC_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_IMR90_IPSC_minus$meanMethCoef_IMR90_IPSC,incorrectPoints_IMR90_IPSC_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_IMR90_IPSC$meanMethCoef_IMR90_IPSC,correctPoints_IMR90_IPSC$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if (numbered=='TRUE'){
# Adding the number of points to the plot
text(0.2,0.8,dim(incorrectPoints_IMR90_IPSC_plus)[1],cex=1.2,col="red")
text(0.8,0.2,dim(incorrectPoints_IMR90_IPSC_minus)[1],cex=1.2,col="red")
text(0.5,0.5,dim(correctPoints_IMR90_IPSC)[1],cex=1.2,col="red")
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


# H9

# lmH9Filtered<-lm(selected_reDinucleotide_Element_H9$E_methCoef~selected_reDinucleotide_Element_H9$meanMethCoef_H9)
# summary(lmH9Filtered)

selected_reDinucleotide_Element_H9Jitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_H9$E_methCoef,factor=25),meanMethCoef_H9=jitter(selected_reDinucleotide_Element_H9$meanMethCoef_H9,factor=25)))

incorrectPoints_H9_plus<-subset(selected_reDinucleotide_Element_H9,selected_reDinucleotide_Element_H9$E_methCoef-selected_reDinucleotide_Element_H9$meanMethCoef_H9>0.25)
#x_H9_plus<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef-selected_reDinucleotide_Element_H9Jitter$meanMethCoef_H9>0.25)
dim(incorrectPoints_H9_plus)
x_H9_plus<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef>=0.75)

incorrectPoints_H9_minus<-subset(selected_reDinucleotide_Element_H9,selected_reDinucleotide_Element_H9$E_methCoef-selected_reDinucleotide_Element_H9$meanMethCoef_H9<(-0.25))
#x_H9_minus<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef-selected_reDinucleotide_Element_H9Jitter$meanMethCoef_H9<(-0.25))
dim(incorrectPoints_H9_minus)
x_H9_minus<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef<=0.25)

correctPoints_H9<-subset(selected_reDinucleotide_Element_H9,abs(selected_reDinucleotide_Element_H9$E_methCoef-selected_reDinucleotide_Element_H9$meanMethCoef_H9)<=0.25)
#y_H9<-subset(selected_reDinucleotide_Element_H9Jitter,abs(selected_reDinucleotide_Element_H9Jitter$E_methCoef-selected_reDinucleotide_Element_H9Jitter$meanMethCoef_H9)<=0.25)
dim(correctPoints_H9)
y_H9<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef>0.25 & selected_reDinucleotide_Element_H9Jitter<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure28_H9FilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
  png(paste(resultsDIR,"figure28_H9FilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_H9$meanMethCoef_H9, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_H9$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_H9$meanMethCoef_H9,selected_reDinucleotide_Element_H9$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n",xlim=c(0,1),ylim=c(0,1))

## Plot the points with some noise (jitter)
points(x_H9_plus$meanMethCoef_H9,x_H9_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_H9_minus$meanMethCoef_H9,x_H9_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_H9$meanMethCoef_H9,y_H9$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# Plot the point without noise
#points(incorrectPoints_H9_plus$meanMethCoef_H9,incorrectPoints_H9_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_H9_minus$meanMethCoef_H9,incorrectPoints_H9_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_H9$meanMethCoef_H9,correctPoints_H9$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if(numbered=='TRUE'){
# Adding the number of points to the plot
text(0.2,0.8,dim(incorrectPoints_H9_plus)[1],cex=1.2,col="red")
text(0.8,0.2,dim(incorrectPoints_H9_minus)[1],cex=1.2,col="red")
text(0.5,0.5,dim(correctPoints_H9)[1],cex=1.2,col="red")
}else{}
abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

groupCPD_H9<-subset(correctPoints_H9,correctPoints_H9$E_methCoef<=0.25)
groupCPM_H9<-subset(correctPoints_H9,correctPoints_H9$E_methCoef>0.25 & correctPoints_H9$E_methCoef<0.75)
groupCPU_H9<-subset(correctPoints_H9,correctPoints_H9$E_methCoef>=0.75)

par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()

# H9
# lmH9Filtered<-lm(selected_reDinucleotide_Element_H9$E_methCoef~selected_reDinucleotide_Element_H9$meanMethCoef_H9)
# summary(lmH9Filtered)

selected_reDinucleotide_Element_H9Jitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_H9$E_methCoef,factor=25),meanMethCoef_H9=jitter(selected_reDinucleotide_Element_H9$meanMethCoef_H9,factor=25)))

incorrectPoints_H9_plus<-subset(selected_reDinucleotide_Element_H9,selected_reDinucleotide_Element_H9$E_methCoef-selected_reDinucleotide_Element_H9$meanMethCoef_H9>0.25)
#x_H9_plus<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef-selected_reDinucleotide_Element_H9Jitter$meanMethCoef_H9>0.25)
dim(incorrectPoints_H9_plus)
x_H9_plus<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef>=0.75)

incorrectPoints_H9_minus<-subset(selected_reDinucleotide_Element_H9,selected_reDinucleotide_Element_H9$E_methCoef-selected_reDinucleotide_Element_H9$meanMethCoef_H9<(-0.25))
#x_H9_minus<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef-selected_reDinucleotide_Element_H9Jitter$meanMethCoef_H9<(-0.25))
dim(incorrectPoints_H9_minus)
x_H9_minus<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef<=0.25)

correctPoints_H9<-subset(selected_reDinucleotide_Element_H9,abs(selected_reDinucleotide_Element_H9$E_methCoef-selected_reDinucleotide_Element_H9$meanMethCoef_H9)<=0.25)
#y_H9<-subset(selected_reDinucleotide_Element_H9Jitter,abs(selected_reDinucleotide_Element_H9Jitter$E_methCoef-selected_reDinucleotide_Element_H9Jitter$meanMethCoef_H9)<=0.25)
dim(correctPoints_H9)
y_H9<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef>0.25 & selected_reDinucleotide_Element_H9Jitter<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure28_H9FilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
  png(paste(resultsDIR,"figure28_H9FilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_H9$meanMethCoef_H9, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_H9$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_H9$meanMethCoef_H9,selected_reDinucleotide_Element_H9$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n",xlim=c(0,1),ylim=c(0,1))

## Plot the points with some noise (jitter)
points(x_H9_plus$meanMethCoef_H9,x_H9_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_H9_minus$meanMethCoef_H9,x_H9_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_H9$meanMethCoef_H9,y_H9$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# Plot the point without noise
#points(incorrectPoints_H9_plus$meanMethCoef_H9,incorrectPoints_H9_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_H9_minus$meanMethCoef_H9,incorrectPoints_H9_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_H9$meanMethCoef_H9,correctPoints_H9$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if(numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_H9_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_H9_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_H9)[1],cex=1.2,col="red")
}else{}
abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

groupCPD_H9<-subset(correctPoints_H9,correctPoints_H9$E_methCoef<=0.25)
groupCPM_H9<-subset(correctPoints_H9,correctPoints_H9$E_methCoef>0.25 & correctPoints_H9$E_methCoef<0.75)
groupCPU_H9<-subset(correctPoints_H9,correctPoints_H9$E_methCoef>=0.75)

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

selected_reDinucleotide_Element_IMR90_IPSC<-subset(reDinucleotide_Element_IMR90_IPSC,
                                           reDinucleotide_Element_IMR90_IPSC$meanMethCoef_IMR90_IPSC!=-1 & 
                                             reDinucleotide_Element_IMR90_IPSC$E_methCoef!=-1)
selected_reDinucleotide_Element_H9<-subset(reDinucleotide_Element_H9,
                                              reDinucleotide_Element_H9$meanMethCoef_H9!=-1 & 
                                                reDinucleotide_Element_H9$E_methCoef!=-1)
selected_reDinucleotide_Element_H9<-subset(reDinucleotide_Element_H9,
                                                    reDinucleotide_Element_H9$meanMethCoef_H9!=-1 & 
                                                      reDinucleotide_Element_H9$E_methCoef!=-1)

# write.table(cbind(selected_reDinucleotide_Element_IMR90_IPSC$id_In_RE,selected_reDinucleotide_Element_IMR90_IPSC$meanMethCoef_IMR90_IPSC,selected_reDinucleotide_Element_IMR90_IPSC$id_In_Type,selected_reDinucleotide_Element_IMR90_IPSC$E_methCoef),file="reducedDataIMR90_IPSC_All.txt",sep="\t",col.names=FALSE,row.names=F)
# write.table(cbind(selected_reDinucleotide_Element_H9$id_In_RE,selected_reDinucleotide_Element_H9$meanMethCoef_H9,selected_reDinucleotide_Element_H9$id_In_Type,selected_reDinucleotide_Element_H9$E_methCoef),file="reducedDataH9_All.txt",sep="\t",col.names=FALSE,row.names=F)
# write.table(cbind(selected_reDinucleotide_Element_H9$id_In_RE,selected_reDinucleotide_Element_H9$meanMethCoef_H9,selected_reDinucleotide_Element_H9$id_In_Type,selected_reDinucleotide_Element_H9$E_methCoef),file="reducedDataH9_All.txt",sep="\t",col.names=FALSE,row.names=F)



#######################
# Distance to extreme
t<-cbind(x=selected_reDinucleotide_Element_IMR90_IPSC$RE_position-selected_reDinucleotide_Element_IMR90_IPSC$E_chromStart,y=selected_reDinucleotide_Element_IMR90_IPSC$E_chromEnd-selected_reDinucleotide_Element_IMR90_IPSC$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_IMR90_IPSC<-cbind(selected_reDinucleotide_Element_IMR90_IPSC,dist2Ext=m)

t<-cbind(x=selected_reDinucleotide_Element_H9$RE_position-selected_reDinucleotide_Element_H9$E_chromStart,y=selected_reDinucleotide_Element_H9$E_chromEnd-selected_reDinucleotide_Element_H9$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_H9<-cbind(selected_reDinucleotide_Element_H9,dist2Ext=m)


t<-cbind(x=selected_reDinucleotide_Element_H9$RE_position-selected_reDinucleotide_Element_H9$E_chromStart,y=selected_reDinucleotide_Element_H9$E_chromEnd-selected_reDinucleotide_Element_H9$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_H9<-cbind(selected_reDinucleotide_Element_H9,dist2Ext=m)

#Obtain O/E ratio of the element

require(RMySQL)
con<-dbConnect(MySQL(),dbname="SolexaDB")
query<-"select id_cpgi_pk,obsExp FROM CpGisland"
rs1<-dbSendQuery(con,query)
OEdata<-fetch(rs1,n=-1) 
colnames(OEdata)=c("id_In_Type","E_OE")

selected_reDinucleotide_Element_IMR90_IPSC<-merge(selected_reDinucleotide_Element_IMR90_IPSC,OEdata,group="id_In_Type")
selected_reDinucleotide_Element_H9<-merge(selected_reDinucleotide_Element_H9,OEdata,group="id_In_Type")
selected_reDinucleotide_Element_H9<-merge(selected_reDinucleotide_Element_H9,OEdata,group="id_In_Type")


selected_reDinucleotide_Element_IMR90_IPSC<-selected_reDinucleotide_Element_IMR90_IPSC[,c(2:7,1,8:17)]
selected_reDinucleotide_Element_H9<-selected_reDinucleotide_Element_H9[,c(2:7,1,8:17)]
selected_reDinucleotide_Element_H9<-selected_reDinucleotide_Element_H9[,c(2:7,1,8:17)]

#Obtain length of the element

selected_reDinucleotide_Element_IMR90_IPSC<-cbind(selected_reDinucleotide_Element_IMR90_IPSC,
                                          cgiLength=selected_reDinucleotide_Element_IMR90_IPSC$E_chromEnd-
                                            selected_reDinucleotide_Element_IMR90_IPSC$E_chromStart)

selected_reDinucleotide_Element_H9<-cbind(selected_reDinucleotide_Element_H9,
                                             cgiLength=selected_reDinucleotide_Element_H9$E_chromEnd-
                                               selected_reDinucleotide_Element_H9$E_chromStart)
selected_reDinucleotide_Element_H9<-cbind(selected_reDinucleotide_Element_H9,
                                                   cgiLength=selected_reDinucleotide_Element_H9$E_chromEnd-
                                                     selected_reDinucleotide_Element_H9$E_chromStart)

############
#IMR90_IPSC
# lmIMR90_IPSCAll<-lm(selected_reDinucleotide_Element_IMR90_IPSC$E_methCoef~selected_reDinucleotide_Element_IMR90_IPSC$meanMethCoef_IMR90_IPSC)
# summary(lmIMR90_IPSCAll)
selected_reDinucleotide_Element_IMR90_IPSCJitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_IMR90_IPSC$E_methCoef,factor=25),meanMethCoef_IMR90_IPSC=jitter(selected_reDinucleotide_Element_IMR90_IPSC$meanMethCoef_IMR90_IPSC,factor=25)))

incorrectPoints_IMR90_IPSC_plus<-subset(selected_reDinucleotide_Element_IMR90_IPSC,selected_reDinucleotide_Element_IMR90_IPSC$E_methCoef-selected_reDinucleotide_Element_IMR90_IPSC$meanMethCoef_IMR90_IPSC>0.25)
#x_IMR90_IPSC_plus<-subset(selected_reDinucleotide_Element_IMR90_IPSCJitter,selected_reDinucleotide_Element_IMR90_IPSCJitter$E_methCoef-selected_reDinucleotide_Element_IMR90_IPSCJitter$meanMethCoef_IMR90_IPSC>0.25)
dim(incorrectPoints_IMR90_IPSC_plus)
x_IMR90_IPSC_plus<-subset(selected_reDinucleotide_Element_IMR90_IPSCJitter,selected_reDinucleotide_Element_IMR90_IPSCJitter$E_methCoef>=0.75)

incorrectPoints_IMR90_IPSC_minus<-subset(selected_reDinucleotide_Element_IMR90_IPSC,selected_reDinucleotide_Element_IMR90_IPSC$E_methCoef-selected_reDinucleotide_Element_IMR90_IPSC$meanMethCoef_IMR90_IPSC<(-0.25))
#x_IMR90_IPSC_minus<-subset(selected_reDinucleotide_Element_IMR90_IPSCJitter,selected_reDinucleotide_Element_IMR90_IPSCJitter$E_methCoef-selected_reDinucleotide_Element_IMR90_IPSCJitter$meanMethCoef_IMR90_IPSC<(-0.25))
dim(incorrectPoints_IMR90_IPSC_minus)
x_IMR90_IPSC_minus<-subset(selected_reDinucleotide_Element_IMR90_IPSCJitter,selected_reDinucleotide_Element_IMR90_IPSCJitter$E_methCoef<=0.25)

correctPoints_IMR90_IPSC<-subset(selected_reDinucleotide_Element_IMR90_IPSC,abs(selected_reDinucleotide_Element_IMR90_IPSC$E_methCoef-selected_reDinucleotide_Element_IMR90_IPSC$meanMethCoef_IMR90_IPSC)<=0.25)
# y_IMR90_IPSC<-subset(selected_reDinucleotide_Element_IMR90_IPSCJitter,abs(selected_reDinucleotide_Element_IMR90_IPSCJitter$E_methCoef-selected_reDinucleotide_Element_IMR90_IPSCJitter$meanMethCoef_IMR90_IPSC)<=0.25)
dim(correctPoints_IMR90_IPSC)
y_IMR90_IPSC<-subset(selected_reDinucleotide_Element_IMR90_IPSCJitter,selected_reDinucleotide_Element_IMR90_IPSCJitter$E_methCoef>0.25 & selected_reDinucleotide_Element_IMR90_IPSCJitter<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure28_IMR90_IPSCAllBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure28_IMR90_IPSCAllBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_IMR90_IPSC$meanMethCoef_IMR90_IPSC, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_IMR90_IPSC$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_IMR90_IPSC$meanMethCoef_IMR90_IPSC,selected_reDinucleotide_Element_IMR90_IPSC$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n",xlim=c(0,1),ylim=c(0,1))
points(x_IMR90_IPSC_plus$meanMethCoef_IMR90_IPSC,x_IMR90_IPSC_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_IMR90_IPSC_minus$meanMethCoef_IMR90_IPSC,x_IMR90_IPSC_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_IMR90_IPSC$meanMethCoef_IMR90_IPSC,y_IMR90_IPSC$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

#points(incorrectPoints_IMR90_IPSC_plus$meanMethCoef_IMR90_IPSC,incorrectPoints_IMR90_IPSC_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_IMR90_IPSC_minus$meanMethCoef_IMR90_IPSC,incorrectPoints_IMR90_IPSC_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_IMR90_IPSC$meanMethCoef_IMR90_IPSC,correctPoints_IMR90_IPSC$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if (numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_IMR90_IPSC_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_IMR90_IPSC_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_IMR90_IPSC)[1],cex=1.2,col="red")
}else{}

abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

groupCPD_IMR90_IPSC<-subset(correctPoints_IMR90_IPSC,correctPoints_IMR90_IPSC$E_methCoef<=0.25)
groupCPM_IMR90_IPSC<-subset(correctPoints_IMR90_IPSC,correctPoints_IMR90_IPSC$E_methCoef>0.25 & correctPoints_IMR90_IPSC$E_methCoef<0.75)
groupCPU_IMR90_IPSC<-subset(correctPoints_IMR90_IPSC,correctPoints_IMR90_IPSC$E_methCoef>=0.75)



par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()

# H9

# lmH9All<-lm(selected_reDinucleotide_Element_H9$E_methCoef~selected_reDinucleotide_Element_H9$meanMethCoef_H9)
# summary(lmH9All)

selected_reDinucleotide_Element_H9Jitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_H9$E_methCoef,factor=25),meanMethCoef_H9=jitter(selected_reDinucleotide_Element_H9$meanMethCoef_H9,factor=25)))

incorrectPoints_H9_plus<-subset(selected_reDinucleotide_Element_H9,selected_reDinucleotide_Element_H9$E_methCoef-selected_reDinucleotide_Element_H9$meanMethCoef_H9>0.25)
#x_H9_plus<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef-selected_reDinucleotide_Element_H9Jitter$meanMethCoef_H9>0.25)
dim(incorrectPoints_H9_plus)
x_H9_plus<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef>=0.75)

incorrectPoints_H9_minus<-subset(selected_reDinucleotide_Element_H9,selected_reDinucleotide_Element_H9$E_methCoef-selected_reDinucleotide_Element_H9$meanMethCoef_H9<(-0.25))
#x_H9_minus<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef-selected_reDinucleotide_Element_H9Jitter$meanMethCoef_H9<(-0.25))
dim(incorrectPoints_H9_minus)
x_H9_minus<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef<=0.25)

correctPoints_H9<-subset(selected_reDinucleotide_Element_H9,abs(selected_reDinucleotide_Element_H9$E_methCoef-selected_reDinucleotide_Element_H9$meanMethCoef_H9)<=0.25)
# y_H9<-subset(selected_reDinucleotide_Element_H9Jitter,abs(selected_reDinucleotide_Element_H9Jitter$E_methCoef-selected_reDinucleotide_Element_H9Jitter$meanMethCoef_H9)<=0.25)
dim(correctPoints_H9)
y_H9<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef>0.25 & selected_reDinucleotide_Element_H9Jitter<0.75)

if (numbered=='TRUE'){
  png(paste(resultsDIR,"figure28_AllH9BNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure28_AllH9BN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_H9$meanMethCoef_H9, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_H9$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_H9$meanMethCoef_H9,selected_reDinucleotide_Element_H9$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n",xlim=c(0,1),ylim=c(0,1))
points(x_H9_plus$meanMethCoef_H9,x_H9_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_H9_minus$meanMethCoef_H9,x_H9_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_H9$meanMethCoef_H9,y_H9$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

#points(incorrectPoints_H9_plus$meanMethCoef_H9,incorrectPoints_H9_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_H9_minus$meanMethCoef_H9,incorrectPoints_H9_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_H9$meanMethCoef_H9,correctPoints_H9$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if(numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_H9_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_H9_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_H9)[1],cex=1.2,col="red")
}else{}
abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

groupCPD_H9<-subset(correctPoints_H9,correctPoints_H9$E_methCoef<=0.25)
groupCPM_H9<-subset(correctPoints_H9,correctPoints_H9$E_methCoef>0.25 & correctPoints_H9$E_methCoef<0.75)
groupCPU_H9<-subset(correctPoints_H9,correctPoints_H9$E_methCoef>=0.75)

par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()


# H9
# lmH9All<-lm(selected_reDinucleotide_Element_H9$E_methCoef~selected_reDinucleotide_Element_H9$meanMethCoef_H9)
# summary(lmH9All)

selected_reDinucleotide_Element_H9Jitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_H9$E_methCoef,factor=25),meanMethCoef_H9=jitter(selected_reDinucleotide_Element_H9$meanMethCoef_H9,factor=25)))

incorrectPoints_H9_plus<-subset(selected_reDinucleotide_Element_H9,selected_reDinucleotide_Element_H9$E_methCoef-selected_reDinucleotide_Element_H9$meanMethCoef_H9>0.25)
#x_H9_plus<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef-selected_reDinucleotide_Element_H9Jitter$meanMethCoef_H9>0.25)
dim(incorrectPoints_H9_plus)
x_H9_plus<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef>=0.75)

incorrectPoints_H9_minus<-subset(selected_reDinucleotide_Element_H9,selected_reDinucleotide_Element_H9$E_methCoef-selected_reDinucleotide_Element_H9$meanMethCoef_H9<(-0.25))
#x_H9_minus<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef-selected_reDinucleotide_Element_H9Jitter$meanMethCoef_H9<(-0.25))
dim(incorrectPoints_H9_minus)
x_H9_minus<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef<=0.25)

correctPoints_H9<-subset(selected_reDinucleotide_Element_H9,abs(selected_reDinucleotide_Element_H9$E_methCoef-selected_reDinucleotide_Element_H9$meanMethCoef_H9)<=0.25)
# y_H9<-subset(selected_reDinucleotide_Element_H9Jitter,abs(selected_reDinucleotide_Element_H9Jitter$E_methCoef-selected_reDinucleotide_Element_H9Jitter$meanMethCoef_H9)<=0.25)
dim(correctPoints_H9)
y_H9<-subset(selected_reDinucleotide_Element_H9Jitter,selected_reDinucleotide_Element_H9Jitter$E_methCoef>0.25 & selected_reDinucleotide_Element_H9Jitter<0.75)

if (numbered=='TRUE'){
  png(paste(resultsDIR,"figure28_AllH9BNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
  png(paste(resultsDIR,"figure28_AllH9BN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_H9$meanMethCoef_H9, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_H9$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_H9$meanMethCoef_H9,selected_reDinucleotide_Element_H9$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n",xlim=c(0,1),ylim=c(0,1))
points(x_H9_plus$meanMethCoef_H9,x_H9_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_H9_minus$meanMethCoef_H9,x_H9_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_H9$meanMethCoef_H9,y_H9$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

#points(incorrectPoints_H9_plus$meanMethCoef_H9,incorrectPoints_H9_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_H9_minus$meanMethCoef_H9,incorrectPoints_H9_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_H9$meanMethCoef_H9,correctPoints_H9$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if(numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_H9_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_H9_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_H9)[1],cex=1.2,col="red")
}else{}
abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

groupCPD_H9<-subset(correctPoints_H9,correctPoints_H9$E_methCoef<=0.25)
groupCPM_H9<-subset(correctPoints_H9,correctPoints_H9$E_methCoef>0.25 & correctPoints_H9$E_methCoef<0.75)
groupCPU_H9<-subset(correctPoints_H9,correctPoints_H9$E_methCoef>=0.75)

par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()


