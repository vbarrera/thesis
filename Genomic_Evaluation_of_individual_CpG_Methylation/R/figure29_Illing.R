##################
### Figure 29 ###
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
WHERE MPA.name=MEA.name AND R.RE_name='HpaII' AND E.type='CpGislandIllingNonOver' 
AND MPA.RE_name='HpaII' AND MEA.type='CpGislandIllingNonOver' AND C.type='CpGislandIllingNonOver'")

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
# Distance to extreme
t<-cbind(x=selected_reDinucleotide_Element_H1$RE_position-selected_reDinucleotide_Element_H1$E_chromStart,y=selected_reDinucleotide_Element_H1$E_chromEnd-selected_reDinucleotide_Element_H1$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_H1<-cbind(selected_reDinucleotide_Element_H1,dist2Ext=m)

t<-cbind(x=selected_reDinucleotide_Element_IMR90$RE_position-selected_reDinucleotide_Element_IMR90$E_chromStart,y=selected_reDinucleotide_Element_IMR90$E_chromEnd-selected_reDinucleotide_Element_IMR90$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_IMR90<-cbind(selected_reDinucleotide_Element_IMR90,dist2Ext=m)


#Obtain O/E ratio of the element

require(RMySQL)
con<-dbConnect(MySQL(),dbname="SolexaDB")
query<-"select id_cpgi_pk,obsExp FROM CpGisland"
rs1<-dbSendQuery(con,query)
OEdata<-fetch(rs1,n=-1) 
colnames(OEdata)=c("id_In_Type","E_OE")

selected_reDinucleotide_Element_H1<-merge(selected_reDinucleotide_Element_H1,OEdata,group="id_In_Type")
selected_reDinucleotide_Element_IMR90<-merge(selected_reDinucleotide_Element_IMR90,OEdata,group="id_In_Type")

selected_reDinucleotide_Element_H1<-selected_reDinucleotide_Element_H1[,c(2:7,1,8:17)]
selected_reDinucleotide_Element_IMR90<-selected_reDinucleotide_Element_IMR90[,c(2:7,1,8:17)]

#Obtain length of the element

selected_reDinucleotide_Element_H1<-cbind(selected_reDinucleotide_Element_H1,
                                          cgiLength=selected_reDinucleotide_Element_H1$E_chromEnd-
                                            selected_reDinucleotide_Element_H1$E_chromStart)

selected_reDinucleotide_Element_IMR90<-cbind(selected_reDinucleotide_Element_IMR90,
                                             cgiLength=selected_reDinucleotide_Element_IMR90$E_chromEnd-
                                               selected_reDinucleotide_Element_IMR90$E_chromStart)

############
#H1
# lmH1Filtered<-lm(selected_reDinucleotide_Element_H1$E_methCoef~selected_reDinucleotide_Element_H1$meanMethCoef_H1)
# summary(lmH1Filtered)


selected_reDinucleotide_Element_H1Jitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_H1$E_methCoef,factor=25),meanMethCoef_H1=jitter(selected_reDinucleotide_Element_H1$meanMethCoef_H1,factor=25)))

incorrectPoints_H1_plus<-subset(selected_reDinucleotide_Element_H1,selected_reDinucleotide_Element_H1$E_methCoef-selected_reDinucleotide_Element_H1$meanMethCoef_H1>0.25)
#x_H1_plus<-subset(selected_reDinucleotide_Element_H1Jitter,selected_reDinucleotide_Element_H1Jitter$E_methCoef-selected_reDinucleotide_Element_H1Jitter$meanMethCoef_H1>0.25)
dim(incorrectPoints_H1_plus)
x_H1_plus<-subset(selected_reDinucleotide_Element_H1Jitter,selected_reDinucleotide_Element_H1Jitter$E_methCoef>=0.75)

incorrectPoints_H1_minus<-subset(selected_reDinucleotide_Element_H1,selected_reDinucleotide_Element_H1$E_methCoef-selected_reDinucleotide_Element_H1$meanMethCoef_H1<(-0.25))
#x_H1_minus<-subset(selected_reDinucleotide_Element_H1Jitter,selected_reDinucleotide_Element_H1Jitter$E_methCoef-selected_reDinucleotide_Element_H1Jitter$meanMethCoef_H1<(-0.25))
dim(incorrectPoints_H1_minus)
x_H1_minus<-subset(selected_reDinucleotide_Element_H1Jitter,selected_reDinucleotide_Element_H1Jitter$E_methCoef<=0.25)

correctPoints_H1<-subset(selected_reDinucleotide_Element_H1,abs(selected_reDinucleotide_Element_H1$E_methCoef-selected_reDinucleotide_Element_H1$meanMethCoef_H1)<=0.25)
# y_H1<-subset(selected_reDinucleotide_Element_H1Jitter,abs(selected_reDinucleotide_Element_H1Jitter$E_methCoef-selected_reDinucleotide_Element_H1Jitter$meanMethCoef_H1)<=0.25)
dim(correctPoints_H1)
y_H1<-subset(selected_reDinucleotide_Element_H1Jitter,selected_reDinucleotide_Element_H1Jitter$E_methCoef>0.25 & selected_reDinucleotide_Element_H1Jitter<0.75)

if (numbered=='TRUE'){
  png(paste(resultsDIR,"figure29_IllingH1FilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
} else {
  png(paste(resultsDIR,"figure29_IllingH1FilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_H1$meanMethCoef_H1, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_H1$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_H1$meanMethCoef_H1,selected_reDinucleotide_Element_H1$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n")

## Plot the points with some noise (jitter)
points(x_H1_plus$meanMethCoef_H1,x_H1_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_H1_minus$meanMethCoef_H1,x_H1_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_H1$meanMethCoef_H1,y_H1$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# Plot the points without noise
#points(incorrectPoints_H1_plus$meanMethCoef_H1,incorrectPoints_H1_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_H1_minus$meanMethCoef_H1,incorrectPoints_H1_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_H1$meanMethCoef_H1,correctPoints_H1$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if (numbered=='TRUE'){
# Adding the number of points to the plot
text(0.2,0.8,dim(incorrectPoints_H1_plus)[1],cex=1.2,col="red")
text(0.8,0.2,dim(incorrectPoints_H1_minus)[1],cex=1.2,col="red")
text(0.5,0.5,dim(correctPoints_H1)[1],cex=1.2,col="red")
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


# IMR90
# lmIMR90Filtered<-lm(selected_reDinucleotide_Element_IMR90$E_methCoef~selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90)
# summary(lmIMR90Filtered)


selected_reDinucleotide_Element_IMR90Jitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_IMR90$E_methCoef,factor=25),meanMethCoef_IMR90=jitter(selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90,factor=25)))

incorrectPoints_IMR90_plus<-subset(selected_reDinucleotide_Element_IMR90,selected_reDinucleotide_Element_IMR90$E_methCoef-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90>0.25)
#x_IMR90_plus<-subset(selected_reDinucleotide_Element_IMR90Jitter,selected_reDinucleotide_Element_IMR90Jitter$E_methCoef-selected_reDinucleotide_Element_IMR90Jitter$meanMethCoef_IMR90>0.25)
dim(incorrectPoints_IMR90_plus)
x_IMR90_plus<-subset(selected_reDinucleotide_Element_IMR90Jitter,selected_reDinucleotide_Element_IMR90Jitter$E_methCoef>=0.75)

incorrectPoints_IMR90_minus<-subset(selected_reDinucleotide_Element_IMR90,selected_reDinucleotide_Element_IMR90$E_methCoef-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90<(-0.25))
#x_IMR90_minus<-subset(selected_reDinucleotide_Element_IMR90Jitter,selected_reDinucleotide_Element_IMR90Jitter$E_methCoef-selected_reDinucleotide_Element_IMR90Jitter$meanMethCoef_IMR90<(-0.25))
dim(incorrectPoints_IMR90_minus)
x_IMR90_minus<-subset(selected_reDinucleotide_Element_IMR90Jitter,selected_reDinucleotide_Element_IMR90Jitter$E_methCoef<=0.25)

correctPoints_IMR90<-subset(selected_reDinucleotide_Element_IMR90,abs(selected_reDinucleotide_Element_IMR90$E_methCoef-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90)<=0.25)
#y_IMR90<-subset(selected_reDinucleotide_Element_IMR90Jitter,abs(selected_reDinucleotide_Element_IMR90Jitter$E_methCoef-selected_reDinucleotide_Element_IMR90Jitter$meanMethCoef_IMR90)<=0.25)
dim(correctPoints_IMR90)
y_IMR90<-subset(selected_reDinucleotide_Element_IMR90Jitter,selected_reDinucleotide_Element_IMR90Jitter$E_methCoef>0.25 & selected_reDinucleotide_Element_IMR90Jitter<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure29_IllingIMR90FilteredBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
  png(paste(resultsDIR,"figure29_IllingIMR90FilteredBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_IMR90$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90,selected_reDinucleotide_Element_IMR90$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n")

## Plot the points with some noise (jitter)
points(x_IMR90_plus$meanMethCoef_IMR90,x_IMR90_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_IMR90_minus$meanMethCoef_IMR90,x_IMR90_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_IMR90$meanMethCoef_IMR90,y_IMR90$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

# Plot the point without noise
#points(incorrectPoints_IMR90_plus$meanMethCoef_IMR90,incorrectPoints_IMR90_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_IMR90_minus$meanMethCoef_IMR90,incorrectPoints_IMR90_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_IMR90$meanMethCoef_IMR90,correctPoints_IMR90$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if(numbered=='TRUE'){
# Adding the number of points to the plot
text(0.2,0.8,dim(incorrectPoints_IMR90_plus)[1],cex=1.2,col="red")
text(0.8,0.2,dim(incorrectPoints_IMR90_minus)[1],cex=1.2,col="red")
text(0.5,0.5,dim(correctPoints_IMR90)[1],cex=1.2,col="red")
}else{}
abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

groupCPD_IMR90<-subset(correctPoints_IMR90,correctPoints_IMR90$E_methCoef<=0.25)
groupCPM_IMR90<-subset(correctPoints_IMR90,correctPoints_IMR90$E_methCoef>0.25 & correctPoints_IMR90$E_methCoef<0.75)
groupCPU_IMR90<-subset(correctPoints_IMR90,correctPoints_IMR90$E_methCoef>=0.75)

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

selected_reDinucleotide_Element_H1<-subset(reDinucleotide_Element_H1,
                                           reDinucleotide_Element_H1$meanMethCoef_H1!=-1 & 
                                             reDinucleotide_Element_H1$E_methCoef!=-1)
selected_reDinucleotide_Element_IMR90<-subset(reDinucleotide_Element_IMR90,
                                              reDinucleotide_Element_IMR90$meanMethCoef_IMR90!=-1 & 
                                                reDinucleotide_Element_IMR90$E_methCoef!=-1)

#######################
# Distance to extreme
t<-cbind(x=selected_reDinucleotide_Element_H1$RE_position-selected_reDinucleotide_Element_H1$E_chromStart,y=selected_reDinucleotide_Element_H1$E_chromEnd-selected_reDinucleotide_Element_H1$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_H1<-cbind(selected_reDinucleotide_Element_H1,dist2Ext=m)

t<-cbind(x=selected_reDinucleotide_Element_IMR90$RE_position-selected_reDinucleotide_Element_IMR90$E_chromStart,y=selected_reDinucleotide_Element_IMR90$E_chromEnd-selected_reDinucleotide_Element_IMR90$RE_position)
m<-apply(t,1,min)
selected_reDinucleotide_Element_IMR90<-cbind(selected_reDinucleotide_Element_IMR90,dist2Ext=m)


#Obtain O/E ratio of the element

require(RMySQL)
con<-dbConnect(MySQL(),dbname="SolexaDB")
query<-"select id_cpgi_pk,obsExp FROM CpGisland"
rs1<-dbSendQuery(con,query)
OEdata<-fetch(rs1,n=-1) 
colnames(OEdata)=c("id_In_Type","E_OE")

selected_reDinucleotide_Element_H1<-merge(selected_reDinucleotide_Element_H1,OEdata,group="id_In_Type")
selected_reDinucleotide_Element_IMR90<-merge(selected_reDinucleotide_Element_IMR90,OEdata,group="id_In_Type")

selected_reDinucleotide_Element_H1<-selected_reDinucleotide_Element_H1[,c(2:7,1,8:17)]
selected_reDinucleotide_Element_IMR90<-selected_reDinucleotide_Element_IMR90[,c(2:7,1,8:17)]

#Obtain length of the element

selected_reDinucleotide_Element_H1<-cbind(selected_reDinucleotide_Element_H1,
                                          cgiLength=selected_reDinucleotide_Element_H1$E_chromEnd-
                                            selected_reDinucleotide_Element_H1$E_chromStart)

selected_reDinucleotide_Element_IMR90<-cbind(selected_reDinucleotide_Element_IMR90,
                                             cgiLength=selected_reDinucleotide_Element_IMR90$E_chromEnd-
                                               selected_reDinucleotide_Element_IMR90$E_chromStart)

############
#H1
# lmH1All<-lm(selected_reDinucleotide_Element_H1$E_methCoef~selected_reDinucleotide_Element_H1$meanMethCoef_H1)
# summary(lmH1All)


selected_reDinucleotide_Element_H1Jitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_H1$E_methCoef,factor=25),meanMethCoef_H1=jitter(selected_reDinucleotide_Element_H1$meanMethCoef_H1,factor=25)))

incorrectPoints_H1_plus<-subset(selected_reDinucleotide_Element_H1,selected_reDinucleotide_Element_H1$E_methCoef-selected_reDinucleotide_Element_H1$meanMethCoef_H1>0.25)
#x_H1_plus<-subset(selected_reDinucleotide_Element_H1Jitter,selected_reDinucleotide_Element_H1Jitter$E_methCoef-selected_reDinucleotide_Element_H1Jitter$meanMethCoef_H1>0.25)
dim(incorrectPoints_H1_plus)
x_H1_plus<-subset(selected_reDinucleotide_Element_H1Jitter,selected_reDinucleotide_Element_H1Jitter$E_methCoef>=0.75)

incorrectPoints_H1_minus<-subset(selected_reDinucleotide_Element_H1,selected_reDinucleotide_Element_H1$E_methCoef-selected_reDinucleotide_Element_H1$meanMethCoef_H1<(-0.25))
#x_H1_minus<-subset(selected_reDinucleotide_Element_H1Jitter,selected_reDinucleotide_Element_H1Jitter$E_methCoef-selected_reDinucleotide_Element_H1Jitter$meanMethCoef_H1<(-0.25))
dim(incorrectPoints_H1_minus)
x_H1_minus<-subset(selected_reDinucleotide_Element_H1Jitter,selected_reDinucleotide_Element_H1Jitter$E_methCoef<=0.25)

correctPoints_H1<-subset(selected_reDinucleotide_Element_H1,abs(selected_reDinucleotide_Element_H1$E_methCoef-selected_reDinucleotide_Element_H1$meanMethCoef_H1)<=0.25)
# y_H1<-subset(selected_reDinucleotide_Element_H1Jitter,abs(selected_reDinucleotide_Element_H1Jitter$E_methCoef-selected_reDinucleotide_Element_H1Jitter$meanMethCoef_H1)<=0.25)
dim(correctPoints_H1)
y_H1<-subset(selected_reDinucleotide_Element_H1Jitter,selected_reDinucleotide_Element_H1Jitter$E_methCoef>0.25 & selected_reDinucleotide_Element_H1Jitter<0.75)

if(numbered=='TRUE'){
  png(paste(resultsDIR,"figure29_IllingH1AllBNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure29_IllingH1AllBN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_H1$meanMethCoef_H1, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_H1$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_H1$meanMethCoef_H1,selected_reDinucleotide_Element_H1$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n")
points(x_H1_plus$meanMethCoef_H1,x_H1_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_H1_minus$meanMethCoef_H1,x_H1_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_H1$meanMethCoef_H1,y_H1$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

#points(incorrectPoints_H1_plus$meanMethCoef_H1,incorrectPoints_H1_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_H1_minus$meanMethCoef_H1,incorrectPoints_H1_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_H1$meanMethCoef_H1,correctPoints_H1$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if (numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_H1_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_H1_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_H1)[1],cex=1.2,col="red")
}else{}

abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

groupCPD_H1<-subset(correctPoints_H1,correctPoints_H1$E_methCoef<=0.25)
groupCPM_H1<-subset(correctPoints_H1,correctPoints_H1$E_methCoef>0.25 & correctPoints_H1$E_methCoef<0.75)
groupCPU_H1<-subset(correctPoints_H1,correctPoints_H1$E_methCoef>=0.75)



par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()

# IMR90
# lmIMR90All<-lm(selected_reDinucleotide_Element_IMR90$E_methCoef~selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90)
# summary(lmIMR90All)

selected_reDinucleotide_Element_IMR90Jitter<-as.data.frame(cbind(E_methCoef=jitter(selected_reDinucleotide_Element_IMR90$E_methCoef,factor=25),meanMethCoef_IMR90=jitter(selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90,factor=25)))

incorrectPoints_IMR90_plus<-subset(selected_reDinucleotide_Element_IMR90,selected_reDinucleotide_Element_IMR90$E_methCoef-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90>0.25)
#x_IMR90_plus<-subset(selected_reDinucleotide_Element_IMR90Jitter,selected_reDinucleotide_Element_IMR90Jitter$E_methCoef-selected_reDinucleotide_Element_IMR90Jitter$meanMethCoef_IMR90>0.25)
dim(incorrectPoints_IMR90_plus)
x_IMR90_plus<-subset(selected_reDinucleotide_Element_IMR90Jitter,selected_reDinucleotide_Element_IMR90Jitter$E_methCoef>=0.75)

incorrectPoints_IMR90_minus<-subset(selected_reDinucleotide_Element_IMR90,selected_reDinucleotide_Element_IMR90$E_methCoef-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90<(-0.25))
#x_IMR90_minus<-subset(selected_reDinucleotide_Element_IMR90Jitter,selected_reDinucleotide_Element_IMR90Jitter$E_methCoef-selected_reDinucleotide_Element_IMR90Jitter$meanMethCoef_IMR90<(-0.25))
dim(incorrectPoints_IMR90_minus)
x_IMR90_minus<-subset(selected_reDinucleotide_Element_IMR90Jitter,selected_reDinucleotide_Element_IMR90Jitter$E_methCoef<=0.25)

correctPoints_IMR90<-subset(selected_reDinucleotide_Element_IMR90,abs(selected_reDinucleotide_Element_IMR90$E_methCoef-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90)<=0.25)
# y_IMR90<-subset(selected_reDinucleotide_Element_IMR90Jitter,abs(selected_reDinucleotide_Element_IMR90Jitter$E_methCoef-selected_reDinucleotide_Element_IMR90Jitter$meanMethCoef_IMR90)<=0.25)
dim(correctPoints_IMR90)
y_IMR90<-subset(selected_reDinucleotide_Element_IMR90Jitter,selected_reDinucleotide_Element_IMR90Jitter$E_methCoef>0.25 & selected_reDinucleotide_Element_IMR90Jitter<0.75)

if (numbered=='TRUE'){
  png(paste(resultsDIR,"figure29_IllingAllIMR90BNNumbered.png",sep=""),height=15,width=15,units="cm",res=300)
}else{
png(paste(resultsDIR,"figure29_IllingAllIMR90BN.png",sep=""),height=15,width=15,units="cm",res=300)
}
def.par <- par(no.readonly = TRUE)
par(lwd=3)
par(cex.axis=1.2)
xhist <- hist(selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90, breaks=seq(0,1,0.1), plot=FALSE)
yhist <- hist(selected_reDinucleotide_Element_IMR90$E_methCoef, breaks=seq(0,1,0.1), plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(0,1)
yrange <- c(0,1)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)

par(mar=c(5,5,1,1))
plot (selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90,selected_reDinucleotide_Element_IMR90$E_methCoef,pch=16,col=rgb(1,0,0,0.2),xlab="",
      ylab="",main="",type="n")
points(x_IMR90_plus$meanMethCoef_IMR90,x_IMR90_plus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(x_IMR90_minus$meanMethCoef_IMR90,x_IMR90_minus$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)
points(y_IMR90$meanMethCoef_IMR90,y_IMR90$E_methCoef,pch=16,col=rgb(0,0,0,0.2),cex=0.8)

#points(incorrectPoints_IMR90_plus$meanMethCoef_IMR90,incorrectPoints_IMR90_plus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(incorrectPoints_IMR90_minus$meanMethCoef_IMR90,incorrectPoints_IMR90_minus$E_methCoef,pch=16,col=rgb(0.5,0,0.5,0.2),cex=1)
#points(correctPoints_IMR90$meanMethCoef_IMR90,correctPoints_IMR90$E_methCoef,pch=16,col=rgb(0,0.5,0.5,0.2),cex=1)

if(numbered=='TRUE'){
  # Adding the number of points to the plot
  text(0.2,0.8,dim(incorrectPoints_IMR90_plus)[1],cex=1.2,col="red")
  text(0.8,0.2,dim(incorrectPoints_IMR90_minus)[1],cex=1.2,col="red")
  text(0.5,0.5,dim(correctPoints_IMR90)[1],cex=1.2,col="red")
}else{}
abline(0.25,1,lty=2,lwd=1.5,col="red")
abline(-0.25,1,lty=2,lwd=1.5,col="red")

groupCPD_IMR90<-subset(correctPoints_IMR90,correctPoints_IMR90$E_methCoef<=0.25)
groupCPM_IMR90<-subset(correctPoints_IMR90,correctPoints_IMR90$E_methCoef>0.25 & correctPoints_IMR90$E_methCoef<0.75)
groupCPU_IMR90<-subset(correctPoints_IMR90,correctPoints_IMR90$E_methCoef>=0.75)

par(mar=c(0,5,1,1))
par(lwd=1)
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,
        xlab="SD")
par(mar=c(5,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0,
        horiz=TRUE)

par(def.par)
dev.off()
