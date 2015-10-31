#################
### Figure 33 ###
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



incorrectPoints_H1_plus<-subset(selected_reDinucleotide_Element_H1,selected_reDinucleotide_Element_H1$E_methCoef-selected_reDinucleotide_Element_H1$meanMethCoef_H1>0.25)
dim(incorrectPoints_H1_plus)

incorrectPoints_H1_minus<-subset(selected_reDinucleotide_Element_H1,selected_reDinucleotide_Element_H1$E_methCoef-selected_reDinucleotide_Element_H1$meanMethCoef_H1<(-0.25))
dim(incorrectPoints_H1_minus)


correctPoints_H1<-subset(selected_reDinucleotide_Element_H1,abs(selected_reDinucleotide_Element_H1$E_methCoef-selected_reDinucleotide_Element_H1$meanMethCoef_H1)<=0.25)
dim(correctPoints_H1)

incorrectPoints_IMR90_plus<-subset(selected_reDinucleotide_Element_IMR90,selected_reDinucleotide_Element_IMR90$E_methCoef-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90>0.25)
dim(incorrectPoints_IMR90_plus)

incorrectPoints_IMR90_minus<-subset(selected_reDinucleotide_Element_IMR90,selected_reDinucleotide_Element_IMR90$E_methCoef-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90<(-0.25))
dim(incorrectPoints_IMR90_minus)


correctPoints_IMR90<-subset(selected_reDinucleotide_Element_IMR90,abs(selected_reDinucleotide_Element_IMR90$E_methCoef-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90)<=0.25)
dim(correctPoints_IMR90)
# 2.2.4.1- Distance to extreme analysis.
#H1

png(paste(resultsDIR,"figure33_H1FilteredDist2Ext.png",sep=""),units="cm",height=8,width=8,res=300)
def.par <- par(no.readonly = TRUE)

par(lwd=1.5)
par(cex.axis=0.8)

plot(density(log2(correctPoints_H1$dist2Ext)),type="n",main="",xlab="",ylab="")
lines(density(log2(incorrectPoints_H1_minus$dist2Ext)),col="blue",lwd=0.8)
lines(density(log2(incorrectPoints_H1_plus$dist2Ext)),col="red",lwd=0.8)
lines(density(log2(correctPoints_H1$dist2Ext)),col="green",lwd=0.8)
dev.off()

#IMR90

png(paste(resultsDIR,"figure33_IMR90FilteredDist2Ext.png",sep=""),units="cm",height=8,width=8,res=300)
def.par <- par(no.readonly = TRUE)

par(lwd=1.5)
par(cex.axis=0.8)
plot(density(log2(correctPoints_IMR90$dist2Ext)),type="n",main="",xlab="",ylab="")
lines(density(log2(incorrectPoints_IMR90_minus$dist2Ext)),col="blue",lwd=0.8)
lines(density(log2(incorrectPoints_IMR90_plus$dist2Ext)),col="red",lwd=0.8)
lines(density(log2(correctPoints_IMR90$dist2Ext)),col="green",lwd=0.8)
dev.off()

# 2.2.4.2- CpG island length

#H1

png(paste(resultsDIR,"figure33_H1FilteredCpgiLength.png",sep=""),units="cm",height=8,width=8,res=300)
def.par <- par(no.readonly = TRUE)

par(lwd=1.5)
par(cex.axis=0.8)
plot(density(log2(correctPoints_H1$cgiLength)),type="n",main="",xlab="",ylab="")
lines(density(log2(incorrectPoints_H1_minus$cgiLength)),col="blue",lwd=0.8)
lines(density(log2(incorrectPoints_H1_plus$cgiLength)),col="red",lwd=0.8)
lines(density(log2(correctPoints_H1$cgiLength)),col="green",lwd=0.8)
dev.off()

#IMR90

png(paste(resultsDIR,"figure33_IMR90FilteredCpgiLength.png",sep=""),units="cm",height=8,width=8,res=300)
def.par <- par(no.readonly = TRUE)

par(lwd=1.5)
par(cex.axis=0.8)
plot(density(log2(correctPoints_IMR90$cgiLength)),type="n",main="",xlab="",ylab="")
lines(density(log2(incorrectPoints_IMR90_minus$cgiLength)),col="blue",lwd=0.8)
lines(density(log2(incorrectPoints_IMR90_plus$cgiLength)),col="red",lwd=0.8)
lines(density(log2(correctPoints_IMR90$cgiLength)),col="green",lwd=0.8)
dev.off()

# 2.2.4.3- CpG island nCpG

#H1
png(paste(resultsDIR,"figure33_H1FilteredCpginCG.png",sep=""),units="cm",height=8,width=8,res=300)
def.par <- par(no.readonly = TRUE)

par(lwd=1.5)
par(cex.axis=0.8)
plot(density(log2(correctPoints_H1$E_nCG)),type="n",main="",xlab="",ylab="")
lines(density(log2(incorrectPoints_H1_minus$E_nCG)),col="blue",lwd=0.8)
lines(density(log2(incorrectPoints_H1_plus$E_nCG)),col="red",lwd=0.8)
lines(density(log2(correctPoints_H1$E_nCG)),col="green",lwd=0.8)
dev.off()

#IMR90

png(paste(resultsDIR,"figure33_IMR90FilteredCpginCG.png",sep=""),units="cm",height=8,width=8,res=300)
def.par <- par(no.readonly = TRUE)

par(lwd=1.5)
par(cex.axis=0.8)
plot(density(log2(correctPoints_IMR90$E_nCG)),type="n",main="",xlab="",ylab="")
lines(density(log2(incorrectPoints_IMR90_minus$E_nCG)),col="blue",lwd=0.8)
lines(density(log2(incorrectPoints_IMR90_plus$E_nCG)),col="red",lwd=0.8)
lines(density(log2(correctPoints_IMR90$E_nCG)),col="green",lwd=0.8)
dev.off()

# 2.2.4.4- CpG island O/E

#H1
png(paste(resultsDIR,"figure33_H1FilteredCpgiOE.png",sep=""),units="cm",height=8,width=8,res=300)
def.par <- par(no.readonly = TRUE)

par(lwd=1.5)
par(cex.axis=0.8)
plot(density(log2(incorrectPoints_H1_minus$E_OE)),type="n",main="",xlab="",ylab="")
lines(density(log2(incorrectPoints_H1_minus$E_OE)),col="blue",lwd=0.8)
lines(density(log2(incorrectPoints_H1_plus$E_OE)),col="red",lwd=0.8)
lines(density(log2(correctPoints_H1$E_OE)),col="green",lwd=0.8)
dev.off()

#IMR90
png(paste(resultsDIR,"figure33_IMR90FilteredCpgiOE.png",sep=""),units="cm",height=8,width=8,res=300)
def.par <- par(no.readonly = TRUE)

par(lwd=1.5)
par(cex.axis=0.8)
plot(density(log2(incorrectPoints_IMR90_minus$E_OE)),type="n",main="",xlab="",ylab="")
lines(density(log2(incorrectPoints_IMR90_minus$E_OE)),col="blue",lwd=0.8)
lines(density(log2(incorrectPoints_IMR90_plus$E_OE)),col="red",lwd=0.8)
lines(density(log2(correctPoints_IMR90$E_OE)),col="green",lwd=0.8)
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

incorrectPoints_H1_plus<-subset(selected_reDinucleotide_Element_H1,selected_reDinucleotide_Element_H1$E_methCoef-selected_reDinucleotide_Element_H1$meanMethCoef_H1>0.25)
dim(incorrectPoints_H1_plus)

incorrectPoints_H1_minus<-subset(selected_reDinucleotide_Element_H1,selected_reDinucleotide_Element_H1$E_methCoef-selected_reDinucleotide_Element_H1$meanMethCoef_H1<(-0.25))
dim(incorrectPoints_H1_minus)


correctPoints_H1<-subset(selected_reDinucleotide_Element_H1,abs(selected_reDinucleotide_Element_H1$E_methCoef-selected_reDinucleotide_Element_H1$meanMethCoef_H1)<=0.25)
dim(correctPoints_H1)

incorrectPoints_IMR90_plus<-subset(selected_reDinucleotide_Element_IMR90,selected_reDinucleotide_Element_IMR90$E_methCoef-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90>0.25)
dim(incorrectPoints_IMR90_plus)

incorrectPoints_IMR90_minus<-subset(selected_reDinucleotide_Element_IMR90,selected_reDinucleotide_Element_IMR90$E_methCoef-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90<(-0.25))
dim(incorrectPoints_IMR90_minus)


correctPoints_IMR90<-subset(selected_reDinucleotide_Element_IMR90,abs(selected_reDinucleotide_Element_IMR90$E_methCoef-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90)<=0.25)
dim(correctPoints_IMR90)
# 2.2.4.1- Distance to extreme analysis.
#H1

png(paste(resultsDIR,"figure33_H1AllDist2Ext.png",sep=""),units="cm",height=8,width=8,res=300)
def.par <- par(no.readonly = TRUE)

par(lwd=1.5)
par(cex.axis=0.8)

plot(density(log2(correctPoints_H1$dist2Ext)),type="n",main="",xlab="",ylab="")
lines(density(log2(incorrectPoints_H1_minus$dist2Ext)),col="blue",lwd=0.8)
lines(density(log2(incorrectPoints_H1_plus$dist2Ext)),col="red",lwd=0.8)
lines(density(log2(correctPoints_H1$dist2Ext)),col="green",lwd=0.8)
dev.off()

#IMR90

png(paste(resultsDIR,"figure33_IMR90AllDist2Ext.png",sep=""),units="cm",height=8,width=8,res=300)
def.par <- par(no.readonly = TRUE)

par(lwd=1.5)
par(cex.axis=0.8)
plot(density(log2(correctPoints_IMR90$dist2Ext)),type="n",main="",xlab="",ylab="")
lines(density(log2(incorrectPoints_IMR90_minus$dist2Ext)),col="blue",lwd=0.8)
lines(density(log2(incorrectPoints_IMR90_plus$dist2Ext)),col="red",lwd=0.8)
lines(density(log2(correctPoints_IMR90$dist2Ext)),col="green",lwd=0.8)
dev.off()

# 2.2.4.2- CpG island length

#H1

png(paste(resultsDIR,"figure33_H1AllCpgiLength.png",sep=""),units="cm",height=8,width=8,res=300)
def.par <- par(no.readonly = TRUE)

par(lwd=1.5)
par(cex.axis=0.8)
plot(density(log2(correctPoints_H1$cgiLength)),type="n",main="",xlab="",ylab="")
lines(density(log2(incorrectPoints_H1_minus$cgiLength)),col="blue",lwd=0.8)
lines(density(log2(incorrectPoints_H1_plus$cgiLength)),col="red",lwd=0.8)
lines(density(log2(correctPoints_H1$cgiLength)),col="green",lwd=0.8)
dev.off()

#IMR90

png(paste(resultsDIR,"figure33_IMR90AllCpgiLength.png",sep=""),units="cm",height=8,width=8,res=300)
def.par <- par(no.readonly = TRUE)

par(lwd=1.5)
par(cex.axis=0.8)
plot(density(log2(correctPoints_IMR90$cgiLength)),type="n",main="",xlab="",ylab="")
lines(density(log2(incorrectPoints_IMR90_minus$cgiLength)),col="blue",lwd=0.8)
lines(density(log2(incorrectPoints_IMR90_plus$cgiLength)),col="red",lwd=0.8)
lines(density(log2(correctPoints_IMR90$cgiLength)),col="green",lwd=0.8)
dev.off()

# 2.2.4.3- CpG island nCpG

#H1
png(paste(resultsDIR,"figure33_H1AllCpginCG.png",sep=""),units="cm",height=8,width=8,res=300)
def.par <- par(no.readonly = TRUE)

par(lwd=1.5)
par(cex.axis=0.8)
plot(density(log2(correctPoints_H1$E_nCG)),type="n",main="",xlab="",ylab="")
lines(density(log2(incorrectPoints_H1_minus$E_nCG)),col="blue",lwd=0.8)
lines(density(log2(incorrectPoints_H1_plus$E_nCG)),col="red",lwd=0.8)
lines(density(log2(correctPoints_H1$E_nCG)),col="green",lwd=0.8)
dev.off()

#IMR90

png(paste(resultsDIR,"figure33_IMR90AllCpginCG.png",sep=""),units="cm",height=8,width=8,res=300)
def.par <- par(no.readonly = TRUE)

par(lwd=1.5)
par(cex.axis=0.8)
plot(density(log2(correctPoints_IMR90$E_nCG)),type="n",main="",xlab="",ylab="")
lines(density(log2(incorrectPoints_IMR90_minus$E_nCG)),col="blue",lwd=0.8)
lines(density(log2(incorrectPoints_IMR90_plus$E_nCG)),col="red",lwd=0.8)
lines(density(log2(correctPoints_IMR90$E_nCG)),col="green",lwd=0.8)
dev.off()

# 2.2.4.4- CpG island O/E

#H1
png(paste(resultsDIR,"figure33_H1AllCpgiOE.png",sep=""),units="cm",height=8,width=8,res=300)
def.par <- par(no.readonly = TRUE)

par(lwd=1.5)
par(cex.axis=0.8)
plot(density(log2(incorrectPoints_H1_minus$E_OE)),type="n",main="",xlab="",ylab="")
lines(density(log2(incorrectPoints_H1_minus$E_OE)),col="blue",lwd=0.8)
lines(density(log2(incorrectPoints_H1_plus$E_OE)),col="red",lwd=0.8)
lines(density(log2(correctPoints_H1$E_OE)),col="green",lwd=0.8)
dev.off()

#IMR90
png(paste(resultsDIR,"figure33_IMR90AllCpgiOE.png",sep=""),units="cm",height=8,width=8,res=300)
def.par <- par(no.readonly = TRUE)

par(lwd=1.5)
par(cex.axis=0.8)
plot(density(log2(incorrectPoints_IMR90_minus$E_OE)),type="n",main="",xlab="",ylab="")
lines(density(log2(incorrectPoints_IMR90_minus$E_OE)),col="blue",lwd=0.8)
lines(density(log2(incorrectPoints_IMR90_plus$E_OE)),col="red",lwd=0.8)
lines(density(log2(correctPoints_IMR90$E_OE)),col="green",lwd=0.8)
dev.off()















