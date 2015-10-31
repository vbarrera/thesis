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
selected_reDinucleotide_Element_ADS_IPSC<-subset(reDinucleotide_Element_ADS_IPSC,
                                                 reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC!=-1 & 
                                                   reDinucleotide_Element_ADS_IPSC$E_methCoef!=-1 & ((reDinucleotide_Element_ADS_IPSC$E_posInf/2)/reDinucleotide_Element_ADS_IPSC$E_nCG)>=0.25)




incorrectPoints_H1_plus<-subset(selected_reDinucleotide_Element_H1,selected_reDinucleotide_Element_H1$E_methCoef-selected_reDinucleotide_Element_H1$meanMethCoef_H1>0.25)
incorrectPoints_H1_minus<-subset(selected_reDinucleotide_Element_H1,selected_reDinucleotide_Element_H1$E_methCoef-selected_reDinucleotide_Element_H1$meanMethCoef_H1<(-0.25))
discordantH1<-c(incorrectPoints_H1_plus$id_In_RE,incorrectPoints_H1_minus$id_In_RE)
incorrectPoints_IMR90_plus<-subset(selected_reDinucleotide_Element_IMR90,selected_reDinucleotide_Element_IMR90$E_methCoef-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90>0.25)
incorrectPoints_IMR90_minus<-subset(selected_reDinucleotide_Element_IMR90,selected_reDinucleotide_Element_IMR90$E_methCoef-selected_reDinucleotide_Element_IMR90$meanMethCoef_IMR90<(-0.25))
discordantIMR90<-c(incorrectPoints_IMR90_plus$id_In_RE,incorrectPoints_IMR90_minus$id_In_RE)
incorrectPoints_ADS_plus<-subset(selected_reDinucleotide_Element_ADS,selected_reDinucleotide_Element_ADS$E_methCoef-selected_reDinucleotide_Element_ADS$meanMethCoef_ADS>0.25)
incorrectPoints_ADS_minus<-subset(selected_reDinucleotide_Element_ADS,selected_reDinucleotide_Element_ADS$E_methCoef-selected_reDinucleotide_Element_ADS$meanMethCoef_ADS<(-0.25))
discordantADS<-c(incorrectPoints_ADS_plus$id_In_RE,incorrectPoints_ADS_minus$id_In_RE)
incorrectPoints_ADS_Adipose_plus<-subset(selected_reDinucleotide_Element_ADS_Adipose,selected_reDinucleotide_Element_ADS_Adipose$E_methCoef-selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose>0.25)
incorrectPoints_ADS_Adipose_minus<-subset(selected_reDinucleotide_Element_ADS_Adipose,selected_reDinucleotide_Element_ADS_Adipose$E_methCoef-selected_reDinucleotide_Element_ADS_Adipose$meanMethCoef_ADS_Adipose<(-0.25))
discordantADS_Adipose<-c(incorrectPoints_ADS_Adipose_plus$id_In_RE,incorrectPoints_ADS_Adipose_minus$id_In_RE)
incorrectPoints_ADS_IPSC_plus<-subset(selected_reDinucleotide_Element_ADS_IPSC,selected_reDinucleotide_Element_ADS_IPSC$E_methCoef-selected_reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC>0.25)
incorrectPoints_ADS_IPSC_minus<-subset(selected_reDinucleotide_Element_ADS_IPSC,selected_reDinucleotide_Element_ADS_IPSC$E_methCoef-selected_reDinucleotide_Element_ADS_IPSC$meanMethCoef_ADS_IPSC<(-0.25))
discordantADS_IPSC<-c(incorrectPoints_ADS_IPSC_plus$id_In_RE,incorrectPoints_ADS_IPSC_minus$id_In_RE)

length(discordantH1)
length(discordantIMR90)
length(discordantADS)
length(discordantADS_Adipose)
length(discordantADS_IPSC)

commonADS_Adipose<-intersect(discordantADS,discordantADS_Adipose)
length(commonADS_Adipose)
commonADS_IPSC<-intersect(discordantADS,discordantADS_IPSC)
length(commonADS_IPSC)
commonAdipose_IPSC<-intersect(discordantADS_Adipose,discordantADS_IPSC)
length(commonAdipose_IPSC)
commonADS_Adipose_IPSC<-intersect(commonADS_Adipose,discordantADS_IPSC)
length(commonADS_Adipose_IPSC)

commonIMR90_ADS_Adipose_IPSC<-intersect(commonADS_Adipose_IPSC,discordantIMR90)
length(commonIMR90_ADS_Adipose_IPSC)
commonH1_IMR90_ADS_Adipose_IPSC<-intersect(commonIMR90_ADS_Adipose_IPSC,discordantH1)
length(commonH1_IMR90_ADS_Adipose_IPSC)

