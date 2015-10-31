getwd()

setwd("input/")
x<-read.table("H1_HMM_States.bed_forRegioneR",header=T)
groups<-data.frame(stateCode=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
                       majorGroup=c(rep("Promoter",3),rep("Enhacer",4),"Insulator",
                                    rep("Transcription",3),"Repressed",
                                    rep("Inactive",3)),
                       minorGroup=c("Active_P","Weak_P","Poised_P",
                                    "Strong_E","Strong_E","Weak_E","Weak_E",
                                    "Insulator","Transcription","Transcription",
                                    "Weak_T","Polycomb_R","Heterochromatin",
                                    "Repetitive_CNV","Repetitive_CNV"))

mX<-merge(x,groups)

mXForBed<-mX[,c(2,3,4,5,1,6,7,8)]

write.table(mXForBed,file="H1_HMM_States_WithGroups.bed_forRegioneR",row.names=F,col.names=T,sep="\t",
            quote=F)