##
# Description: Code to generate the graphs used in the thesis:
# Metilación del ADN en posiciones individuales del genoma humano
# dentro de contextos epignenéticos regionales
# Section (and paper):
# Evaluation of single CpG sites as proxies of CpG island methylation
# states at the genome scale
# 
# It contains the analysis for CpGi with and without the filter of
# informativeness.

#########################################################################
# Used Folders:

# Function to create folders

createDIR <-function (DIR){
  if(file.exists(DIR)){
    
  }else{
    dir.create(file.path(DIR))
  }
}

baseDIR<-paste(Sys.getenv("HOME"),"/projects/Evaluation_of_single_CpG_sites_as_proxies_of_CpG_island_methylation_states_at_the_genome_scale/",sep="")
resultsDIR<-paste(baseDIR,"results/20120817/",sep="")
resultsDIRTMP<-resultsDIR

createDIR(resultsDIR)

numbered='FALSE'


# subDIR <- paste(resultsDIRTMP,"fig22/",sep="")
# createDIR(subDIR)
# resultsDIR<-subDIR
# source("figure22_H1_IMR90.R")
# source("figure22_rest.R")
# 
# subDIR <- paste(resultsDIRTMP,"fig23/",sep="")
# createDIR(subDIR)
# resultsDIR<-subDIR
# source("figure23_H1_IMR90.R")
# source("figure23_ADS_ADSC.R")
#
# subDIR <- paste(resultsDIRTMP,"fig24/",sep="")
# createDIR(subDIR)
# resultsDIR<-subDIR
# source("figure24.R")
# 
# subDIR <- paste(resultsDIRTMP,"fig25/",sep="")
# createDIR(subDIR)
# resultsDIR<-subDIR
# source("figure25.R")
# 
# subDIR <- paste(resultsDIRTMP,"fig26/",sep="")
# createDIR(subDIR)
# resultsDIR<-subDIR
# source("figure26.R")
#  
# subDIR <- paste(resultsDIRTMP,"fig27_fig28/",sep="")
# createDIR(subDIR)
# resultsDIR<-subDIR
# source("figure27_H1_IMR90.R")
# source("figure27_H1_IMR90_means.R")
# source("figure27_28_rest.R")
# source("figure27_ADS_ADSC_ADSI_means.R")
#  
# subDIR <- paste(resultsDIRTMP,"fig29/",sep="")
# createDIR(subDIR)
# resultsDIR<-subDIR
# source("figure29_Markov.R")
# source("figure29_Illing.R")
# # 
# subDIR <- paste(resultsDIRTMP,"fig30/",sep="")
# createDIR(subDIR)
# resultsDIR<-subDIR
# source("figure30_H1_IMR90.R")
# source("figure30_ADS_ADSC.R")
# 
# subDIR <- paste(resultsDIRTMP,"fig31/",sep="")
# createDIR(subDIR)
# resultsDIR<-subDIR
# source("figure31_H1_IMR90.R")
# source("figure31_ADS_ADSC_ADSCI.R")
# 
# subDIR <- paste(resultsDIRTMP,"fig32/",sep="")
# createDIR(subDIR)
# resultsDIR<-subDIR
# No graph generated, execute the code manually to get numbers
#
# 
# subDIR <- paste(resultsDIRTMP,"fig33/",sep="")
# createDIR(subDIR)
# resultsDIR<-subDIR
# source("figure33.R")
# 
