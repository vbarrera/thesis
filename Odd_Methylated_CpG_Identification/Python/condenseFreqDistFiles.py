'''
@author: Victor Barrera Burgos
Created on 19th May 2014
Description: This script processes the Freq and Dist 
files (two for each combination Chrom-CellLine) and provides 
condensed Freq and Dist files with the counting of all the patterns 
that appear among all the chromosomes (X and Y excluded) of a cellLine. 
(in a pattern file) is given.
'''

chrList=["chr1","chr2MOD"]
for val in range(3,23):
    chrList.append("chr"+str(val))

# Process the Frequency Files
for cellLine in ["H1","IMR90"]:
    freqDict={}
    outputFileName="condensedFreq_"+str(cellLine)+".txt"
    for chromosome in chrList:
        chromosomeFileName="methProf_CGi_hg18_"+str(cellLine)+"_CG_"+str(chromosome)+"_Freq.txt"
        chromoFH=open(chromosomeFileName,"r")
        for line in chromoFH:
            line=line.strip("\r\n")
            if str(line.split("\t")[0]) in freqDict.keys():
                freqDict[str(line.split("\t")[0])][0]=int(freqDict[str(line.split("\t")[0])][0])+int(line.split("\t")[1])
            else:
                freqDict[str(line.split("\t")[0])]=[line.split("\t")[1],line.split("\t")[2]]
       
        chromoFH.close()
    outputFH=open(outputFileName,"w")    
    # Print results (one file per cell line)
    for sKey in sorted(freqDict.keys()):
        outputFH.write((sKey)+"\t"+str(freqDict[sKey][0])+"\t"+str(freqDict[sKey][1])+"\n")
    outputFH.close()

# Process the size distribution Files
for cellLine in ["H1","IMR90"]:
    distDict={}
    outputFileName="condensedDist_"+str(cellLine)+".txt"
    for chromosome in chrList:
        chromosomeFileName="methProf_CGi_hg18_"+str(cellLine)+"_CG_"+str(chromosome)+"_Dist.txt"
        chromoFH=open(chromosomeFileName,"r")
        for line in chromoFH:
            line=line.strip("\r\n")
            if str(line.split("\t")[0]) in distDict.keys():
                distDict[str(line.split("\t")[0])]=int(distDict[str(line.split("\t")[0])])+int(line.split("\t")[1])
            else:
                distDict[str(line.split("\t")[0])]=line.split("\t")[1]
       
        chromoFH.close()
    # Print results (one file per cell line)
    outputFH=open(outputFileName,"w")
    for sKey in sorted(distDict.keys(),key=lambda x: int(x)):
        outputFH.write(sKey+"\t"+str(distDict[sKey])+"\n")
    outputFH.close()