'''
@author: Victor Barrera Burgos
Created on 10 March 2014
Description: This script scans a genomic methylation
profile (coded with 1,0,N,X) for uninterrupted sequences
composed only by 1s and/or 0s. It has two output files:
- a frequency file where each uninterrupted sequence is 
indicated along with the number of times that it appears and 
the size of the sequence.
- a size distribution file that summarizes the number of times
a certain size of sequences appears.
'''

import sys

# Main

# Obtain the files
profileFile=sys.argv[1]
profileFH=open(profileFile,"r")
sizeDistFile=sys.argv[2]
sizeDistFH=open(sizeDistFile,"w")
freqProfFile=sys.argv[3]
freqProfFH=open(freqProfFile,"w")


chromosomeProfile=profileFH.readline().strip("\r\n")

profilesDict={}
profilesSizeDistDict={}

i=0
j=0
endText=(i==len(chromosomeProfile) or (j==len(chromosomeProfile)))

# Scan the profile
while not endText:
    if chromosomeProfile[i] == '1' or chromosomeProfile[i] =='0':
        j=i+1
        # Finish scanning when it arrives to the end of the profile
        endText=(i==len(chromosomeProfile) or (j==len(chromosomeProfile)))
        # Store the values until it arrives to a value different of 1 or 0
        while ((not endText) and (chromosomeProfile[j] == '1' or chromosomeProfile[j] =='0')):
            j+=1
            endText=(i==len(chromosomeProfile) or (j==len(chromosomeProfile))) 
        prof=chromosomeProfile[i:j]
        # Add it to the dictionaries profileDict and profileSizeDistDict. If it exists, sum 1. 
        # If not, create a new entry in the dictionaries
        if prof in profilesDict.keys():
            profilesDict[prof][0]+=1
        else:
            profilesDict[prof]=[1,len(prof)]
        if str(len(prof)) in profilesSizeDistDict.keys():
            profilesSizeDistDict[str(len(prof))]+=1
        else:
            profilesSizeDistDict[str(len(prof))]=1
        i=j
    else:
        i+=1
    endText=(i==len(chromosomeProfile) or (j==len(chromosomeProfile)))   


# Report results
# Frequency files: sequence \t number of times \t size
for sKey in sorted(profilesDict.keys()):
    freqProfFH.write((sKey)+"\t"+str(profilesDict[sKey][0])+"\t"+str(profilesDict[sKey][1])+"\n")
# Size distribution file: size \t number of times
for sKey in sorted(profilesSizeDistDict.keys(),key=lambda x: int(x)):
    sizeDistFH.write(sKey+"\t"+str(profilesSizeDistDict[sKey])+"\n")

profileFH.close()
freqProfFH.close()
sizeDistFH.close()
