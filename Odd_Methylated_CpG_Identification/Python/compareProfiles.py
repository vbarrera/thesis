'''
@author: Victor Barrera Burgos
Created on 28th May 2014
Description: This script allows the comparison of two 
methylation profiles (A and B) looking for specific sequences changes.
In the patternFile all the patterns to look for are specified (one per line).
desirePatFile specifies which pattern changes must be reported.
When a pattern from patternsFile is found in profileA, the sequence from profile B
corresponding to the same genomic coordinates is obtain. If this sequence is among 
the desired patterns in desirePatFile, this change is reported 
'''

#imports

import sys
import re
# SubRoutines

# Pre: the function receives a string to scan and a 
# a pattern to search (sub). 
def allindices(string, sub, listindex=[], offset=0):
    i = string.find(sub, offset)
    while i >= 0:
        listindex.append(i)
        i = string.find(sub, i + 1)
    return listindex
# Post: a list of all the start positions 
# of the pattern (sub) 


# variables.

profileAFile=sys.argv[1]
profileBFile=sys.argv[2]
patternsFile=sys.argv[3]
desirePatFile=sys.argv[4]
compareFile=sys.argv[5]
desPatPosFile=sys.argv[6]


profAFH=open(profileAFile,"r")
profBFH=open(profileBFile,"r")
patternsFH=open(patternsFile,"r")
desiredPatFH=open(desirePatFile,"r")
compareFH=open(compareFile,"w")
desPatPosFH=open(desPatPosFile,"w")


# Read the profiles
profileA=str(profAFH.readline().strip("\r\n"))
profileB=str(profBFH.readline().strip("\r\n"))
# Create the list of desired changes
desiredPatList=[]
for desPat in desiredPatFH:
    desPat=desPat.strip("\r\n")
    desiredPatList.append(desPat)
hitsProfB={}
desiredPatPos=[]

# Scan the profile A for each pattern in patternsFile
for pattern in patternsFH:
    pattern=pattern.strip("\r\n")
    lengthPattern=len(pattern)
    starts=allindices(profileA,str(pattern),[],0)
    for start in starts:
        # Obtain the corresponding sequence in profile B and add it
        # to the hitsProfB dictionary
        profBpatt=profileB[start:start+lengthPattern]
        if profBpatt in hitsProfB.keys():
            hitsProfB[profBpatt][0]=hitsProfB[profBpatt][0]+1
        else:
            hitsProfB[profBpatt]=[1,len(profBpatt)]
        # If the sequence in profile B in the desired changes, report the 
        # the position 
        if profBpatt in desiredPatList:
            desiredPatPos.append(str(start))
            
# Print results: 
for sKey in sorted(hitsProfB.keys()):
    compareFH.write((sKey)+"\t"+str(hitsProfB[sKey][0])+"\t"+str(hitsProfB[sKey][1])+"\n")

desPatPosFH.write(";".join(desiredPatPos))    

profAFH.close()
profBFH.close()
patternsFH.close()
compareFH.close()
desPatPosFH.close()



