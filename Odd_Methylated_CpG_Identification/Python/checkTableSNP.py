'''
@author: Victor Barrera Burgos
Created on 9th August 2015
Description: This script uses the "mpileup" output from samtools
to check if a position corresponds to a SNP or not. 
'''
# Imports

import sys
import re

# Subroutines
# Pre: The function receives a string
def countNucleotides(nuclString):
        seq=nuclString.upper()
        numA=int(seq.count("A"))
        numT=int(seq.count("T"))
        numC=int(seq.count("C"))
        numG=int(seq.count("G"))
        return [numA,numT,numC,numG]
# Post: A tuple with the count number of each nucleotide
# is returned

# Pre: The function receives a string corresponding to 
# a column (a position in the genome).
def processReadsColumn(readCol):
    readCol=readCol.strip("\r\n")
    # Eliminate ^~$ symbols
    readCol=readCol.translate(None,'^~$')
    # Separate Forward (Major case) and reverse (minor case)
    forwardVals="".join(re.findall('[ATCG]', readCol))
    reverseVals="".join(re.findall('[atcg]', readCol))
    return [countNucleotides(forwardVals),countNucleotides(reverseVals)]
# Post: The function returns a count of each nucleotide for 
# the forward and reverse strand for a position

# Main

snpPosReadsFile=sys.argv[1]
snpPosReadsFH=open(snpPosReadsFile,"r")

while snpPosReadsFH:
    firstLine=snpPosReadsFH.readline().strip("\r\n")
    firstReadCol=firstLine.split()[4]
    secondLine=snpPosReadsFH.readline().strip("\r\n")
    secondReadCol=secondLine.split()[4]
    pos=str(firstLine.split()[0])+"_"+str(firstLine.split()[1])
    # Get the values for each position and direction and process them
    firstColFwd=(processReadsColumn(firstReadCol))[0]
    firstColRvs=(processReadsColumn(firstReadCol))[1]
    secondColFwd=(processReadsColumn(secondReadCol))[0]
    secondColRvs=(processReadsColumn(secondReadCol))[1]
    
    # Is it a SNP?
    SNP="NoSNP"
    if ((firstColRvs[0] > firstColRvs[3]) or (secondColFwd[0] > secondColFwd[3])):
        SNP="SNP"
    # Is it trustable
    trustable="T"
    if ((firstColRvs[0]==0 and firstColRvs[3]==0) and (secondColFwd[0]==0 and secondColFwd[3]==0)):
        trustable="NT" 
    # Any Warning
    warning="NW"
    if ((abs(firstColRvs[0] - firstColRvs[3])<10) or (abs(secondColFwd[0] - secondColFwd[3]))<10):
        warning="W"

# Print the results       
    firstColString=str(firstColFwd[0])+"\t"+str(firstColFwd[1])+"\t"+str(firstColFwd[2])+"\t"+str(firstColFwd[3])+"\t\t"+ str(firstColRvs[0])+"\t"+str(firstColRvs[1])+"\t"+str(firstColRvs[2])+"\t"+str(firstColRvs[3])
    secondColString=str(secondColFwd[0])+"\t"+str(secondColFwd[1])+"\t"+str(secondColFwd[2])+"\t"+str(secondColFwd[3])+"\t\t"+ str(secondColRvs[0])+"\t"+str(secondColRvs[1])+"\t"+str(secondColRvs[2])+"\t"+str(secondColRvs[3])

    print str(pos)+"\t"+str(firstColString)+"\t\t"+str(secondColString)+"\t\t"+str(SNP)+"\t"+str(trustable)+"\t"+str(warning)
    


