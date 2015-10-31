'''
@author: Victor Barrera Burgos
Created on 09 Feb 2014
Description: This script permits the obtention of the
methylation profile of a CpGRegion indicating the 
methylation status of each CpG dinucleotide.

Addition on 02 March 2014
Description: permits the obtention of the
methylation profile of the whole genome using methylationMap.
'''

# Imports

import sys
import pysam
import re

# Defining types
# Structures from the first level of abstraction
class CpGRegion:
    def __init__(self,id,chrom,start,end,sequence):
        self.id=id
        self.chrom=chrom
        self.start=start
        self.end=end
        self.sequence=sequence
        self.methCoef=-1
        self.nCG=0
        self.cpgList=[]

# Structures from the second level of abstraction

class CpGdinucleotide:
    def __init__(self,chrom,firstPos,secondPos):
        self.chrom=chrom
        self.firstPos=firstPos
        self.secondPos=secondPos
        self.firstPosMethCoef=-1.0
        self.secondPosMethCoef=-1.0
        self.meanMethCoef=-1.0
        

# Defining functions
# Functions from the first level of abstraction

def methylationMap(cgR,filter):

# Pre: The function receives an object from the class CpGRegion and a filter value

    upSequence=(cgR.sequence).upper()

# Look for CG positions
    starts = [match.start() for match in re.finditer('CG',upSequence)]
    for i in starts:
        cpgDin=CpGdinucleotide(cgR.chrom,int(i)+cgR.start-1,int(i)+cgR.start)
# Call the methCG function
        methCG(cpgDin,filter)
        cgR.nCG=cgR.nCG+1
        (cgR.cpgList).append(cpgDin)
    cgRPositions=""
    for j in cgR.cpgList:
        if (j.meanMethCoef>=0):
            if (j.meanMethCoef<=0.2):
                cgRPositions=cgRPositions+str(j.firstPos)+";"+"0"+"\n"
            elif (j.meanMethCoef>=0.8):
                cgRPositions=cgRPositions+str(j.firstPos)+";"+"1"+"\n"
            else:
                cgRPositions=cgRPositions+str(j.firstPos)+";"+"X"+"\n"
        else:
            cgRPositions=cgRPositions+str(j.firstPos)+";"+"N"+"\n"
    print "%s" % (cgRPositions)
# Post: 

def methylationProfile(cgR,filter):

# Pre: The function receives an object from the class CpGRegion and a filter value

    upSequence=(cgR.sequence).upper()

# Look for CG positions
    starts = [match.start() for match in re.finditer('CG',upSequence)]
    for i in starts:
        cpgDin=CpGdinucleotide(cgR.chrom,int(i)+cgR.start-1,int(i)+cgR.start)
# Call the methCG function
        methCG(cpgDin,filter)
        cgR.nCG=cgR.nCG+1
        (cgR.cpgList).append(cpgDin)
# Generate the profile using the code (0,1,X,N)
# 0 For unmethylated, 1 for methylated
# X for intermediate methylation, N for not informative
    cgRProf=""
    infCpG=0
    cgRAcumMethCoef=0
    for j in cgR.cpgList:
        if (j.meanMethCoef>=0):
            infCpG=infCpG+1
            cgRAcumMethCoef=cgRAcumMethCoef+j.meanMethCoef
            if (j.meanMethCoef<=0.2):
                cgRProf=cgRProf+"0"
            elif (j.meanMethCoef>=0.8):
                cgRProf=cgRProf+"1"
            else:
                cgRProf=cgRProf+"X"
        else:
            cgRProf=cgRProf+"N"
    if (infCpG>0):
        cgR.methCoef=cgRAcumMethCoef/infCpG
    
    print "%s;%s;%i;%i;%i;%i;%f;%s" % (cgR.id,cgR.chrom,cgR.start,cgR.end,cgR.nCG,infCpG,cgR.methCoef,cgRProf)
# Post: The id, chrom, start, end, total number of CG, number of informative CpG
# and a ternary profile for each of its CpG corresponding to the CpGRegion object
# have been printed

# Functions from the second level of abstraction

def methCG(cpgDin,filter):
# Pre: The function receives an object from the class CpGdinucleotide and a filter value
    seq=""
    for pileupcolumn in samfile.pileup(cpgDin.chrom,cpgDin.firstPos,cpgDin.firstPos+1):
        if not (pileupcolumn.pos==cpgDin.firstPos and pileupcolumn.n>=filter):
            continue
        for pileupread in pileupcolumn.pileups:
            seq+=pileupread.alignment.seq[pileupread.qpos]
        seq=seq.upper()
        numA=int(seq.count("A"))
        numT=int(seq.count("T"))
        numC=int(seq.count("C"))
        numG=int(seq.count("G"))
        reads=numA+numT+numC+numG
        if ((numT+numC)>=filter):
            cpgDin.firstPosMethCoef=(float(numC)/float(numC+numT))

    seq=""
    for pileupcolumn in samfile.pileup(cpgDin.chrom,cpgDin.secondPos,cpgDin.secondPos+1):
        if not (pileupcolumn.pos==cpgDin.secondPos and pileupcolumn.n>=filter):
            continue
        for pileupread in pileupcolumn.pileups:
            seq+=pileupread.alignment.seq[pileupread.qpos]
        seq=seq.upper()
        numA=int(seq.count("A"))
        numT=int(seq.count("T"))
        numC=int(seq.count("C"))
        numG=int(seq.count("G"))
        reads=numA+numT+numC+numG
        if ((numT+numC)>=filter):
            cpgDin.secondPosMethCoef=(float(numC)/float(numC+numT))
  
    if (((cpgDin.firstPosMethCoef)!=-1) & ((cpgDin.secondPosMethCoef)==-1)):
        cpgDin.meanMethCoef=cpgDin.firstPosMethCoef

    elif (((cpgDin.firstPosMethCoef)==-1) & ((cpgDin.secondPosMethCoef)!=-1)):
        cpgDin.meanMethCoef=cpgDin.secondPosMethCoef

    else:
        cpgDin.meanMethCoef=float(cpgDin.firstPosMethCoef+cpgDin.secondPosMethCoef)/2.0
# Post: The object is returned with its methylation Coefficients recalculated according
# to the data present in the alignment file  and using a minimum read filter.
 

####################
### Main ###########
####################

# Obtain the files
cpgr_sec_path=sys.argv[1]
sam_path=sys.argv[2]
filter=int(sys.argv[3])

# Load the files
cpgr_sec_file=open(cpgr_sec_path,'r')
samfile = pysam.Samfile(sam_path, "rb" )

for cpgr in cpgr_sec_file:
    cgRTuple=cpgr.split()
    cgR=CpGRegion(cgRTuple[0],str(cgRTuple[1]),int(cgRTuple[2]),int(cgRTuple[3]),str(cgRTuple[4]))
# We can use methylationMap or methylationProfile
    methylationMap(cgR,filter)
