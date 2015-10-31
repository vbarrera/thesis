'''
@author: Victor Barrera
Description: This scripts makes use of the siq module to obtain the methylation value for each CpG island region.

'''
import sys
import pysam
import re
from siq import *

# Obtain the CpG island sequence file,
# the sam-file and the minimum read filter
cpgi_sec_path=sys.argv[1]
sam_path=sys.argv[2]
filter=int(sys.argv[3])


cpgi_sec_file=open(cpgi_sec_path,'r')
samfile = pysam.Samfile(sam_path, "rb" )

for cpgi in cpgi_sec_file:

# For each CpG island structural data is obtained and the relative
# positions for the CG dinucleotide are obtained.

    id=cpgi.split()[0]
    chr=str(cpgi.split()[1])
    startPosition=int(cpgi.split()[2])
    endPosition=int(cpgi.split()[3])
    cpgisec=(cpgi.split()[4])
    cpgisec=cpgisec.upper();
    starts = [match.start() for match in re.finditer('CG',cpgisec)]
    CG_coord=[]
    nCG=0
    for i in starts:
# The absolute positions for the CG dinucleotide are obtain.
        c_pos=str(chr)+"\t"+str(int(i)+startPosition-1)
        g_pos=str(chr)+"\t"+str(int(i)+startPosition)
        CG_coord.append(c_pos)
        CG_coord.append(g_pos)
        nCG+=1
# An object of the class Meth_region is generated.
    cgi=Meth_region(CG_coord,samfile)
    print str(id)+"\t"+str(nCG)+"\t"+str(chr)+"\t"+str(startPosition)+"\t"+str(endPosition)+"\t",
# A call to obtain the mean and the standard deviation is done.
    print "%i\t%.2f\t%.2f\t%i" %(cgi.methcoef_sd(filter))
