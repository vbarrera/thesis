'''
@author: Victor Barrera
Description: This scripts makes use of the siq module to obtain the methylation value for the two positions of
a restriction site.
'''

import sys
import pysam
import re
from siq import *

# Obtain the Restriction Enzime positions file,
# the sam-file and the minimum read filter

resEnzPos_path=sys.argv[1]
sam_path=sys.argv[2]
filter=int(sys.argv[3])

resEnzPos_file=open(resEnzPos_path,'r')
samfile = pysam.Samfile(sam_path, "rb" )

for resEnzPos in resEnzPos_file:
    id=resEnzPos.split()[0]
    chr=str(resEnzPos.split()[1])
    position=int(resEnzPos.split()[2])
    CG_coord_wat=[]
    CG_coord_crick=[]
    c_pos=str(chr)+"\t"+str(position-1)
    g_pos=str(chr)+"\t"+str(position)
    CG_coord_wat.append(c_pos)
    CG_coord_crick.append(g_pos)
# We generate a Region object with the coordinates for each of the two nucleotides of the site
    cg_restEnzPos_wat=Meth_region(CG_coord_wat,samfile)
    cg_restEnzPos_crick=Meth_region(CG_coord_crick,samfile)
# We obtain the methylation value for each position.
    print str("C")+"\t"+str(id)+"\t",
    print "%i\t%.2f\t%i" %(cg_restEnzPos_wat.methcoef(filter))
    print str("G")+"\t"+str(id)+"\t",
    print "%i\t%.2f\t%i" %(cg_restEnzPos_crick.methcoef(filter))
    
