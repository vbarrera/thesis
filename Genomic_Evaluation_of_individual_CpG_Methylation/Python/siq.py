'''
Created on 15 Oct, 2010

@author: vbarrera
'''

### siq (Sam Index Queries) is an auxiliary module for 
### queries made to a Sam Index.
### This module will be the template for object creation and
### related methods.
############################################################
### v1.0 ():
import pysam
import re
import math

class Region:
    '''
    Object containing general information about the region
    with methods to obtain the mapped sequences from a sam
    index
    '''

# __init__ method: This method is called during the creation
# of the object from the class Region. It just assigns the  
# desired coordinates to the associated samfile to work
# with the provided coordinates 

    def __init__(self,array_coord,samfile):
        self.array_coordinates=array_coord
        self.samfile=samfile


# table method: with a number of reads as threshold filter, this
# method scans the positions listed in the atribute "coordinates" and returns
# the number of each nucleotide in that position acording to the samfile
# supplied during the class creation.                
# This method is mainly used to check the calculation of the region mean.


    def table(self,nreads_filter):
        self.nreads_filter=nreads_filter
        self.data=[]
        for self.specific_coord in self.array_coordinates:
            self.chr=(self.specific_coord.split()[0])
            self.start_position=int(self.specific_coord.split()[1])
            self.seq=""
            for pileupcolumn in self.samfile.pileup(self.chr,self.start_position,self.start_position+1):
                
                if not (pileupcolumn.pos==self.start_position and pileupcolumn.n>=self.nreads_filter):
                    continue

                for pileupread in pileupcolumn.pileups:
                    self.seq+=pileupread.alignment.seq[pileupread.qpos]
            
            self.numA=int(self.seq.count("A"))
            self.numT=int(self.seq.count("T"))
            self.numC=int(self.seq.count("C"))
            self.numG=int(self.seq.count("G"))
            
            
            self.sub_data=[str(self.chr),self.start_position,\
                           self.numA,self.numT,self.numC,self.numG]
            self.data.append(self.sub_data)
        return self.data
    





# We generate a new class, Meth_region wich inherits
# methods and atributes from region (including __init__)
# This new class is used for Methylation parameters obtention

class Meth_region(Region):
    
# methcoef method:  This method is used to calculate the 
# methylation coeficient, defined as number of Cytosines
# divided by the sum of the number of cytosines and the
# number of thymines. This methylation coeficient
# is calculated for each position indicated in the 
# coordinates atribute. If the region contains multiple 
# coordinates, the value for the methcoef is the 
# mean value for the informative positions.
# Useful for point position (not region). 

    def methcoef(self,nreads_filter):
        self.nreads_filter=nreads_filter
        self.region_methcoef=0
        self.pos_inf=0
        self.total_reads=0       
        for self.specific_coord in self.array_coordinates:
            self.chr=(self.specific_coord.split()[0])
            self.start_position=int(self.specific_coord.split()[1])
            self.seq=""
            for pileupcolumn in self.samfile.pileup(self.chr,self.start_position,self.start_position+1):
                
                if not (pileupcolumn.pos==self.start_position and pileupcolumn.n>=self.nreads_filter):
                    continue
                
                for pileupread in pileupcolumn.pileups:
                    
                    self.seq+=pileupread.alignment.seq[pileupread.qpos]
            self.seq=self.seq.upper()
            self.numA=int(self.seq.count("A"))
            self.numT=int(self.seq.count("T"))
            self.numC=int(self.seq.count("C"))
            self.numG=int(self.seq.count("G"))
            self.reads=self.numA+self.numT+self.numC+self.numG
            if (self.numT+self.numC>=self.nreads_filter):
                self.pos_inf+=1
                self.methcoef=float(self.numC)/float(self.numC+self.numT)
            else:
                self.methcoef=-1
                #IMPORTANT! FOR REGIONS WORK WITH METHCOEF_SD METHOD
            self.total_reads+=self.reads    
            self.region_methcoef+=self.methcoef
        if (self.pos_inf>0):
            self.region_methcoef=self.region_methcoef/self.pos_inf    
        return (self.pos_inf,self.region_methcoef,self.total_reads)    


# method methcoef_sd: this method is a variant from 
# methcoef. It computes the mean value of the methylation coeficient
# which corresponds exactly with the methcoef method return value. The
# new feature is the calculus of the standard deviation of the methcoef
# between the different coordinates that define the region.           
    def methcoef_sd(self,nreads_filter):
        self.nreads_filter=nreads_filter
        self.region_methcoef=0
        self.pos_inf=0
        self.total_reads=0       
        self.regionMethcoeflist=[]
        for self.specific_coord in self.array_coordinates:
            self.chr=(self.specific_coord.split()[0])
            self.start_position=int(self.specific_coord.split()[1])
            self.seq=""
            
            for pileupcolumn in self.samfile.pileup(self.chr,self.start_position,self.start_position+1):
                
                if not (pileupcolumn.pos==self.start_position and pileupcolumn.n>=self.nreads_filter):
                    continue
                
                for pileupread in pileupcolumn.pileups:
                    
                    self.seq+=pileupread.alignment.seq[pileupread.qpos]
            self.seq=self.seq.upper()
            self.numA=int(self.seq.count("A"))
            self.numT=int(self.seq.count("T"))
            self.numC=int(self.seq.count("C"))
            self.numG=int(self.seq.count("G"))
            self.reads=self.numA+self.numT+self.numC+self.numG
            
            if (self.numT+self.numC>=self.nreads_filter):
                self.pos_inf+=1
                self.methcoef=float(self.numC)/float(self.numC+self.numT)
                self.regionMethcoeflist.append(self.methcoef)
            self.total_reads+=self.reads    
            
        if (self.pos_inf>1):
            self.sum=0.0
            
            for i in self.regionMethcoeflist:
                self.sum+=i
            
            self.mean=self.sum/(float(self.pos_inf))
            
            self.dev=0.0
            for j in self.regionMethcoeflist:
                self.dev+=(j-self.mean)*(j-self.mean)
                
            self.var=self.dev/(self.pos_inf-1)
            self.sd=math.sqrt(self.var)
        else:
            self.mean=-1.0
            self.sd=25.0
         
        return (self.pos_inf,self.mean,self.sd, self.total_reads)
 
