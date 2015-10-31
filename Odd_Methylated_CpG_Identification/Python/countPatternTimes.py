'''
@author: Victor Barrera Burgos
Created on 21th August 2015
Description: The frequency files report each uninterrupted string.
However, that string can contain more than one ocurrence of the 
desired sequence.
For example, the uninterrupted string:
11111011111011111 contains two times the sequence 11111011111 (with an overlap).
This string allows to obtain, for each string, the real number of appearances of 
the desired sequence.
'''

import sys
import re

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

freqFile=sys.argv[1]
freqFH=open(freqFile,"r")

for line in freqFH:
    line = line.strip("\t\n")
    array=line.split("\t")
    
    timesAppear=len(allindices(array[0],"11111011111",[],0))
    #timesAppear=len(allindices(array[0],"00000100000",[],0))
    if timesAppear>0:
        print int(timesAppear)*int(array[1])

freqFH.close()