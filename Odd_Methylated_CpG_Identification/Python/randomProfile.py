'''
@author: Victor Barrera Burgos
Created on 19th May 2014
Description: This script processes takes a profile file
and making use of random.sample and join functions generates 
a reordered new one.

Change on 20th May 2014: Now we count the 0 and 1 from the profile 
and generate a list that we shuffle (using random.shuffle). 
Afterwards and using the original profile as a template 
we generate a random profile. It works setting a pointer 
that travels across the original profile. When the pointer 
is on a 0 or 1, it appends to the random profile being generated 
a value from the shuffled list. If it is on a X o N it appends 
that value. With this, we only change the 0 and 1 regions while 
the X and N remain the same as in the original profile.
'''

import random
import sys

profileFile=sys.argv[1]
profileFH=open(profileFile,"r")

# Obtain number of 1s and 0s for a given profile
profileString=profileFH.readline().strip("\r\n")
numZeros=int(profileString.count("0"))
numOnes=int(profileString.count("1"))

#generate False list
valuesList=([0]*numZeros)+([1]*numOnes)
random.shuffle(valuesList)
# Change on 20th May, see description
#randomProfile=''.join(random.sample(profileString,len(profileString)))


i=0
j=0
newProfile=""
endText=(i==len(profileString) or (j==len(valuesList)))

while not endText:
    if profileString[i] == '1' or profileString[i] =='0':
        newProfile=newProfile+str(valuesList[j])
        i+=1
        j+=1
        endText=(i==len(profileString) or (j==len(valuesList)))
    else:
        newProfile=newProfile+str(profileString[i])
        i+=1
    endText=(i==len(profileString))
    
    
        
profileFH.close()

# Print the random profile
print newProfile