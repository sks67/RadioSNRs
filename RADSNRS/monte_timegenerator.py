#######################################################################
# PURPOSE: Generates Times based on star formation rate probability.
# AUTHOR: Sumit 
#######################################################################


from numpy import *
import random as rn
import sys


#prints the program name during runtime

progname = sys.argv[0]

tmax = 6.0e6
def timegen(snrate):
   time = []
   j = 0
   count = 0
   for i in range(tmax):
       prob = rn.random()
       if prob > snrate:
           continue
       else:
           count=count+1
           time.append(i)
           j=j+1
   
   return np.array(time)    
  

           
