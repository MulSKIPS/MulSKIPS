import sys
import os
import string
import numpy
import xml.etree.ElementTree
import shutil
import time
import math
import subprocess
filein=open("Test_NiSiC_SiC_t_max.csv","r")
fileout=open("GR_NiSiC_SiC_t_max.csv","w")
timeinter=100.*1.62799E-7
fileout.write('<Pulse xmlns="http://www.screen-lasse.com/LIAB/Pulse">'+'\n')
time=[]
powperc=[]
for line in filein:
    words=line.split()
    time.append(float(words[0]))
    powperc.append(float(words[1]))
linenum=len(time)
print linenum
timeinter=0.
i=1
while i <= linenum-1 :
      timeinter=timeinter+powperc[i]*(time[i]-time[i-1])
      i+=1
i=0
while i <= linenum-1 :
      power=powperc[i]/timeinter
      fileout.write('  <Point Coordinate="'+str(time[i])+'" Value="'+str(power)+'" />'+"\n")
      i+=1

fileout.write('</Pulse>'+'\n')
fileout.close()
