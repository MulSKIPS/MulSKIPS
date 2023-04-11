import sys
import os
import string
import numpy
import xml.etree.ElementTree
import shutil
import time
import math
import subprocess
moddir=os.path.join(os.getcwd(),"modules")
filein=open(moddir+"/defsystem_templ.f","r")
fileout=open(moddir+"/defsystem.f","w")
filebox=open("box.dat","r")
for line in filebox :
    box=line.split()
    lenx=int(box[0])
    leny=int(box[1])
    lenz=int(box[2])
    print(lenx,leny,lenz)
    
#timeinter=100.*1.62799E-7
#fileout.write('<Pulse xmlns="http://www.screen-lasse.com/LIAB/Pulse">'+'\n')
#time=[]
#powperc=[]
ind=0
for line in filein:
    if ind==6 :
        fileout.write('       INTEGER, PARAMETER :: LenX='+str(lenx)+', LenY='+str(leny)+', LenZ='+str(lenz)+"\n")
    else :
        fileout.write(line)
    ind+=1
fileout.close()
filein.close()
filebox.close()

#subprocess.call(['make clean'])    
#subprocess.call(['make'])    
            

