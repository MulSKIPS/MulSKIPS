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

#subprocess.call(['make clean'],cwd=os.getcwd())    
#subprocess.call(['make'],cwd=os.getcwd())    
            

#    words=line.split()
#    time.append(float(words[0]))
#    powperc.append(float(words[1]))
#linenum=len(time)
#print linenum
#timeinter=0.
#i=1
#while i <= linenum-1 :
#      timeinter=timeinter+powperc[i]*(time[i]-time[i-1])
#      i+=1
#i=0
#while i <= linenum-1 :
#      power=powperc[i]/timeinter
#      fileout.write('  <Point Coordinate="'+str(time[i])+'" Value="'+str(power)+'" />'+"\n")
#      i+=1

#fileout.write('</Pulse>'+'\n')
#fileout.close()
#subprocess.call(['/bin/bash', '-c', "gmsh -"+ str(self.__dimension) + " " + self.__GEO_FILENAME + " -o " + self.__MSH_FILENAME + " -format msh2"],
#                        cwd=self.__MESH_DIR)

