# -*- coding: utf-8 -*-
"""
Created on Sat Sep 13 19:59:46 2014

@author: pedro
"""
import os
from libphys import *

dname = '/home/pedro/Downloads/LAB/DATA/2014/Aug/PF/20_08_14/pf02/'

files = []
for file in os.listdir(dname):
    if file.endswith(".bmp"):
        files = np.append(files,file)
files.sort()
print 'Found %d files' %len(files)

for i in range(10,len(files)):
    if (len(files[i])==32):
        os.rename(dname+files[i],dname+files[i][0:3]+'0'+files[i][3:])