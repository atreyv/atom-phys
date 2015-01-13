# -*- coding: utf-8 -*-
"""
Created on Sat Sep 13 19:59:46 2014

@author: pedro
"""
import os
from libphys import *

dname = '/home/pedro/LAB/DATA/2014/Nov/PF/28_11_14/pf01/'

files = []
for file in os.listdir(dname):
    if file.endswith(".bmp"):
        files = np.append(files,file)
files.sort()
print 'Found %d files' %len(files)

j=0
for i in range(1380,1400):
    if (len(files[i])==42):
        if (files[i][11:14]=='0.7'):
            os.rename(dname+files[i],dname+files[i][0:11]+'0.70'+files[i][14:])
            print 'success'
            j+=1
print j