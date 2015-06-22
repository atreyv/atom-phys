# -*- coding: utf-8 -*-
"""
Created on Sat Sep 13 19:59:46 2014

@author: pedro
"""
import os
from libphys import *

for k in range(14,15):
    dname = '/home/pedro/LAB/Nice stuff/paper data and some extra/Fig 3a/'+str(k)+'0ms/'
    
    files = []
    for file in os.listdir(dname):
        if file.endswith(".bmp"):
            files = np.append(files,file)
    files.sort()
    print 'Found %d files' %len(files)
    
    j=0
    for i in range(0,10):
        if (len(files[i])>=0):
            if (files[i][0:1]=='s'):
                os.rename(dname+files[i],dname+files[i][0:1]+str(k)+'0-'+files[i][1:])
                print 'success'
                j+=1
    print j