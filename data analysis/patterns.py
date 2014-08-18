# -*- coding: utf-8 -*-
"""
Created on Sat Aug 16 18:29:09 2014

@author: pedro
"""

import os
from libphys import *
dname = '/home/pedro/Downloads/LAB/DATA/2014/Aug/PF/14_08_14/pf02/'
raw_image = 0
frac = 0.8


files = []
for file in os.listdir(dname):
    if file.endswith(".bmp"):
        files = np.append(files,file)
files.sort()
print 'Found %d files' %len(files)
# Do a fit to probe using a reference image. Also correct for jitter in
# first row from Chameleon
ref = np.array(plt.imread(dname+files[raw_image]),dtype=float64)
ref = ref[1:]
#centre = np.unravel_index(ref.argmax(), ref.shape)
param = fitgaussian2d(ref)
centre = (param[1],param[2])
dx = int(param[4]*frac)
dy = int(param[3]*frac)

ref1 = np.zeros((2*dy,2*dx))
ref2 = np.zeros((2*dy,2*dx))
signal = np.zeros((2*dy,2*dx))

ratio = np.array([])

for i in range(0,10):
    ref = np.array(plt.imread(dname+files[i]),dtype=float64)
    ref = ref[1:]
    ref = ref[centre[0]-dy:centre[0]+dy, centre[1]-dx:centre[1]+dx]
    ref1 += ref / 10.
for i in range(10,20):
    ref = np.array(plt.imread(dname+files[i]),dtype=float64)
    ref = ref[1:]
    ref = ref[centre[0]-dy:centre[0]+dy, centre[1]-dx:centre[1]+dx]
    ref2 += ref / 10.
for i in range(20,390,10):
    for j in range(0,9):
        ref = np.array(plt.imread(dname+files[i]),dtype=float64)
        ref = ref[1:]
        ref = ref[centre[0]-dy:centre[0]+dy, centre[1]-dx:centre[1]+dx]
        signal += ref / 10.
    
        #res = normalize_by_division(signal,ref2)
        res = signal
        res = prepare_for_fft_padding(res)
        
        resft = ft.fft2(res)
        resft = ft.fftshift(resft)
        resft = np.absolute(resft)
        #imshow(resft,interpolation='none',norm=LogNorm())
    print i, 'done'
    
    radial_plot = np.array([])
    for i in range(1,int(len(resft)/2)):
        radial_plot = np.append(radial_plot,circle_line_integration(resft,i)[0])
    p1 = fitgaussian1d(None,radial_plot[20:40])
#    plt.plot(radial_plot)
    p = fitgaussian2d(resft[len(resft)/2-10:len(resft)/2+10,len(resft)/2-10:len(resft)/2+10])
    ratio = np.append(ratio, (p1[0]+p1[3])/(p[0]+p[5]))
plt.plot(ratio)