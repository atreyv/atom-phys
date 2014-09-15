# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 12:32:45 2014

@author: pedro
"""

import os
from libphys import *
dname = '/home/pedro/Downloads/LAB/DATA/2014/Aug/PF/19_08_14/pf02/'
raw_image = 0
frac = 1
half_size = 30

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
signal = np.zeros((2*dy,2*dx))

for i in range(0,10):
    ref = np.array(plt.imread(dname+files[i]),dtype=float64)
    ref = ref[1:]
    ref = ref[centre[0]-dy:centre[0]+dy, centre[1]-dx:centre[1]+dx]
    ref1 += ref / 10.
for i in range(10,11):
    ref = np.array(plt.imread(dname+files[i]),dtype=float64)
    ref = ref[1:]
    ref = ref[centre[0]-dy:centre[0]+dy, centre[1]-dx:centre[1]+dx]
    signal += ref / 1.
centre,dx,dy
#imshow(ref1,interpolation='none')#,norm=LogNorm())
imshow(signal,interpolation='none')#,norm=LogNorm())
plt.axis('off')
#savefig(dname+'red_patterns.png',dpi=100,bbox_inches='tight',frameon=None)
#imsave(dname+'red_patterns.png',signal,dpi=100)
#res = normalize_by_division(signal,ref2)
#imshow(res)

res = signal
res = prepare_for_fft_padding(res)
resft = ft.fft2(res)
resft = ft.fftshift(resft)
resft = np.absolute(resft)
#figure()

#imshow(resft[len(resft)/2-half_size:len(resft)/2+half_size,len(resft)/2-half_size:len(resft)/2+half_size],
#       interpolation='none',norm=LogNorm())
#plt.axis('off')
#savefig(dname+'red_patterns_fft.png',dpi=100,bbox_inches='tight', frameon=False)
#imsave(dname+'red_patterns_fft.png',resft[len(resft)/2-half_size:len(resft)/2+half_size,
#                                          len(resft)/2-half_size:len(resft)/2+half_size],dpi=100,norm=LogNorm())

radial_plot = np.array([])
for i in range(1,int(len(resft)/2)):
    radial_plot = np.append(radial_plot,circle_line_integration(resft,i)[0])
#plt.plot(radial_plot)
argmax(radial_plot)