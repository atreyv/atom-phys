# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 19:20:50 2014

@author: pedro
"""

from pylab import *
import scipy as sp
from scipy import ndimage
from libphys import *
mpl.rcParams['mathtext.fontset'] = 'stix'

dpump = '/home/pedro/Downloads/LAB/Nice stuff/paper data and some extra/Fig 2 data 20112012/pump/'
dprobe = '/home/pedro/Downloads/LAB/Nice stuff/paper data and some extra/Fig 2 data 20112012/probe/'
number = '8'
signalpump = np.array(plt.imread(dpump + 's' + number + '.bmp'), dtype='float64')
pumppump = np.array(plt.imread(dpump + 'p' + number + '.bmp'), dtype='float64')
signalprobe = np.array(plt.imread(dprobe + 's' + number + '.bmp'), dtype='float64')
pumpprobe = np.array(plt.imread(dprobe + 'p' + number + '.bmp'), dtype='float64')

pixel_size = 4.4 #micrometer
I_sat = 3.58
aI_0 = 129
det = 7

plot_angle = -27.5
img_slice_y = 162
img_slice_xi = 115#95
img_slice_xf = 215#270

# Define the matrix from intensity 2D profile for FFT. Image normally 1600*1200
matrix_size = 1024
row_start = (1100 - matrix_size)/2
row_end = row_start + matrix_size
col_start = (1650 - matrix_size)/2 #1650 for pump, 1645 for probe
col_end = col_start + matrix_size
col_start_probe = (1645 - matrix_size)/2

# Define the window size for the image plotting
img_size = 256
img_ini = matrix_size/2 - img_size/2
img_end = matrix_size/2 + img_size/2

signalpump = signalpump[row_start:row_end,col_start:col_end,0]
pumppump = pumppump[row_start:row_end,col_start:col_end,0]

signalprobe = signalprobe[row_start:row_end,col_start_probe:col_end,0]
pumpprobe = pumpprobe[row_start:row_end,col_start_probe:col_end,0]

npump = signalpump[img_ini:img_end,img_ini:img_end] / pumppump[img_ini:img_end,img_ini:img_end]
npump = ndimage.rotate(npump,plot_angle)
npump = ndimage.gaussian_filter(npump,sigma=.0)

nprobe = signalprobe[img_ini:img_end,img_ini:img_end] / pumpprobe[img_ini:img_end,img_ini:img_end]
nprobe = ndimage.rotate(nprobe,plot_angle)
nprobe = ndimage.gaussian_filter(nprobe,sigma=.0)

plt.close()
imsave(dpump+'npump8.tif',npump,origin='upper',cmap='gray')
imsave(dprobe+'nprobe8.tif',nprobe,origin='upper',cmap='gray')

fig = figure()
#fig.subplots_adjust(left=0.1,right=0.9)
ax1 = fig.add_subplot(211)
ax1.set_ylabel('Transmission\nof pump',fontsize=16)
ax2 = fig.add_subplot(212)
ax2.set_xlabel(r'$x$ ($\mathrm{\mu m}$)',fontsize=18)#/\Lambda$")
ax2.set_ylabel('Transmission\nof probe',fontsize=16)

x = np.linspace(0,(img_slice_xf-img_slice_xi)*pixel_size,img_slice_xf-img_slice_xi)#/87
y1 = npump[img_slice_y,img_slice_xi:img_slice_xf]# * aI_0 / I_sat * 1 / (1+4*det**2)
ax1.plot(x,y1)
y2 = nprobe[img_slice_y,img_slice_xi:img_slice_xf]
#y2normalisation = x[len(x)-1]/(pixel_size*y2.sum())
#y2 = y2normalisation * nprobe[img_slice_y,img_slice_xi:img_slice_xf]
ax2.plot(x,y2)
#ax1.plot(pump[0],pump[1])
#savefig('experimental_intensity_density.eps',bbox_inches='tight')
