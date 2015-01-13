# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 16:22:00 2014

@author: ivor
"""
import time
import scipy.ndimage
import matplotlib.patches 
import pickle
from pypeaks import Data, Intervals 
from libphys import *
#from mayavi import mlab
import numpy as np
import matplotlib.pyplot as plt
#import sys
#sys.path.insert(1, r'./../functions')  # add to pythonpath
#from detect_peaks import detect_peaks
import os

image=[]
im ='/home/pedro/Downloads/LAB/DATA/2014/Sept/10_09_14/pf01/pump_2014-09-10-163807-0000.bmp'
image=np.array(plt.imread(im),dtype=float64)

plt.figure(figsize=(10,10))
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)

t = time.time() #Start timing the peak-finding algorithm
threshold=200

im_thresh = np.copy(image)
im_thresh[im_thresh<threshold] = 0
print 'Time after thresholding: %.5f seconds'%(time.time()-t)

#now find the objects
labeled_image, number_of_objects = scipy.ndimage.label(im_thresh)
print 'Time after labeling: %.5f seconds'%(time.time()-t)

peak_slices = scipy.ndimage.find_objects(labeled_image)
print 'Time after finding objects: %.5f seconds'%(time.time()-t)

def centroid(data):
    h,w = np.shape(data)   
    x = np.arange(0,w)
    y = np.arange(0,h)

    X,Y = np.meshgrid(x,y)

    cx = np.sum(X*data)/np.sum(data)
    cy = np.sum(Y*data)/np.sum(data)

    return cx,cy

centroids = []

for peak_slice in peak_slices:
    dy,dx  = peak_slice
    x,y = dx.start, dy.start
    cx,cy = centroid(im_thresh[peak_slice])
    centroids.append((x+cx,y+cy))

print 'Total time: %.5f seconds\n'%(time.time()-t)

###########################################
#Now make the plots:
for ax in (ax1,ax2,ax3,ax4): ax.clear()
ax1.set_title('Original image')
ax1.imshow(image,origin='lower')

ax2.set_title('Thresholded image')
ax2.imshow(im_thresh,origin='lower')

ax3.set_title('Labeled image')
ax3.imshow(labeled_image,origin='lower') #display the color-coded regions

for peak_slice in peak_slices:  #Draw some rectangles around the objects
    dy,dx  = peak_slice
    xy     = (dx.start, dy.start)
    width  = (dx.stop - dx.start + 1)
    height = (dy.stop - dy.start + 1)
    rect = matplotlib.patches.Rectangle(xy,width,height,fc='none',ec='red')
    ax3.add_patch(rect,)

ax4.set_title('Centroids on original image')
ax4.imshow(image,origin='lower')

for x,y in centroids:
    ax4.plot(x,y,'kx',ms=10)
#ax4.set_xlim(0,size)
#ax4.set_ylim(0,size)

plt.tight_layout
plt.show()