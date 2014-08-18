# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 01:03:49 2013

@author: pedro
"""

# General module import
#from libphys import *

# This program specific modules
from numpy import *
from scipy import optimize
import mayavi.mlab

# Specific functions used/tested in this program
def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*exp(
                 -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = sqrt(abs((arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = sqrt(abs((arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y
    
def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: ravel(gaussian(*p)(*indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

# Define the working directory
dname = '/home/pedro/Downloads/LAB/Nice stuff/paper data and some extra/nice patterns with no parameters 03042013/'

# Define the matrix from intensity 2D profile for FFT. Image normally 1600*1200
matrix_size = 1024
row_start = (1250 - matrix_size)/2
row_end = row_start + matrix_size
col_start = (1650 - matrix_size)/2 #1650 for pump, 1645 for probe
col_end = col_start + matrix_size

# Define the limits where to look for peaks in FT amplitude spectrum
# p0 is pump, 1 to 6 represent the diffracted beams CW
#p_box_width = 20
#p0_min, p0_max = 505, 520
#p1 = [470, 500]
#p2 = [490, 525]
#p3 = [515, 530]
#p4 = [530, 505]
#p5 = [515, 480]
#p6 = [490, 475]

# Define the window size for the FFT plotting
fft_size = 256
fft_ini = matrix_size/2. - fft_size/2.
fft_end = matrix_size/2. + fft_size/2.

# Define the window size for the image plotting
img_size = 750
img_ini = matrix_size/2. - img_size/2.
img_end = matrix_size/2. + img_size/2.


pump_cut = 15000000
noise_cut = 0.

#angle = [] #angles list
#dist = [] #radial distances list
#images = [] #np.zeros(20*120*120).reshape(4*120,5*120)
#normals = []

# counters
file_i = 0 # 

for i in range(16,17):
    plt.clf()
    signal = np.array(plt.imread(dname + 's'+str(i)+'.bmp'), dtype='float64')
    pump = np.array(plt.imread(dname + 'p'+str(i)+'.bmp'), dtype='float64')
#    backg = np.array(plt.imread(dname + 'b'+str(i)+'.bmp'), dtype='float64')
    #normal = Image.open('/home/pedro/Dropbox/PhD at Strathclyde/experiment at INLN/growth_observed_in_pattern/t600.png')
    #normal = np.array(normal)*1.
#   normal = np.sqrt(signal)
    signal = signal[row_start:row_end,col_start:col_end,0]
    pump = pump[row_start:row_end,col_start:col_end,0]
    #normal = signal
    #normal[fft_ini:fft_end,fft_ini:fft_end] = signal[fft_ini:fft_end,fft_ini:fft_end] / pump[fft_ini:fft_end,fft_ini:fft_end]
    normal = (signal[img_ini:img_end,img_ini:img_end]) / (pump[img_ini:img_end,img_ini:img_end])
#    pump = np.sqrt(pump)
#    for x in np.nditer(normal, op_flags=['readwrite']):
#        if (x[...] > pump_cut):
#            x[...] = pump_cut
#        if (x[...] < noise_cut):
#            x[...] = 0.        
    plt.axis('off')
    plt.imshow(normal,interpolation='none',origin='upper',cmap='gray')
    plt.colorbar(shrink=1.,format='%.2f')
    plt.annotate('', xy = (-.02*img_size, 0), xytext=(-.02*img_size,img_size), arrowprops=dict(arrowstyle='<->'), annotation_clip=False)
    plt.text(-.1*img_size,img_size/2, '1 mm', fontsize=32,rotation=90)
    plt.text(10,30,'a)', fontsize=32, color='white')
    plt.savefig(dname + 's'+str(i)+'.tiff', bbox_inches='tight')
#    normals += [normal]
    #normal = normal[300:400,300:400,0]
    #normal = signal[:,:,0]
    #normal = normal[:,:,0]

# Generating Picture for covers
#    mayavi.mlab.clf()
#    probe3 = normal *6
#    pump3 = normal *3
#    probe5=sp.ndimage.gaussian_filter(probe3, sigma=0.7)
#    pump5=sp.ndimage.gaussian_filter(pump3, sigma=0.7)
#    surf(probe5[400:580,410:590],colormap='Reds')
#    surf(pump5[400:580,410:590],colormap='cool',opacity=0.5)
#    mayavi.mlab.savefig('titlepage.tiff', size=(1500,2500), magnification='auto')
#    plt.clf()
#    profile = []
#    for px in range(0,img_size-80):
#        profile += [normal[img_ini+40+int(px*cos(30*pi/180)),img_ini+50+int(px*sin(30*pi/180))]]
#    plt.plot(profile)



#    ff = ft.fft2(normal)
#    ff = ft.fftshift(ff)
#    ff = np.absolute(ff) # same as ff = np.sqrt(np.real(ff)**2 + np.imag(ff)**2)
#    #ff = np.real(ff)
#    #print unravel_index(ff.argmax(), ff.shape)
#    for x in np.nditer(ff, op_flags=['readwrite']):
#        if (x[...] > pump_cut):
#            x[...] = pump_cut
#        if (x[...] < noise_cut):
#            x[...] = 0
##    plt.imsave('/home/pedro/LAB/Nice stuff/paper data and some extra/nice patterns with no parameters 03042013/res'+str(i)+'.png',ff[fft_ini:fft_end,fft_ini:fft_end],format='png',origin='upper')
#    plt.clf()
#    plt.axis('off')
##    plt.imshow(ff[pl_ini:pl_end,pl_ini:pl_end],cmap='jet',norm=LogNorm(), interpolation='none',origin='upper')
#    plt.imshow(ff[fft_ini:fft_end,fft_ini:fft_end],cmap='jet',interpolation='none',origin='upper')
#    plt.colorbar(shrink=1.,format='%1.1e')
#    plt.text(10,30,'b)', fontsize=32, color='white')
#    plt.savefig(dname + 'res'+str(i)+'.pdf', bbox_inches='tight')
#    #plt.imshow(ff[450:570,450:570],origin='upper')
##    images += [ff[fft_ini:fft_end,fft_ini:fft_end]]

##    p0pos = unravel_index(ff[p0_min:p0_max,p0_min:p0_max].argmax(), ff[p0_min:p0_max,p0_min:p0_max].shape)
##    p1pos = unravel_index(ff[p1[0]:p1[0]+p_box_width,p1[1]:p1[1]+p_box_width].argmax(), ff[p1[0]:p1[0]+p_box_width,p1[1]:p1[1]+p_box_width].shape)
##    p2pos = unravel_index(ff[p2[0]:p2[0]+p_box_width,p2[1]:p2[1]+p_box_width].argmax(), ff[p2[0]:p2[0]+p_box_width,p2[1]:p2[1]+p_box_width].shape)
##    p3pos = unravel_index(ff[p3[0]:p3[0]+p_box_width,p3[1]:p3[1]+p_box_width].argmax(), ff[p3[0]:p3[0]+p_box_width,p3[1]:p3[1]+p_box_width].shape)
##    p4pos = unravel_index(ff[p4[0]:p4[0]+p_box_width,p4[1]:p4[1]+p_box_width].argmax(), ff[p4[0]:p4[0]+p_box_width,p4[1]:p4[1]+p_box_width].shape)
#    p5pos2 = unravel_index(ff[p5[0]:p5[0]+p_box_width,p5[1]:p5[1]+p_box_width].argmax(), ff[p5[0]:p5[0]+p_box_width,p5[1]:p5[1]+p_box_width].shape)
##    p6pos = unravel_index(ff[p6[0]:p6[0]+p_box_width,p6[1]:p6[1]+p_box_width].argmax(), ff[p6[0]:p6[0]+p_box_width,p6[1]:p6[1]+p_box_width].shape)
#
#    p0pos = fitgaussian(ff[p0_min:p0_max,p0_min:p0_max])[1:3]
#    p1pos = fitgaussian(ff[p1[0]:p1[0]+p_box_width,p1[1]:p1[1]+p_box_width])[1:3]
#    p2pos = fitgaussian(ff[p2[0]:p2[0]+p_box_width,p2[1]:p2[1]+p_box_width])[1:3]
#    p3pos = fitgaussian(ff[p3[0]:p3[0]+p_box_width,p3[1]:p3[1]+p_box_width])[1:3]
#    p4pos = fitgaussian(ff[p4[0]:p4[0]+p_box_width,p4[1]:p4[1]+p_box_width])[1:3]
#    p5pos = fitgaussian(ff[p5[0]:p5[0]+p_box_width,p5[1]:p5[1]+p_box_width])[1:3]
#    p6pos = fitgaussian(ff[p6[0]:p6[0]+p_box_width,p6[1]:p6[1]+p_box_width])[1:3]
#
#    d0 = [p0_min + p0pos[0], p0_min + p0pos[1]]
#
#    d1 = [d0[0] - (p1[0] + p1pos[0]), p1[1] + p1pos[1] - d0[1]]
#    d2 = [d0[0] - (p2[0] + p2pos[0]), p2[1] + p2pos[1] - d0[1]]
#    d3 = [d0[0] - (p3[0] + p3pos[0]), p3[1] + p3pos[1] - d0[1]]
#    d4 = [d0[0] - (p4[0] + p4pos[0]), p4[1] + p4pos[1] - d0[1]]
#    d5 = [d0[0] - (p5[0] + p5pos[0]), p5[1] + p5pos[1] - d0[1]]
#    d6 = [d0[0] - (p6[0] + p6pos[0]), p6[1] + p6pos[1] - d0[1]]
#
#    angle += [np.arctan(float(d1[0])/d1[1])*180/pi,np.arctan(float(d2[0])/d2[1])*180/pi,np.arctan(float(d3[0])/d3[1])*180/pi,np.arctan(float(d4[0])/d4[1])*180/pi,np.arctan(float(d5[0])/d5[1])*180/pi,np.arctan(float(d6[0])/d6[1])*180/pi]
#    dist += [np.sqrt(d1[0]**2+d1[1]**2),np.sqrt(d2[0]**2+d2[1]**2),np.sqrt(d3[0]**2+d3[1]**2),np.sqrt(d4[0]**2+d4[1]**2),np.sqrt(d5[0]**2+d5[1]**2),np.sqrt(d6[0]**2+d6[1]**2)]
#
#    if(file_i==0):
#        av0 = np.average(normals[file_i][fft_ini:fft_end,fft_ini:fft_end])
#        n0 = normals[file_i][fft_ini:fft_end,fft_ini:fft_end] - av0
#        for j in range(0,len(n0)):
#            for k in range(0,len(n0[j])):
#                if(n0[j,k] < 0):
#                    n0[j,k] = 0.
#        t0 = correlate2d(n0, n0, mode='full')
#    if(file_i>0):
#        print file_i
#        av = np.average(normals[file_i][fft_ini:fft_end,fft_ini:fft_end])
#        n = normals[file_i][fft_ini:fft_end,fft_ini:fft_end] - av
#        for j in range(0,len(n)):
#            for k in range(0,len(n[j])):
#                if(n[j,k] < 0):
#                    n[j,k] = 0.
#        t = correlate2d(n, n0, mode='full')        
#        print np.sqrt((unravel_index(np.argmax(t), shape(t))[0] - unravel_index(np.argmax(t0), shape(t0))[0])**2 + (unravel_index(np.argmax(t), shape(t))[1] - unravel_index(np.argmax(t0), shape(t0))[1])**2)
##    plt.xticks(())
##    plt.yticks(())
#    file_i += 1
##plt.imshow(t0,interpolation='none')
##plt.contour(t0,levels=[10000])
##plt.contour(t,levels=[10000],linestyles='dashed')
##plt.text(10,10,'s11',color='w',fontsize=20)
#
#np.savetxt(dname + 'angle.txt',angle, fmt='%2.3f')
#np.savetxt(dname + 'dist.txt',dist, fmt='%2.1f')
