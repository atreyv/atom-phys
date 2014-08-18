# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 10:14:02 2014

@author: pedro
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 01:03:49 2013

@author: pedro
"""

# General module import
from libphys import *

# This program specific modules

# Specific functions used/tested in this program


# Define the working directory
dname = '/home/pedro/Dropbox/PhD at Strathclyde/experiment at INLN/Guillaume exp/newexperiments-20140402/'

# Initial conditions
plt.clf()

# Define the images size
#imgw = 1600
#imgh = 1200

# Define the matrix from intensity 2D profile for FFT. Image normally 1600*1200
#matrix_size = 1024
#row_start = (1200 - matrix_size)/2
#row_end = row_start + matrix_size
#col_start = (1600 - matrix_size)/2
#col_end = col_start + matrix_size

# Define the window size for the FFT plot
#fft_size = 512
#fft_ini = matrix_size/2. - fft_size/2.
#fft_end = matrix_size/2. + fft_size/2.


# Define parameters

fft_size = 256
fit_window_radius = 10
pump_cut = 80000000
noise_cut = -1

pump_peaks = np.array([])
q1_peaks = np.array([])
q2_peaks = np.array([])
#pump = np.array(plt.imread(dname + 'pump.bmp'), dtype='float64')
#pump = pump[:,:,0]
#pumpff = prepare_for_fft(pump,fft_size,0)
#centre = unravel_index(pump.argmax(), pump.shape)    

for image in range(0,1):
#    plt.clf()
    print image
    signal = np.array(plt.imread(dname + 'optomeca_blue.png'), dtype='float64')
    signal = signal[:,:,0]
#    pump = np.array(plt.imread(dname + 'p'+str(image)+'.bmp'), dtype='float64')
#    pump = pump[:,:,0]
#    backg = np.array(plt.imread(dname + 'b'+str(image)+'.bmp'), dtype='float64')
#    backg = backg[:,:,0]
#    centre = unravel_index(pump.argmax(), pump.shape)
#    centre = unravel_index(signal.argmax(), signal.shape)

#    pumpff = prepare_for_fft(pump,fft_size,0)    
#    signalff = prepare_for_fft(signal,fft_size,centre)
    signalff=signal
#    signalff = normalize_by_divison(signalff,pumpff)
#    signalff = np.sqrt(signalff)
#    signalff = normalize_by_divison(signalff,signalff)
#    plt.axis('off')
#    plt.imshow(signalff,interpolation='none',cmap='gray',origin='upper')
#    plt.imshow(pumpff,interpolation='none',norm=LogNorm(),cmap='jet',origin='upper')
#    plt.savefig(dname + 's'+str(image)+'.eps', bbox_inches='tight')

    ff = ft.fft2(signalff)
    ff = ft.fftshift(ff)
    ff = np.absolute(ff) # same as ff = np.sqrt(np.real(ff)**2 + np.imag(ff)**2)
#    ff=ff**2
#    ff = np.real(ff)
#    for x in np.nditer(ff, op_flags=['readwrite']):
#        if (x[...] > pump_cut):
#            x[...] = pump_cut
#        if (x[...] < noise_cut):
#            x[...] = noise_cut
#    plt.clf()
#    plt.colorbar(shrink=.92)
#    plt.savefig(dname + 'res'+str(image)+'.eps', bbox_inches='tight')
#    plt.axis('off')
#    circle_line_integration(ff,50)
#    plt.imshow(ff,cmap='jet', interpolation='none',norm=LogNorm(),origin='upper')
#    plt.imshow(ff,cmap='jet', interpolation='none',origin='upper')
    ff_radial = np.array([])
    pix_count = np.array([])
    ff_radial_av = np.array([])
    for i in range(0,int(len(ff)/2)):
        ff_radial = np.append(ff_radial,circle_line_integration(ff,i)[0])
        pix_count = np.append(pix_count,circle_line_integration(ff,i)[1])
    ff_radial_av = ff_radial / pix_count
    plt.scatter(np.arange(len(ff)/2),ff_radial_av)
    plt.plot(np.arange(len(ff)/2),ff_radial_av)
#    print pix_count.sum()
    ff_radial_copy = ff_radial_av
#    # 1
#    pump_peak = fitgaussian1d(ff_radial_copy[:fit_window_radius])
#    pump_peaks = np.append(pump_peaks, pump_peak)
#    g=gaussian1d(*pump_peak)
#    x=np.linspace(0,fit_window_radius,fit_window_radius*10)
#    plt.plot(x,g(x))
#    # 2
#    pump_peak = fitgaussian1d(ff_radial_copy[2*fit_window_radius:4*fit_window_radius])
#    pump_peak[1] += 2*fit_window_radius
#    q1_peaks = np.append(q1_peaks, pump_peak)
#    g=gaussian1d(*pump_peak)
#    x=np.linspace(2*fit_window_radius,4*fit_window_radius,fit_window_radius*10)
#    plt.plot(x,g(x))
#    # 3
#    pump_peak = fitgaussian1d(ff_radial_copy[4*fit_window_radius:5.5*fit_window_radius])
#    pump_peak[1] += 4*fit_window_radius
#    q2_peaks = np.append(q2_peaks, pump_peak)
#    g=gaussian1d(*pump_peak)
#    x=np.linspace(4*fit_window_radius,5.5*fit_window_radius,fit_window_radius*10)
#    plt.plot(x,g(x))
#pump_peaks = pump_peaks.reshape((image+1,pump_peak.size))
#q1_peaks = q1_peaks.reshape((image+1,pump_peak.size))
#q2_peaks = q2_peaks.reshape((image+1,pump_peak.size))
    
