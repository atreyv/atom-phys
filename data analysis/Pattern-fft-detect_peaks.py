# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 01:03:49 2013

@author: pedro
"""

# General module import
from libphys import *

# This program specific modules
matplotlib.rc('font', size=16)


# Specific functions used/tested in this program


# Define the working directory
dname = '/home/pedro/Downloads/LAB/Nice stuff/paper data and some extra/17102012_contrast vs intensity and detuning/'

# Initial conditions
fig = figure()
points = np.zeros(103, dtype='int64')
a = np.arange(11-1,21)
points[a] = a+1
a = np.arange(31-1,34)
points[a] = a+1
a = np.arange(38-1,46)
points[a] = a+1
a = np.arange(68-1,74)
points[a] = a+1
a = np.arange(100-1,103)
points[a] = a+1
full_list = np.arange(1,104)

#int_det_file = np.loadtxt(dname+'17102012.dat',skiprows=2,comments='--')
int_det_file = np.loadtxt(dname+'int_det.dat')
int_det_file = int_det_file[full_list!=points]
#int_det_file = int_det_file[np.argsort(int_det_file[:,2])]
images = full_list[full_list!=points]

#det = int_det_file[0,2]
#images = images[0:10]
#int_det = int_det_file[0:10,1]

#det = int_det_file[10,2]
#images = images[10:19]
#int_det = int_det_file[10:19,1]

det = int_det_file[22,2]
images = images[22:32]
int_det = int_det_file[22:32,1]

#det = int_det_file[32,2]
#images = images[32:43]
#int_det = int_det_file[32:43,1]

#det = int_det_file[43,2]
#images = images[43:56]
#int_det = int_det_file[43:56,1]

#det = int_det_file[56,2]
#images = images[56:68]
#int_det = int_det_file[56:68,1]


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

fft_size = 400
fit_window_i = 15
fit_window_f = 45
fit_window_radius = fit_window_f - fit_window_i
pixel_size = 8.8 # um
#pump_cut = 80000000
noise_cut = -1
a = fft_size*pixel_size # um
fourier_pixel = 1./a # um^-1
theta_unit = 0.78 * fourier_pixel *1000 # mrad
def Lambda_cal(theta):
    return 0.78 / theta *1000


pump_peaks = np.array([])
q1_peaks = np.array([])
q2_peaks = np.array([])
pump = np.array(plt.imread(dname + 'pump.bmp'), dtype='float64')
pump = pump[:,:,0]
pumpff = prepare_for_fft(pump,fft_size,0)
centre = unravel_index(pump.argmax(), pump.shape)    

for image in images[6:7]:
#    plt.clf()
    print image
    signal = np.array(plt.imread(dname + 's'+str(image)+'.bmp'), dtype='float64')
    signal = signal[:,:,0]
#    pump = np.array(plt.imread(dname + 'p'+str(image)+'.bmp'), dtype='float64')
#    pump = pump[:,:,0]
#    backg = np.array(plt.imread(dname + 'b'+str(image)+'.bmp'), dtype='float64')
#    backg = backg[:,:,0]
#    centre = unravel_index(pump.argmax(), pump.shape)

#    pumpff = prepare_for_fft(pump,fft_size,0)    
    signalff = prepare_for_fft(signal,fft_size,centre)
#    signalff = normalize_by_divison(signalff,pumpff)
#    signalff = np.sqrt(signalff)
#    signalff = normalize_by_divison(signalff,signalff)
#    plt.axis('off')
#    plt.imshow(signalff,interpolation='none',norm=LogNorm(),cmap='jet',origin='upper')
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
#    plt.imshow(ff[150:250,150:250],cmap='jet', interpolation='none',origin='upper',vmax=20000)

    ff_radial = np.array([])
    pix_count = np.array([])
    ff_radial_av = np.array([])
    for i in range(1,int(len(ff)/2)):
        ff_radial = np.append(ff_radial,circle_line_integration(ff,i)[0])
        pix_count = np.append(pix_count,circle_line_integration(ff,i)[1])
    ff_radial_av = ff_radial# / pix_count
    ylabel('Normalized amplitude',fontsize=24)    
    xlabel(r'$\mathrm{\theta}$ (mrad)',fontsize=24) 
    plt.plot(ff_radial_av[0:],color='green')
##    print pix_count.sum()
    
    # 1
    ff_centre = np.unravel_index(ff.argmax(), ff.shape)
    pump_peak = fitgaussian2d(ff[ff_centre[0]-fit_window_radius:ff_centre[0]+fit_window_radius,\
                              ff_centre[0]-fit_window_radius:ff_centre[0]+fit_window_radius])
#    pump_peak = np.amax(ff_radial_av[0:5])
##    pump_peak[1] -= fit_window_radius
    pump_peaks = np.append(pump_peaks, pump_peak)
#    g=gaussian1d(pump_peak[0],pump_peak[2],pump_peak[4],pump_peak[5])
#    x=np.linspace(0,2*fit_window_radius,2*fit_window_radius*100)
#    plt.plot(x,g(x))

    # 2
    q1_peak = fitgaussian1d(None,ff_radial_av[fit_window_i:fit_window_f])
    q1_peak[1] += fit_window_i
    g=gaussian1d(*q1_peak)    
#    q1_peak[1] *= theta_unit
    q1_peaks = np.append(q1_peaks, q1_peak)
    x=np.linspace(fit_window_i,fit_window_f,fit_window_radius*100)
#    plt.plot(x,g(x))

    # 3
#    pump_peak = fitgaussian1d(ff_radial_copy[4*fit_window_radius:6*fit_window_radius])
#    pump_peak[1] += 4*fit_window_radius
#    q2_peaks = np.append(q2_peaks, pump_peak)
#    g=gaussian1d(*pump_peak)
#    x=np.linspace(4*fit_window_radius,5.5*fit_window_radius,fit_window_radius*10)
#    plt.plot(x,g(x))
   

#pump_peaks = pump_peaks.reshape((len(images),pump_peak.size))
#q1_peaks = q1_peaks.reshape((len(images),q1_peak.size))
##q2_peaks = q2_peaks.reshape((image,pump_peak.size))
#
## Plot 4 graphs
#ax0 = fig.add_subplot(211)
#ax1 = fig.add_subplot(212)
##ax2 = fig.add_subplot(223)
##ax3 = fig.add_subplot(224)
#ax0.text(500,0.8,'$\Delta$ = %1.1f $\Gamma$'%det)
#
##x0 = np.linspace(1,len(q1_peaks[:,0])+1,len(q1_peaks[:,0]))
#ax0.set_ylabel('Normalized energy')
#ax0.set_xlabel('I (mW/cm$^2$)')
##q1_peaks[:,0]=q1_peaks[:,0]/np.average(pump_peaks[:,0])
#q1_peaks_height = (q1_peaks[:,0]+q1_peaks[:,3]) / (pump_peaks[:,0]+pump_peaks[:,5])
#ax0.scatter(int_det,q1_peaks_height,color='blue')
#
##x1 = np.linspace(1,len(q1_peaks[:,1])+1,len(q1_peaks[:,1]))
#ax1.set_ylabel('$\Lambda$ ($\mathrm{\mu}$m)')
#ax1.set_xlabel('I (mW/cm$^2$)')
##q1_peaks[:,1]=Lambda_cal(q1_peaks[:,1]*theta_unit)
#ax1.scatter(int_det,q1_peaks[:,1],color='red')
#


#x2 = np.linspace(1,len(q1_peaks[:,2])+1,len(q1_peaks[:,2]))
#ax2.set_ylabel('hwhm')
#ax2.set_xlabel('I (mW/cm$^2$)')
#q1_peaks[:,2]=q1_peaks[:,2]/q1_peaks[9,2]
#ax2.scatter(int_det,q1_peaks[:,2],color='blue')
#
##x3 = np.linspace(1,len(q1_peaks[:,3])+1,len(q1_peaks[:,3]))
#ax3.set_ylabel('background')
#ax3.set_xlabel('I (mW/cm$^2$)')
#q1_peaks[:,3]=q1_peaks[:,3]/q1_peaks[9,3]
#ax3.scatter(int_det,q1_peaks[:,3],color='blue')








#    plt.imshow(ff,cmap='jet', interpolation='none',norm=LogNorm(),origin='upper')    
    
    
#    if(i==0):
#        av0 = np.average(normals[i][fft_ini:fft_end,fft_ini:fft_end])
#        n0 = normals[i][fft_ini:fft_end,fft_ini:fft_end] - av0
#        for j in range(0,len(n0)):
#            for k in range(0,len(n0[j])):
#                if(n0[j,k] < 0):
#                    n0[j,k] = 0.
#        t0 = correlate2d(n0, n0, mode='full')
#    if( i>0):
#        print i
#        av = np.average(normals[i][fft_ini:fft_end,fft_ini:fft_end])
#        n = normals[i][fft_ini:fft_end,fft_ini:fft_end] - av
#        for j in range(0,len(n)):
#            for k in range(0,len(n[j])):
#                if(n[j,k] < 0):
#                    n[j,k] = 0.
#        t = correlate2d(n, n0, mode='full')        
#        print np.sqrt((unravel_index(np.argmax(t), shape(t))[0] - unravel_index(np.argmax(t0), shape(t0))[0])**2 + (unravel_index(np.argmax(t), shape(t))[1] - unravel_index(np.argmax(t0), shape(t0))[1])**2)
##    plt.xticks(())
##    plt.yticks(())
