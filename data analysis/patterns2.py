# -*- coding: utf-8 -*-
"""
Created on Wed Dec 10 10:34:08 2014

@author: pedro
"""

#%matplotlib inline
import os
from libphys import *
from numpy import *
from matplotlib import animation
mpl.rcParams['figure.figsize'] = (16.0, 6.0)
mpl.rc('font', size=16)

dname = '/home/pedro/LAB/DATA/2014/Nov/PF/28_11_14/pf01/'
files = load_files(dname,".bmp")

pixel_size = 3.75 # micrometers
magnification = 2.
raw_image = 1361
frac = 1
half_size = 35

### Do a fit to a gaussian envelope using a reference image (can have patterns). Also correct for jitter in first row from Chameleon.

param = use_ref_to_locate_centre(dname+files[raw_image])
refft, ref = do_fft_with_ref(dname+files[raw_image],param,frac)
fft_size = np.shape(ref)[0]
int_total = ref.sum()
fft_size,int_total,plt.imshow(ref)

mpl.rcParams['figure.figsize'] = (16.0, 6.0)
peaks_pos = np.array([])
peaks_values = np.array([])
B = np.array([])
int_channel = np.array([])
for i in range(0,10):
    resft, signal = do_fft_with_ref(dname+files[i],param,frac)
    fig = plt.figure()
    plot1 = fig.add_subplot(121)
    plot2 = fig.add_subplot(122)
    imshowfft(plot1,resft,0.2)
    plot1.text(5,5,'B=%.2f'%float(files[i][9:13]),fontsize=20,color='black')
    radial_plot = np.array([])
    for j in range(1,int(len(resft)/5)):
        radial_plot = np.append(radial_plot,circle_line_integration(resft,j)[0])
    correction = gets_integration_noise_on_fourier_space(radial_plot,start_pos=len(radial_plot)/2)
    radial_plot -= correction(np.arange(0,len(radial_plot)))
    #radial_plot = remove_backg_on_fourier_space(radial_plot,start_pos=len(radial_plot)/2)
    plot2.plot(radial_plot)
    #plot2.imshow(signal,interpolation='none',origin='upper')
    peaks_temp = find_peaks(radial_plot,1000,20,plot=True)
    if (peaks_temp != 0):
        sorting = np.argsort(peaks_temp.peaks['peaks'][0])
        if (peaks_temp.peaks['peaks'][1][sorting[0]] > 0.04):
            pos_guess = peaks_temp.peaks['peaks'][0][sorting[0]]
            pos_actual = np.argmax(radial_plot[pos_guess-5:pos_guess+5]) + pos_guess-5
            #correction = gets_integration_noise_on_fourier_space(radial_plot,start_pos=len(radial_plot)/2)
            Lambda = 1. / (pos_actual / (pixel_size/magnification*fft_size))
            peaks_pos = np.append(peaks_pos, Lambda)
            peaks_values = np.append(peaks_values,(np.amax(radial_plot[pos_actual])#-correction(pos_actual)
                                                   )/radial_plot[0])
#            peaks_values = np.append(peaks_values,peaks_temp.peaks['peaks'][1][sorting[0]])
            B = np.append(B,(float(files[i][9:13])-0.63) * 1000 * 3)
            int_channel = np.append(int_channel, signal.sum()/int_total)

value_mean = np.array([])
value_errors = np.array([])
pos_mean = np.array([])
pos_errors = np.array([])
int_mean = np.array([])
int_errors = np.array([])
B_value = np.array([])
j = np.inf
for i in B:
    if i != j:
        temp = B==i
        B_value = np.append(B_value,i)
        value_mean = np.append(value_mean, np.average(peaks_values[temp==True]))
        value_errors = np.append(value_errors, np.std(peaks_values[temp==True]))
        pos_mean = np.append(pos_mean, np.average(peaks_pos[temp==True]))
        pos_errors = np.append(pos_errors, np.std(peaks_pos[temp==True]))
        int_mean = np.append(int_mean, np.average(int_channel[temp==True]))
        int_errors = np.append(int_errors, np.std(int_channel[temp==True]))
        j = i
    else:
        pass

mpl.rcParams['figure.figsize'] = (16.0, 16.0)
fig = plt.figure()
plot1 = fig.add_subplot(221)
plot2 = fig.add_subplot(222)
plot3 = fig.add_subplot(224)
plot1.plot(B_value,pos_mean,linewidth=0,marker='o')
plot1.errorbar(B_value,pos_mean,yerr=pos_errors,ls='none')
plot1.set_xlabel(r'$B_{transverse}$ (mG)')
plot1.set_xlim(xmin=-400,xmax=400)
plot1.set_ylabel(r'$\Lambda$ ($\mathrm{\mu}$m)')
#Bchar = ['normal', 'y 1.9G', 'y 5.7G', 'z 1.9G', 'z 5.7G', 'z 11.4G']
#plt.xticks(range(1,len(Bchar)+1), Bchar)
plot2.plot(B_value,value_mean,linewidth=0,marker='o')
plot2.errorbar(B_value,value_mean,yerr=value_errors,ls='none')
plot2.set_xlim(xmin=-400,xmax=400)
plot2.set_xlabel(r'$B_{transverse}$ (mG)')
plot2.set_ylabel(r'$A_{q_1}$/$A_{q_0}$')
plot3.plot(B_value,int_mean,linewidth=0,marker='o')
plot3.errorbar(B_value,int_mean,yerr=int_errors,ls='none')
plot3.set_xlim(xmin=-400,xmax=400)
plot3.set_xlabel(r'$B_{transverse}$ (mG)')
plot3.set_ylabel(r'transmission')
plot3.set_ylim(ymin=0)
#fig.savefig(dname+'ortho-postmirror-pol.pdf')

np.savetxt(dname+'pattern_strength_parallel_post-pol.dat',np.transpose(np.array([B_value,value_mean,value_errors])),
           fmt='%1.f %.4f %.4f', header='B(mG) Aq1/Aq0 std', comments='')