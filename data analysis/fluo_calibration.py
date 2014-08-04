# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 17:25:50 2014

@author: lab
"""

from libphys import *
#import mayavi.mlab
import os

caldname = 'Z:/Data/2014/July/MOT/25_07_14/fluo_calibration/'
averages = 10
raw_image = 1
pixel_size = .000375 # centimeters per pixel for the bare camera without telescope
power = np.array([0.000385, 0.00062, 0.00101, 0.00149, 0.0017, 0.0019, 0.00221, 0.00247, 0.00297, 0.00338, 0.00365, 0.004]) #errorbars 0.2 uW 

def load_and_average_cal(caldname,raw_image,averages):
    files = []
    for file in os.listdir(caldname):
        if file.endswith(".bmp"):
            files = np.append(files,file)
    files.sort()
    print 'Found %d files' %len(files)
    
    ref = np.array(plt.imread(caldname+files[raw_image]),dtype=float64)
    ref = ref[1:]
    #centre = np.unravel_index(ref.argmax(), ref.shape)
    param = fitgaussian2d(ref)
    centre = (param[1],param[2])
    dx = int(param[4]*2.9) 
    dy = int(param[3]*2.9)
    res = np.zeros((2*dy, 2*dx))
    Data = np.array([res],dtype=float64) 
    
    for i in range(0,len(files),averages):
#        pattern2 = files[i][len(prefix):len(prefix)+number_of
        for j in range(0,averages,1):
            # Load images for each iteration
            signal = np.array(plt.imread(caldname+files[j+i]),dtype=float64)
            #ref = np.array(plt.imread(dname+files[j+i+1]),dtype=float64)
            #backg = np.array(plt.imread(dname+files[j+i+1]),dtype=float64)
            # Correct for some weird first row jitter from Chameleon    
            signal = signal[1:]
            #ref = ref[1:]
            #backg = backg[1:]
            # Crop around probe beam
            signal = signal[centre[0]-dy:centre[0]+dy,centre[1]-dx:centre[1]+dx]
            #ref = ref[centre[0]-dy:centre[0]+dy,centre[1]-dx:centre[1]+dx]
            #backg = backg[centre[0]-dy_sigma:centre[0]+dy_sigma,centre[1]-dx_sigma:centre[1]+dx_sigma]
            # Sum and average
            res += signal / averages
            print "Signal image", j
        # Load final matrix and reset res
        Data = np.append(Data,[res],axis=0)
        res = np.zeros((2*dy, 2*dx))
        print "Phase", i/(j+1)+1 , "complete"
    # Clean the Data points matrix of zeros in the first slot 
    Data = np.delete(Data,0,axis=0)
    
    return Data, dx/2.5, dy/2.5
    
def calibration(Data,power):
    
    power_pixel = np.zeros(len(power))
    intensity_pixel = np.zeros(len(power))
    sig_x=np.zeros(len(power))
    sig_y=np.zeros(len(power))
    intensity=np.zeros(len(power)) 
    
    for i in range(0,len(power),1):
        
        p = fitgaussian2d(Data[i])      
        centre = (p[1],p[2])
        #get the beamsize for different points (powers)
        sig_y[i] = np.abs(p[3])
        sig_x[i] = np.abs(p[4])
        #get the power and intensities in pixel units, sumation over 4*sigma to 
        #get the power (should be over 6*sigma but image is too small)
        power_pixel[i]=np.sum(Data[i][centre[0]-2.5*sig_y[i]:centre[0]+2.5*sig_y[i],centre[1]-2.5*sig_x[i]:centre[1]+2.5*sig_x[i]])
        power_pixel[i]=power_pixel[i]+0.05*power_pixel[i]
        intensity_pixel[i]=0.66*power_pixel[i]/(4*sig_x[i]*sig_y[i])
        #get the intensity in SI units for each point
        intensity[i]=0.66*power[i]/(4*sig_x[i]*sig_y[i]*pixel_size**2) 
        print "sig_x = %.2f"%(sig_x[i])
        
    #calculate fitting parameters
    powerfit=np.polyfit(power_pixel,power,1)
    intensityfit=np.polyfit(intensity_pixel,intensity,1)
    
    #plot the points and fitted line on two windows    
    
    #power
    figure()
    plt.scatter(power_pixel,power,label='Power')
    plt.xlabel('P')
    plt.ylabel('P (mW)')
    x_pow=linspace(power_pixel[0],power_pixel[len(power)-1],num=1000)
    plt.plot(x_pow,powerfit[1]+powerfit[0]*x_pow)
    plt.text(power_pixel[len(power)-12],power[len(power)-2],'eta_pow=%.2e mW/pixelheight (200us, gain=0)'%powerfit[0])    
    plt.legend(loc=0)
    savepath_pow = os.path.join(caldname, 'calibration_pow.pdf')  
    plt.savefig(savepath_pow) 
    
    #intensity
    figure()
    x_int=linspace(intensity_pixel[0],intensity_pixel[len(intensity)-1],num=1000) 
    plt.scatter(intensity_pixel,intensity,label='Intensity')
    plt.plot(x_int,intensityfit[1]+intensityfit[0]*x_int)
    plt.text(intensity_pixel[len(intensity)-12],intensity[len(intensity)-2],'eta_int=%.2e (mW/cm^2)/(pixelheight/pixelsize)'%intensityfit[0])
    plt.legend(loc=0)
    plt.xlabel('I')
    plt.ylabel('I (mW/cm^2)')
    savepath_int = os.path.join(caldname, 'calibration_int.pdf')    
    plt.savefig(savepath_int)    

    print "eta_power = %.12f mW/pixelheight (200 us, gain=0)"%powerfit[0]
    print "\neta_intensity = %f (mW/cm^2)/(pixelheight/pixelsize)"%intensityfit[0]

    return powerfit, intensityfit  