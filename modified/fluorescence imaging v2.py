# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 14:05:59 2014

@author: lab
"""

from numpy import *
from libphys import *
import os


# Parameters to set for load_and_average
dname = '/home/pedro/LAB/DATA/2014/Nov/MOT/25_11_14/tof01/'
#prefix = 'tof'
#number_of_digits_for_files = 1
cycles = 10 # number of averages
phases = 13 # number of interesting cases
# Image to use for 2D gaussian fit. From this fit extract the gaussian zero
# and FWHM. FWHM is used as size of frame for operation: (sig-backg)/(ref-backg)
raw_image = 240

# Parameters for time of glight:
time = np.array([2,3,4,5,6,7,8,9,10,11,12,13,14])
#time = np.array([2,3,4,5,6,7,8,9,10,11,12,13,14])
#time = np.array([2,3,4,5,6,7,8,9,10])

pixel_size = 17.4 # microm
atomic_mass = 87


# FWHM fraction of gaussian fit to a probe profile to use for integration
# of transmission intensity 
fraction = .4
# On resonace unsorted data point (default = 5)
unsorted_resonance_point = 3
#
#
#
def load_and_average_fluo(dname,raw_image,cycles):
    ## Load all file names in folder
    files = []
    for file in os.listdir(dname):
        if file.endswith(".bmp"):
            files = np.append(files,file)
    files.sort()
    print 'Found %d files' %len(files)
#        
#    for i in range(1,n_cycles+1):
#        elem=i*n_avg-1
#        cond_list.append(elem)
    # Do a fit to probe using a reference image. Also correct for jitter in
    # first row from Chameleon
    ref = np.array(plt.imread(dname+files[raw_image]),dtype=float64)
    ref = ref[1:]
    #centre = np.unravel_index(ref.argmax(), ref.shape)
    param = fitgaussian2d(ref)
    centre = (param[1],param[2])
    # Choose a size for the images good to perform division
    # Choose either the waist of the probe or the FWHM
    # The waist:
    #dx = int(gaussian_hwhm_to_radius(param[4]))
    #dy = int(gaussian_hwhm_to_radius(param[3]))
    # The FWHM:
    dx = int(param[4]) *3
    dy = int(param[3]) *3
    
    # Start the count, assign filename pattern and set the containers
    # for data: res for loading each iteration as transmission
    # and Data to load each data point.
   # count = 0
   # pattern = files[0][len(prefix):len(prefix)+number_of_digits_for_files]
    print "dy",dy
    print "dx",dx
    res = np.zeros((2*dy, 2*dx))
    Data = np.array([res],dtype=float64)
    Data_ind = np.array([res],dtype=float64)
    
    # Start main loop
    for i in range(0,len(files),cycles*2):
#        pattern2 = files[i][len(prefix):len(prefix)+number_of_digits_for_files]
        # print pattern2
        for j in range(0,cycles*2,2):
            # Load images for each iteration
            signal = np.array(plt.imread(dname+files[j+i]),dtype=float64)
            #ref = np.array(plt.imread(dname+files[j+i+1]),dtype=float64)
            backg = np.array(plt.imread(dname+files[j+i+1]),dtype=float64)
            # Correct for some weird first row jitter from Chameleon    
            signal = signal[1:]
            #ref = ref[1:]
            backg = backg[1:]
            # Crop around probe beam
            signal = signal[centre[0]-dy:centre[0]+dy,centre[1]-dx:centre[1]+dx]
            #ref = ref[centre[0]-dy:centre[0]+dy,centre[1]-dx:centre[1]+dx]
            backg = backg[centre[0]-dy:centre[0]+dy,centre[1]-dx:centre[1]+dx]
            # Sum and average
            res += (signal-backg) / cycles
            Data_ind = np.append(Data_ind,[signal-backg],axis=0)
            print "Signal image", j
        # Load final matrix and reset res
        Data = np.append(Data,[res],axis=0)
        res = np.zeros((2*dy, 2*dx))
        print "Phase", i/(j+1)+1 , "complete"
            
    # Clean the Data points matrix of zeros in the first slot 
    Data = np.delete(Data,0,axis=0)
    Data_ind = np.delete(Data_ind,0,axis=0)

    return Data, Data_ind
    
#    # Start main loop
#    for i in range(0,len(files),2):
#        pattern2 = files[i][len(prefix):len(prefix)+number_of_digits_for_files]
#        #print pattern2
#        # Load images for each iteration
#        signal = np.array(plt.imread(dname+files[i]),dtype=float64)
#        #ref = np.array(plt.imread(dname+files[i+1]),dtype=float64)
#        backg = np.array(plt.imread(dname+files[i+1]),dtype=float64)
#        # Correct for some weird first row jitter from Chameleon    
#        signal = signal[1:]
#        #ref = ref[1:]
#        backg = backg[1:]
#        # Crop around probe beam
#        signal = signal[centre[0]-dy:centre[0]+dy,centre[1]-dx:centre[1]+dx]
#    #    ref = ref[centre[0]-dy:centre[0]+dy,centre[1]-dx:centre[1]+dx]
#        backg = backg[centre[0]-dy:centre[0]+dy,centre[1]-dx:centre[1]+dx]
#        # Check data point    
#
#        if (i not in cond_list):
#            # Continue averaging if is same data point
#            res += signal-backg
#            count += 1
#            print "Signal image", i
#        else:  
#            # Finish averaging and pass to point slot
#            Data = np.append(Data,[res],axis=0)
#            # Adjust the pattern in filename
#            pattern = pattern2
#            # Reset container for averaging each data poitn, res
#            res = np.zeros((2*dy, 2*dx))
#            count = 0
#            print 'next data point: %s'%pattern
#            # Start new data point averaging
#            res += signal-backg
#            count += 1
#            print "Signal image", i
#    
#    # Clean the Data points matrix of zeros in the first slot 
#    # and append last point that skipped "else" condition 
#    Data = np.delete(Data,0,axis=0)
#    Data = np.append(Data,[res],axis=0)
#    for i in range(0,len(Data)):
#        Data[i] = Data[i]/n_avg
#    return Data


def tof(time,Data,Data_ind,atomic_mass):  
    plt.clf()
    tof_x = np.array([])
    tof_y = np.array([])
    cm_x = np.array([])
    cm_y = np.array([])
    yavg_data=np.array([])
    yavg=np.float
    
 #binning and plotting
    for j in range(0, len(Data)): 
        
        x=linspace(0,len(Data[j][1,:]),10000)
        y=linspace(0,len(Data[j][:,1]),10000)

        
        for i in range(0,len(Data[j][1,:])):
            yavg=Data[j][:,i]
            yavg_data=np.append(yavg_data,np.sum(yavg)/ len(Data[j][:,1])) #two different ways of summing were used just because
        
        fit_x=fitgaussian1d(None,yavg_data)
        
        xavg_data=np.zeros(len(Data[j][:,1]))
    
        for k in range(0,len(Data[j][:,1])):
            xavg_data[k]=np.sum(Data[j][k,:])/ len(Data[j][1,:]) 
        
        fit_y=fitgaussian1d(None,xavg_data)    
 
        plt.figure(1)        
        
        plt.plot(Data[j][fit_y[1],:],label='x')
        plt.plot(yavg_data,label='y_avg')
        plt.plot(x,gaussian1d(fit_x[0],fit_x[1],fit_x[2],fit_x[3])(x),label='xfit')
        plt.xlabel('x (pixel)')
        plt.ylabel('Fluorescence (a.u.)')
        
        plt.figure(2)

        plt.plot(Data[j][:,fit_x[1]],label='y')
        plt.plot(xavg_data,label='x_avg')
        plt.plot(y,gaussian1d(fit_y[0],fit_y[1],fit_y[2],fit_y[3])(y),label='yfit')
        plt.xlabel('z (pixel)')
        plt.ylabel('Fluorescence (a.u.)')
        
        tof_y = np.append(tof_y,np.abs(fit_y[2]))
        tof_x = np.append(tof_x,np.abs(fit_x[2]))
        cm_x = np.append(cm_x,np.abs(fit_x[1]))
        cm_y = np.append(cm_y,np.abs(fit_y[1]))
        yavg_data=np.array([])

        print "Data point", j       

    tofx_ind=np.array([])
    tofy_ind=np.array([])
    yind_data = np.array([])
    
    for j in range(0, len(Data_ind)): 
        
        for i in range(0,len(Data_ind[j][1,:])):
            yind=Data_ind[j][:,i]
            yind_data=np.append(yind_data,np.sum(yind)/ len(Data_ind[j][:,1])) #two different ways of summing were used just because
        
        fitx_ind=fitgaussian1d(None,yind_data)
        yind_data = np.array([])
        
        xind_data=np.zeros(len(Data_ind[j][:,1]))
    
        for k in range(0,len(Data_ind[j][:,1])):
            xind_data[k]=np.sum(Data_ind[j][k,:])/ len(Data_ind[j][1,:]) 
            
        fity_ind=fitgaussian1d(None,xind_data)  
        
        tofy_ind = np.append(tofy_ind,np.abs(fity_ind[2]))
        tofx_ind = np.append(tofx_ind,np.abs(fitx_ind[2]))

    tofx_error=np.array([])
    tofy_error=np.array([])
    tofx_test=np.array([])
    tofy_test=np.array([])

    for i in range(0,len(tofx_ind),cycles):
#        pattern2 = files[i][len(prefix):len(prefix)+number_of_digits_for_files]
        # print pattern2
        for j in range(0,cycles,1):
           tofx_test=np.append(tofx_test,tofx_ind[i+j])
           
        tofx_error = np.append(tofx_error,0.5*(max(tofx_test)-min(tofx_test)))
        tofx_test = np.array([])   
    
    for i in range(0,len(tofy_ind),cycles):
#        pattern2 = files[i][len(prefix):len(prefix)+number_of_digits_for_files]
        # print pattern2
        for j in range(0,cycles,1):
           tofy_test=np.append(tofy_test,tofy_ind[i+j])
           
        tofy_error = np.append(tofy_error,0.5*(max(tofy_test)-min(tofy_test)))
        tofy_test = np.array([])    

    
    plt.figure()
    times = np.array([])
    tof_ys = np.array([])
    tof_xs = np.array([])
    cm_ys = np.array([])
    cm_xs = np.array([])
    error_x = np.array([])
    error_y = np.array([])
    
    for i in time.argsort():
        tof_ys = np.append(tof_ys,tof_y[i]*pixel_size)
        tof_xs = np.append(tof_xs,tof_x[i]*pixel_size)
        cm_ys = np.append(cm_ys,cm_y[i]*pixel_size)
        cm_xs = np.append(cm_xs,cm_x[i]*pixel_size)
        error_x = np.append(error_x,tofx_error[i]*pixel_size)
        error_y = np.append(error_y,tofy_error[i]*pixel_size)
    times = np.sort(time)
    Temp_x, sigma0x = fit_tof(times,tof_xs)
    Temp_y, sigma0y = fit_tof(times,tof_ys)
    plt.scatter(times,tof_xs,label='x')
    plt.scatter(times,tof_ys,label='z', marker='o', color='g')
    plt.errorbar(times,tof_xs,yerr=error_x,xerr=None,ls='none', color='blue')
    plt.errorbar(times,tof_ys,yerr=error_y,xerr=None,ls='none', color='green')
    t = np.linspace(times[0],times[len(times)-1],len(times)*10)
    plt.plot(t,kinetic_expansion(Temp_x,sigma0x)(t))
    plt.plot(t,kinetic_expansion(Temp_y,sigma0y)(t))
    plt.xlabel('t (ms)')
    plt.ylabel('$\sigma$ ($\mu$m)')
    u=1.660538921e-27
    k=1.3806488e-23
    Temp_x = Temp_x * atomic_mass * u / k
    Temp_y = Temp_y * atomic_mass * u / k
    plt.text(times[1],tof_xs[len(tof_xs)-2],'Temp_x=%d uK'%Temp_x)    
    plt.text(times[1],tof_xs[len(tof_xs)-3],'Temp_z=%d uK'%Temp_y)
    plt.legend(loc=0)
    print "Temp_x = %d uK " %Temp_x
    print "Temp_z = %d uK"%Temp_y
    savepath = os.path.join(dname, 'tof.pdf') 
    plt.savefig(savepath)    
    return Temp_x, sigma0x, Temp_y, sigma0y, times, tof_xs, tof_ys, cm_xs, cm_ys
    
def atom_number(image): 
        
        power_pixel = np.float      
        N = np.float
        E = 2.54*10**(-16) #photon energy in milijoules
        eta = 1.19*10**(-10) #conversion of P from pixelheight to mW
        omega = 1.008*10**(-3) #solid angle
        
        p = fitgaussian2d(image) 
        centre = (p[1],p[2])
        sig_y = np.abs(p[3])
        sig_x = np.abs(p[4])
        #sum over 6*sigma to calculate the fluoresced power of the cloud in pixel units
        power_pixel=np.sum(image[centre[0]-3.5*sig_y:centre[0]+3.5*sig_y,centre[1]-3.5*sig_x:centre[1]+3.5*sig_x])
        #intensity of each beam is 3.8 mW/cm^2 for 205 mW MOT beam, detuning=-2.8G
        scatt_rate=(10**6)*scattering_rate(52.8,2.8)
        
        N = eta*power_pixel/(E*scatt_rate*omega)
        
        print "N = %.2e "%N
        print "sigma = %f mm"%(sig_x*0.0174)

 
        return N

data0, dataind0 = load_and_average_fluo(dname,raw_image,cycles)
Temp_x, sigma0x, Temp_y, sigma0y, times, tof_xs, tof_ys, cm_xs1, cm_ys = tof(time,data0,dataind0,atomic_mass)
atom_number(data0[0])

#t = np.linspace(times[0],times[len(times)-1],len(times)*10)
#A = np.vstack([times**2, np.ones(len(times))]).T
#a2,a0=np.linalg.lstsq(A,cm_xs1)[0]
#plot(times,cm_xs1,'o')    
#plot(t,a2*t**2+a0)
#plt.xlabel('t (ms)')
#plt.ylabel('x$_{cm}$ ($\mu$m)')