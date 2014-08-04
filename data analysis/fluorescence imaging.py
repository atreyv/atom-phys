# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 14:05:59 2014

@author: lab
"""

from libphys import *
#import mayavi.mlab
import os


# Parameters to set for load_and_average
#dname = 'D:/Data/2014/June/MOT/16_06_14/tof04/'
dname = '/home/pedro/Downloads/LAB/DATA/2014/June/23_06_14/tof04/'
prefix = 'fluo'
number_of_digits_for_files = 1
cycles = 5 # number of averages
phases = 9 # number of interesting cases
# Image to use for 2D gaussian fit. From this fit extract the gaussian zero
# and FWHM. FWHM is used as size of frame for operation: (sig-backg)/(ref-backg)
raw_image = 120

# Parameters for time of glight:
time = np.array([2,3,4,5,6,7,8,9,10,11,12,13,14])
pixel_size = 17.4 # microm
atomic_mass = 87

# Parameters for optical thickness:
# Create array with frequencies

# Quick run
#freq=np.array([-2.19,0.03,1.95])

# Standard points
freq = np.array([-3.04, -2.19, -1.26, -0.24, 
                0.43,  0.91, 1.95, 2.99 ])

# For absorption data 15 Apr2014:
#freq = np.array([ -3.97, -2.74, -2.19, -1.65, -1.26, -0.87, -0.36, -0.24, 
             #    0.03, 0.28, 0.91, 1.43, 1.95, 2.47, 3.53])
                                            
#freq = np.array([0.099,0.198,0.330,0.462,0.561,0.758,0.857,0.989,1.121,\
#1.352,-0.132,-1.220,-1.319,-1.484,1.616,1.846,2.110,2.407,-1.714,-2.011,\
#-2.275,-0.297,2.671,2.901,3.165,3.693,-2.572,-2.901,-3.132,-3.808,-0.396,\
#-0.528,-0.659,-0.791,-0.923,-1.088,0.000])
# FWHM fraction of gaussian fit to a probe profile to use for integration
# of transmission intensity 
fraction = .2
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
    dx = int(param[4]) *2
    dy = int(param[3]) *2
    
    # Start the count, assign filename pattern and set the containers
    # for data: res for loading each iteration as transmission
    # and Data to load each data point.
   # count = 0
   # pattern = files[0][len(prefix):len(prefix)+number_of_digits_for_files]
    print "dy",dy
    print "dx",dx
    res = np.zeros((2*dy, 2*dx))
    Data = np.array([res],dtype=float64)

    
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
            print "Signal image", j
        # Load final matrix and reset res
        Data = np.append(Data,[res],axis=0)
        res = np.zeros((2*dy, 2*dx))
        print "Phase", i/(j+1)+1 , "complete"
            
    # Clean the Data points matrix of zeros in the first slot 
    Data = np.delete(Data,0,axis=0)
    return Data    
    
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


def tof(time,Data,atomic_mass):  
    plt.clf()
    tof_x = np.array([])
    tof_y = np.array([])
    
    for i in range(0, len(Data)):
        p = fitgaussian2d(Data[i])
        tof_y = np.append(tof_y,np.abs(p[3]))#gaussian_hwhm_to_stdev(np.abs(p[3])))
        tof_x = np.append(tof_x,np.abs(p[4]))#gaussian_hwhm_to_stdev(np.abs(p[4])))
        cm_x = np.append(tof_y,np.abs(p[2]))
        cm_y = np.append(tof_x,np.abs(p[1]))
        print "Data point", i
    
    times = np.array([])
    tof_ys = np.array([])
    tof_xs = np.array([])
    cm_ys = np.array([])
    cm_xs = np.array([])
    for i in time.argsort():
        tof_ys = np.append(tof_ys,tof_y[i]*pixel_size)
        tof_xs = np.append(tof_xs,tof_x[i]*pixel_size)
        cm_ys = np.append(cm_ys,cm_y[i]*pixel_size)
        cm_xs = np.append(cm_xs,cm_x[i]*pixel_size)
    times = np.sort(time)
    Tempx, sigma0x = fit_tof(times,tof_xs)
    Tempy, sigma0y = fit_tof(times,tof_ys)
    plt.scatter(times,tof_xs,label='x')
    plt.scatter(times,tof_ys,label='z', marker='x', color='g')
    x = np.linspace(times[0],times[len(times)-1],len(times)*10)
    plt.plot(x,kinetic_expansion(Tempx,sigma0x)(x))
    plt.plot(x,kinetic_expansion(Tempy,sigma0y)(x))
    plt.xlabel('t (ms)')
    plt.ylabel('$\sigma$ ($\mu$m)')
    Tempx = Tempx * atomic_mass * u / k
    Tempy = Tempy * atomic_mass * u / k
    plt.text(times[1],tof_xs[len(tof_xs)-2],'Temp_x=%d uK'%Tempx)    
    plt.text(times[1],tof_xs[len(tof_xs)-3],'Temp_z=%d uK'%Tempy)
    plt.legend(loc=0)
    print "Temp_x =", Tempx
    print "\nTemp_z =", Tempy
    savepath = os.path.join(dname, 'tof.pdf') 
    plt.savefig(savepath)    
    return Tempx, sigma0x, Tempy, sigma0y, times, tof_xs, tof_ys, cm_xs, cm_ys


data = load_and_average_fluo(dname,raw_image,cycles)
Tempx, sigma0x, Tempy, sigma0y, times, tof_xs, tof_ys, cm_xs, cm_ys = tof(time,data,atomic_mass)