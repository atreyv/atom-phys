# -*- coding: utf-8 -*-
"""
Created on Mon Apr 7 10:02:56 2014

@author: pedro
"""
  
from libphys import *
import mayavi.mlab
import os


# Parameters to set for load_and_average
dname = 'D:/Data/2014/June/MOT/12_06_14/abs_test6/'
prefix = 'abs'
number_of_digits_for_files = 2
# Image to use for 2D gaussian fit. From this fit extract the gaussian zero
# and FWHM. FWHM is used as size of frame for operation: (sig-backg)/(ref-backg)
raw_image = 1


# Parameters for time of glight:
time = np.array([1,2,3,4,5,6,7,8,9,10,11])
pixel_size = 17.4 # microm
atomic_mass = 87

# Parameters for optical thickness:
# Create array with frequencies
# For absorption data 15 Apr2014:

freq=np.array([-2.19,0.03,1.95])

# Standard points
#freq = np.array([-3.04, -2.19, -1.26, -0.24, 
#                0.43,  0.91, 1.95, 2.99 ])

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
unsorted_resonance_point = 1
#
#
#
def load_and_average(dname,raw_image,prefix,number_of_digits_for_files):
    ## Load all file names in folder
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
    # Choose a size for the images good to perform division
    # Choose either the waist of the probe or the FWHM
    # The waist:
    #dx = int(gaussian_hwhm_to_radius(param[4]))
    #dy = int(gaussian_hwhm_to_radius(param[3]))
    # The FWHM:
    dx = int(param[4])
    dy = int(param[3])
    
    # Start the count, assign filename pattern and set the containers
    # for data: res for loading each iteration as transmission
    # and Data to load each data point.
    count = 0
    pattern = files[0][len(prefix):len(prefix)+number_of_digits_for_files]
    res = np.zeros((2*dy, 2*dx))
    Data = np.array([res],dtype=float64)
    
    # Start main loop
    for i in range(0,len(files),3):
        pattern2 = files[i][len(prefix):len(prefix)+number_of_digits_for_files]
    #    print pattern2
        # Load images for each iteration
        signal = np.array(plt.imread(dname+files[i]),dtype=float64)
        ref = np.array(plt.imread(dname+files[i+1]),dtype=float64)
        backg = np.array(plt.imread(dname+files[i+2]),dtype=float64)
        # Correct for some weird first row jitter from Chameleon    
        signal = signal[1:]
        ref = ref[1:]
        backg = backg[1:]
        # Crop around probe beam
        signal = signal[centre[0]-dy:centre[0]+dy,centre[1]-dx:centre[1]+dx]
        ref = ref[centre[0]-dy:centre[0]+dy,centre[1]-dx:centre[1]+dx]
        backg = backg[centre[0]-dy:centre[0]+dy,centre[1]-dx:centre[1]+dx]
        # Check data point    
        if (pattern == pattern2):
            # Continue averaging if is same data point
            res += normalize_by_division((signal-backg), (ref-backg))
            count += 1
            print "Signal image", i
        else:
            # Finish averaging and pass to point slot
            Data = np.append(Data,[res],axis=0)
            # Adjust the pattern in filename
            pattern = pattern2
            # Reset container for averaging each data poitn, res
            res = np.zeros((2*dy, 2*dx))
            count = 0
            print 'next data point: %s'%pattern
            # Start new data point averaging
            res += normalize_by_division((signal-backg), (ref-backg))
            count += 1
            print "Signal image", i
    
    # Clean the Data points matrix of zeros in the first slot 
    # and append last point that skipped "else" condition 
    Data = np.delete(Data,0,axis=0)
    Data = np.append(Data,[res],axis=0)
    for i in range(0,len(Data)):
        Data[i] = Data[i]/count
    return Data


def tof(time,Data,atomic_mass):  
    plt.clf()
    tof_x = np.array([])
    tof_y = np.array([])
    
    for i in range(0, len(Data)):
        p = fitgaussian2d(Data[i])
        tof_y = np.append(tof_y,gaussian_hwhm_to_stdev(np.abs(p[3])))
        tof_x = np.append(tof_x,gaussian_hwhm_to_stdev(np.abs(p[4])))
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
    plt.scatter(times,tof_ys,label='y', marker='x', color='g')
    x = np.linspace(times[0],times[len(times)-1],len(times)*10)
    plt.plot(x,kinetic_expansion(Tempx,sigma0x)(x))
    plt.plot(x,kinetic_expansion(Tempy,sigma0y)(x))
    plt.xlabel('t (ms)')
    plt.ylabel('$\sigma$ ($\mu$m)')
    Tempx = Tempx * atomic_mass * u / k
    Tempy = Tempy * atomic_mass * u / k
    plt.text(times[1],tof_xs[len(tof_xs)-2],'Temp_x=%d uK'%Tempx)    
    plt.text(times[1],tof_xs[len(tof_xs)-3],'Temp_y=%d uK'%Tempy)
    plt.legend(loc=0)
    print "Temp_x =", Tempx
    print "\nTemp_y =", Tempy
    
    return Tempx, sigma0x, Tempy, sigma0y, times, tof_xs, tof_ys, cm_xs, cm_ys


def optical_thickness(freq,Data,unsorted_resonance_point):
    # Create area to look at by using the on-resonance data point
    clf()
    p = fitgaussian2d(Data[unsorted_resonance_point])
    py = np.abs(p[3])*fraction
    px = np.abs(p[4])*fraction
    T = np.array([],dtype=float64)
    # Extract transmission through MOT
    for i in range(0, len(Data)):
#        p = fitgaussian2d(Data[i])
#        p[1] = np.abs(p[1])
#        p[2] = np.abs(p[2])
#        Data_roi = Data[i][p[1]-py:p[1]+py,p[2]-px:p[2]+px]
        Data_roi = Data[i][p[1]-py:p[1]+py,p[2]-px:p[2]+px]
        T = np.append(T,np.sum(Data_roi) / (4*px*py))
        print "Data point", i
#        plt.imshow(Data_roi)
#        savefig(str(i)+'jpg')
#        plt.colorbar()
    
    # Prepare arrays for sorted data
    freqs = np.array([])
    Ts = np.array([])
    for i in freq.argsort():
        Ts = np.append(Ts,T[i])
    freqs = np.sort(freq)
    
    # Perform fit do T versus detuning curve using T = exp(-b0/(1+4*nu**2/Gamma**2))
    b0, nu0 = fit_extinction_lorentz(freqs,Ts)
    fit = extinction_lorentz(b0,nu0)
    plt.text(0.,0.15,'b0 =%.1f'%b0)
    plt.scatter(freqs,Ts)
    x = np.linspace(freqs[0],freqs[len(freqs)-1],len(freqs)*10)
    plt.plot(x,fit(x))
    plt.xlabel('$\Delta$ ($\Gamma$)')
    plt.ylabel('T')
    print "b0 =", b0
    return b0, freqs, Ts



# Run functions
data = load_and_average(dname,raw_image,prefix,number_of_digits_for_files)
b0, freqs, Ts = optical_thickness(freq,data,unsorted_resonance_point)
