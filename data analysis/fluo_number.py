# -*- coding: utf-8 -*-
"""
Created on Wed Jul 02 19:00:36 2014

@author: lab
"""

from libphys import *
#import mayavi.mlab
import os


# Parameters to set for load_and_average
dname = 'Z:/Data/2014/July/MOT/31_07_14/fluo1/'
#prefix = 'tof'
#number_of_digits_for_files = 1
averages = 5 # number of averages
#phases = 9 # number of interesting cases
# Image to use for 2D gaussian fit. From this fit extract the gaussian zero
# and FWHM. FWHM is used as size of frame for operation: (sig-backg)/(ref-backg)

# Parameters for time of glight:
#time = np.array([6,7,8,9,10,11,12,13,14])
#time = np.array([2,3,4,5,6,7,8,9,10])
pixel_size = .00174 # milimeter
atomic_mass = 87

# Parameters for optical thickness:
# Create array with frequencies

# Quick run
#freq=np.array([-2.19,0.03,1.95])

#loading times
loading_times = np.array([100,300,500,700,900,1200,2000,3000,4000,5000,6000,7000,8000])
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
def load_and_average_num(dname,averages):
    ## Load all file names in folder
    files = []
    for file in os.listdir(dname):
        if file.endswith(".bmp"):
            files = np.append(files,file)
    files.sort()
    print 'Found %d files' %len(files)
    raw_image = len(files)-2
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
    dx = int(param[4]*3.5) 
    dy = int(param[3]*3.5)
    
    # Start the count, assign filename pattern and set the containers
    # for data: res for loading each iteration as transmission
    # and Data to load each data point.
   # count = 0
   # pattern = files[0][len(prefix):len(prefix)+number_of_digits_for_files]
    print "dy",dy
    print "dx",dx
    res = np.zeros((2*dy,2*dx))
    Data = np.array([res],dtype=float64)    

    for i in range(0,len(files),2*averages):
#        pattern2 = files[i][len(prefix):len(prefix)+number_of
        for j in range(0,2*averages,2):
        #for j in range(0,averages*2,2):
            # Load images for each iteration
            signal = np.array(plt.imread(dname+files[i+j]),dtype=float64)
            #ref = np.array(plt.imread(dname+files[j+i+1]),dtype=float64)
            backg = np.array(plt.imread(dname+files[i+j+1]),dtype=float64)
            # Correct for some weird first row jitter from Chameleon    
            signal = signal[1:]
            #ref = ref[1:]
            backg = backg[1:]
            # Crop around atom cloud to select only the data you need
            signal = signal[centre[0]-dy:centre[0]+dy,centre[1]-dx:centre[1]+dx]
            #ref = ref[centre[0]-dy:centre[0]+dy,centre[1]-dx:centre[1]+dx]
            backg = backg[centre[0]-dy:centre[0]+dy,centre[1]-dx:centre[1]+dx]
            # Sum and average
            res += (signal-backg)/averages
            print "Signal image", j
        # Load final matrix and reset res
        Data = np.append(Data,[res],axis=0)
        res = np.zeros((2*dy, 2*dx))
        print "Phase", i/(j+1)+1 , "complete"

    # Clean the Data points matrix of zeros in the first slot (due to jitter)
    Data = np.delete(Data,0,axis=0)
    
    return Data   
    
    
    
def atom_number(image): 
        
        power_pixel = np.float      
        N = np.float
        E = 2.54*10**(-16) #photon energy in milijoules
        eta = 1.19*10**(-10) #conversion of P from pixelheight to mW
        omega = 1.008*10**(-3) #solid angle
        
        p = fitgaussian2d(image) 
        centre = (p[1],p[2])
        sig_y = np.abs(p[3])
        sig_x = np.abs(p[K. Lindquist, M. Stephens, and C. Wieman, Phys. Rev.
A 46, 4082 (1992).4])
        #sum over 6*sigma to calculate the fluoresced power of the cloud in pixel units
        power_pixel=np.sum(image[centre[0]-3.5*sig_y:centre[0]+3.5*sig_y,centre[1]-3.5*sig_x:centre[1]+3.5*sig_x])
        #intensity of each beam is 3.8 mW/cm^2 for 205 mW MOT beam, detuning=-2.8G
        scatt_rate=(10**6)*scattering_rate(52.8,2.8)
        
        N = eta*power_pixel/(E*scatt_rate*omega)
        
        print "N = %.2e "%N
        print "sigma = %f "%(sig_x*0.0174)

 
        return N
    
def number_time(Data,loading_times):
    numbers = np.zeros(len(Data))
    
    #load numbers for all of the points    
    for i in range(0,len(Data),1):
            numbers[i] = atom_number(Data[i])

    #plot all the points
    plt.scatter(loading_times,numbers/10**8,label='N')
    plt.xlabel('t (ms)')
    plt.ylabel('N x 10$^8$')
    plt.xlim(-250,8500)
    #x_pow=linspace(power_pixel[0],power_pixel[len(power)-1],num=1000)
    #plt.plot(x_pow,powerfit[1]+powerfit[0]*x_pow)
    #plt.text(power_pixel[len(power)-12],power[len(power)-2],'eta_pow=%.2e mW/pixelheight (200us, gain=0)'%powerfit[0])    
    #plt.legend(loc=0)
    savepath = os.path.join(dname, 'N_fluo.pdf')  
    plt.savefig(savepath) 

    return numbers


fluodata=load_and_average_num(dname,averages)
atomnum=atom_number(fluodata[0])