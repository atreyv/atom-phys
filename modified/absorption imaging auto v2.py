# -*- coding: utf-8 -*-
"""
Created on Mon Apr 7 10:02:56 2014

"""
  
from libphys import *
import mayavi.mlab
import os

font = {'family' : 'normal','weight' : 0,'size'  : 17}
matplotlib.rc('font', **font)
# Parameters to set for load_and_average
dname = 'D:/Data/2014/November/MOT/04_11_14/abs01/'
#dname = '/home/pedro/Downloads/LAB/DATA/2014/June/13_06_14/abs_test2/'

# Image to use for 2D gaussian fit. From this fit extract the gaussian zero
# and FWHM. FWHM is used as size of frame for operation: (sig-backg)/(ref-backg)
raw_image = 1
cycles = 5
phases = 6

# Parameters for time of glight:
time = np.array([2,3,4,5,6,7,8,9,10])
pixel_size = 17.4 # microm
atomic_mass = 87

# Parameters for optical thickness:
# Create array with frequencies

# Quick run
#freq=np.array([-2.19,0.03,1.95])

## Standard points

freq = np.array([0, -0.66, -1.29, -1.98, -2.9, -4.32])

#freq = np.array([-3.99,-2.19, -1.02, -0.24, 
#                -0.03,  0.43, 1.19])        
                #, 2.47, 3.92])
#freq = np.array([-3.04, -2.19, -1.26, -0.24, 0.43,  0.91, 1.95, 2.99])
# Full range                                            
#freq = np.array([-3.99,-3.04,-2.60,-2.19,-1.95,-1.68,-1.26,-1.02,
        #        -0.63,-0.24,0.03,0.43,0.66,0.91,1.19,1.42,1.71,
      #          1.95,2.47,2.99,3.92])

# For absorption data 15 Apr2014:
#freq = np.array([ -3.97, -2.74, -2.19, -1.65, -1.26, -0.87, -0.36, -0.24, 
             #    0.03, 0.28, 0.91, 1.43, 1.95, 2.47, 3.53])

# FWHM fraction of gaussian fit to a probe profile to use for integration
# of transmission intensity 
#fraction = .2
# On resonace unsorted data point (default = 5)
unsorted_resonance_point = 1
#
#
#
def load_and_average(dname,raw_image,cycles,fraction):
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
    dx = np.abs(int(param[4]*fraction))
    dy = np.abs(int(param[3]*fraction))
    
    # Set the containers
    # for data: res for loading each iteration as transmission
    # and Data to load each data point.
    res_avg = np.zeros((2*dy, 2*dx))
    Data_avg = np.array([res_avg],dtype=float64)
    Data_ind = np.array([res_avg],dtype=float64)
   
    # Start main loop
    for i in range(0,len(files),cycles*3):
#        pattern2 = files[i][len(prefix):len(prefix)+number_of_digits_for_files]
        # print pattern2
        for j in range(0,cycles*3,3):
            # Load images for each iteration
            signal = np.array(plt.imread(dname+files[j+i]),dtype=float64)
            ref = np.array(plt.imread(dname+files[j+i+1]),dtype=float64)
            backg = np.array(plt.imread(dname+files[j+i+2]),dtype=float64)
            # Correct for some weird first row jitter from Chameleon    
            signal = signal[1:]
            ref = ref[1:]
            backg = backg[1:]
            # Crop around probe beam
            signal = signal[centre[0]-dy:centre[0]+dy,centre[1]-dx:centre[1]+dx]
            ref = ref[centre[0]-dy:centre[0]+dy,centre[1]-dx:centre[1]+dx]
            backg = backg[centre[0]-dy:centre[0]+dy,centre[1]-dx:centre[1]+dx]
            # Sum and average
            res_avg += normalize_by_division((signal-backg), (ref-backg)) / cycles
            Data_ind = np.append(Data_ind,[normalize_by_division((signal-backg), (ref-backg))],axis=0)
            print "Signal image", j
        # Load final matrix and reset res
        Data_avg = np.append(Data_avg,[res_avg],axis=0)
        res_avg = np.zeros((2*dy, 2*dx))
        print "Phase", i/(j+1)+1 , "complete"
            
    # Clean the Data points matrix of zeros in the first slot 
    Data_avg = np.delete(Data_avg,0,axis=0)
    Data_ind = np.delete(Data_ind,0,axis=0)
    
    return Data_avg, Data_ind


def tof(time,Data,atomic_mass):  
    plt.clf()
    tof_x = np.array([])
    tof_y = np.array([])
    
    for i in range(0, len(Data)):
        p = fitgaussian2d(Data[i])
        tof_y = np.append(tof_y,np.abs(p[3]))
        tof_x = np.append(tof_x,np.abs(p[4]))
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


def optical_thickness(freq,Data_avg,Data_ind,size_point,fractionsum):
    # Create area to look at by using the on-resonance data point
    #p1 = fitgaussian2d(Data[size_point])
    
    yavg_data=np.array([])
    yavg=np.float
    
    for i in range(0,len(Data_avg[size_point][1,:])):
        yavg=-log(Data_avg[size_point][:,i])
        yavg_data=np.append(yavg_data,np.sum(yavg)/ len(Data_avg[size_point][:,1])) #two different ways of summing were used just because
        
    fit_x=fitgaussian1d(None,yavg_data)
    
    xavg_data=np.zeros(len(Data_avg[size_point][:,1]))
    
    for i in range(0,len(Data_avg[size_point][:,1])):
        xavg_data[i]=np.sum(-log(Data_avg[size_point][i,:]))/ len(Data_avg[size_point][1,:]) 
        
    fit_y=fitgaussian1d(None,xavg_data)
    
    center_x= fit_x[1]
    center_y= fit_y[1]
    px = int(np.abs(fractionsum*fit_x[2]))
    py = int(np.abs(fractionsum*fit_y[2]))
    sigx = int(np.abs(fit_x[2]))
    sigy = int(np.abs(fit_y[2]))
    T_avg = np.array([],dtype=float64)
    T_ind = np.array([],dtype=float64)  
    T_fixed_avg = np.array([],dtype=float64)  
    T_fixed_ind = np.array([],dtype=float64)  
    
    # Extract transmission through MOT
    for i in range(0, len(Data_avg)):
#        p = fitgaussian2d(Data[i])
#        p[1] = np.abs(p[1])
#        p[2] = np.abs(p[2])
#        Data_roi = Data[i][p[1]-py:p[1]+py,p[2]-px:p[2]+px]
        Data_roi_avg = Data_avg[i][center_y-py:center_y+py,center_x-px:center_x+px]
        T_avg = np.append(T_avg,np.sum(Data_roi_avg) / (4*px*py))
        print "Data_avg point", i
#        plt.imshow(Data_roi)
#        savefig(str(i)+'jpg')
#        plt.colorbar()

    for i in range(0, len(Data_avg)):
#calc an average transmission for integration along 0.4x0.4 mm^2 (size of pump beam)    
        Data_roi_fix = Data_avg[i][center_y-23:center_y+23,center_x-23:center_x+23]
        T_fixed_avg = np.append(T_fixed_avg,np.sum(Data_roi_fix) / (4*23*23))
        print "Data_avg_fixed point", i
        
    for i in range(0, len(Data_ind)):
#        p = fitgaussian2d(Data[i])
#        p[1] = np.abs(p[1])
#        p[2] = np.abs(p[2])
#        Data_roi = Data[i][p[1]-py:p[1]+py,p[2]-px:p[2]+px]
        Data_roi_ind = Data_ind[i][center_y-py:center_y+py,center_x-px:center_x+px]
        T_ind = np.append(T_ind,np.sum(Data_roi_ind) / (4*px*py))
        print "Data_ind point", i

    for i in range(0, len(Data_ind)):
#        p = fitgaussian2d(Data[i])
#        p[1] = np.abs(p[1])
#        p[2] = np.abs(p[2])
#        Data_roi = Data[i][p[1]-py:p[1]+py,p[2]-px:p[2]+px]
        Data_roi_fixed = Data_ind[i][center_y-23:center_y+23,center_x-23:center_x+23]
        T_fixed_ind = np.append(T_fixed_ind,np.sum(Data_roi_fixed) / (4*23*23))
        print "Data__ind_fixed point", i
        
    # Prepare arrays for sorted data
    freqs = np.array([])
    Ts_avg = np.array([])
    for i in freq.argsort():
        Ts_avg = np.append(Ts_avg,T_avg[i])
    freqs = np.sort(freq)

    # Prepare arrays for sorted data
    freqs = np.array([])
    Ts_fixed_avg = np.array([])
    for i in freq.argsort():
        Ts_fixed_avg = np.append(Ts_fixed_avg,T_fixed_avg[i])
    freqs = np.sort(freq)
    
    T_error = np.array([])
    T_test = np.array([])
    T_fixed_error = np.array([])
    T_fixed_test = np.array([])

    for i in range(0,len(T_ind),cycles):
#        pattern2 = files[i][len(prefix):len(prefix)+number_of_digits_for_files]
        # print pattern2
        for j in range(0,cycles,1):
           T_test=np.append(T_test,T_ind[i+j])
           
        T_error = np.append(T_error,0.5*(max(T_test)-min(T_test)))
        T_test = np.array([])

    for i in range(0,len(T_fixed_ind),cycles):
#        pattern2 = files[i][len(prefix):len(prefix)+number_of_digits_for_files]
        # print pattern2
        for j in range(0,cycles,1):
           T_fixed_test=np.append(T_fixed_test,T_fixed_ind[i+j])
           
        T_fixed_error = np.append(T_fixed_error,0.5*(max(T_fixed_test)-min(T_fixed_test)))
        T_fixed_test = np.array([])        
    
    # Perform fit do T versus detuning curve using T = exp(-b0/(1+4*nu**2/Gammasize_point**2))
    figure()
    b0, nu0 = fit_extinction_lorentz(freqs,Ts_avg)
    print "b0 = ", b0
    sigma_z = sigy*pixel_size
    sigma_x = sigx*pixel_size
    print "sigma_x = %.f um"%sigma_x
    print "sigma_z = %.f um"%sigma_z
    lambda_ = 0.78024
    sig_res = 7./15. * (3*lambda_**2) / (2.*np.pi)
    N = 2*np.pi * sigma_z * sigma_x * b0 / sig_res
    n0=b0/(sig_res*1e-08)
    print "N = %.1e"%N
    print "n0 = %.1e 1/cm^2"%n0
    fit = extinction_lorentz(b0,nu0)
    plt.text(-1.65,0.4,'b0 =%.1f'%b0)
    plt.text(-1.65,0.48,'n$_0$ =%.1e cm$^{-2}$'%n0)
    plt.text(-1.65,0.56,'N =%.1e'%N)
    plt.scatter(freqs,Ts_avg)
    plt.errorbar(freqs,Ts_avg,yerr=T_error,xerr=None,ls='none', color='blue')
    x = np.linspace(freqs[0],freqs[len(freqs)-1],len(freqs)*10)
    plt.plot(x,fit(x),)
    plt.xlabel('$\Delta$ ($\Gamma$)')
    plt.ylabel('T')
    size_point
     # The final path to save to
    savepath = os.path.join(dname, 'b0.pdf') 
    plt.savefig(savepath)
    
    # Perform fit do T versus detuning curve using T = exp(-b0/(1+4*nu**2/Gammasize_point**2)) for a fixed integration of cloud
    figure()
    b0_fixed, nu0_fixed = fit_extinction_lorentz(freqs,Ts_fixed_avg)
    print "b0 fixed = ", b0_fixed
    sigma_z = sigy*pixel_size
    sigma_x = sigx*pixel_size
    print "sigma_x = %.f um"%sigma_x
    print "sigma_z = %.f um"%sigma_z
    lambda_ = 0.78024
    sig_res = 7./15. * (3*lambda_**2) / (2.*np.pi)
    N_fixed = 2*np.pi * sigma_z * sigma_x * b0_fixed / sig_res
    n0_fixed=b0_fixed/(sig_res*1e-08)
    print "N = %.1e"%N_fixed
    print "n0 = %.1e 1/cm^2"%n0_fixed
    fit = extinction_lorentz(b0_fixed,nu0_fixed)
    plt.text(-1.65,0.4,'b0 =%.1f'%b0_fixed)
    plt.text(-1.65,0.48,'n$_0$ =%.1e cm$^{-2}$'%n0_fixed)
    plt.text(-1.65,0.56,'N =%.1e'%N_fixed)
    plt.scatter(freqs,Ts_fixed_avg)
    plt.errorbar(freqs,Ts_fixed_avg,yerr=T_fixed_error,xerr=None,ls='none', color='blue')
    x = np.linspace(freqs[0],freqs[len(freqs)-1],len(freqs)*10)
    plt.plot(x,fit(x),)
    plt.xlabel('$\Delta$ ($\Gamma$)')
    plt.ylabel('T')
    size_point
     # The final path to save to
    savepath = os.path.join(dname, 'b0_fixed.pdf') 
    plt.savefig(savepath)

    
    
    return b0, N, n0, Ts_avg

def atom_number_abs(Data,size_point,fractionfit,fractionsum):

    yavg_data=np.array([])
    yavg=np.float
    x=linspace(0,len(Data[size_point][1,:]),10000)
    
    for i in range(0,len(Data[size_point][1,:])):
        yavg=-log(Data[size_point][:,i])
        yavg_data=np.append(yavg_data,np.sum(yavg)/ len(Data[size_point][:,1])) #two different ways of summing were used just because
        
    fit_x=fitgaussian1d(None,yavg_data)
    
    xavg_data=np.zeros(len(Data[size_point][:,1]))
    y=linspace(0,len(Data[size_point][:,1]),10000)
    
    for i in range(0,len(Data[size_point][:,1])):
        xavg_data[i]=np.sum(-log(Data[size_point][i,:]))/ len(Data[size_point][1,:]) 
    

    
    fit_y=fitgaussian1d(None,xavg_data)
    
    center_x= fit_x[1]
    center_y= fit_y[1]
    peakb_x = int(np.abs(fractionfit*fit_x[2]))
    peakb_y = int(np.abs(fractionfit*fit_y[2]))
    peakb = np.sum(-log(Data[size_point][center_y-peakb_y:center_y+peakb_y,center_x-peakb_x:center_x+peakb_x]))/(4*peakb_x*peakb_y)
    
    N1_calc=2*np.pi*np.abs(fit_x[2]*fit_y[2])*peakb*pixel_size**2
        
    figure()
    plot(-log(Data[size_point][center_y,:]),label='x')
    plot(yavg_data,label='y_avg')
    plot(x,gaussian1d(fit_x[0],fit_x[1],fit_x[2],fit_x[3])(x),label='xfit')
    plot(-log(Data[size_point][:,center_x]),label='y')
    plot(xavg_data,label='x_avg')
    plot(y,gaussian1d(fit_y[0],fit_y[1],fit_y[2],fit_y[3])(y),label='yfit')

    #p1 = fitgaussian2d(-log(Data[size_point])) 

    px= int(fractionsum*np.abs(fit_x[2]))
    py= int(fractionsum*np.abs(fit_y[2]))

    N2_calc=np.sum(-log(Data[size_point][center_y-py:center_y+py,center_x-px:center_x+px]))*pixel_size**2

    lambda_ = 0.78024
    sig_res = 7./15.*(3*lambda_**2)/(2.*np.pi)
    sig_det = sig_res/(1+4*freq[size_point]**2)
    sigma_x= np.abs(fit_x[2])*pixel_size
    sigma_z= np.abs(fit_y[2])*pixel_size  
    N_factor=1 / sig_det
    
    N1=N_factor*N1_calc
    N2=N_factor*N2_calc
    
    print "N1 = %.1e"%N1
    print "N2 = %.1e"%N2
    print "sigma_x = %.1e um"%sigma_x
    print "sigma_z = %.1e um"%sigma_z
 
    return N1, N2
    
## Run functions for OD
size_point=3

#data_avg, data_ind = load_and_average(dname,raw_image,cycles,0.9)
N1, N2 = atom_number_abs(data_avg,size_point,0.5,2.7)
b, N, n0, Ts = optical_thickness(freq,data_avg,data_ind,size_point,0.5)

#b1, freqs, Ts = optical_thickness(freq,data,unsorted_resonance_point,1.25)

#N3, N4 = atom_number_abs(data,size_point,0.3,2.5)
#N5, N6 = atom_number_abs(data,size_point,0.4,2.5)
#N7, N8 = atom_number_abs(data,size_point,0.5,2.5)
#N9, N10 = atom_number_abs(data,size_point,0.6,2.5)
#N11, N12 = atom_number_abs(data,size_point,0.7,2.5)
#numbers = N1, N3, N5, N7, N9, N11
#fractions = 0.2, 0.3, 0.4, 0.5, 0.6, 0.7
#
#plot(fractions,numbers)


# Runa functions for ToF
#prefix = 'tof'
#data = load_and_average(dname,raw_image,cycles)
#Tempx, sigma0x, Tempy, sigma0y, times, tof_xs, tof_ys, cm_xs, cm_ys = tof(time,data,atomic_mass)
