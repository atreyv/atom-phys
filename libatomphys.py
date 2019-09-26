# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 15:47:18 2015

@author: pedro
"""

from libphys import *
import matplotlib


def double_gaussian1d(A1, x1, sig1, A0):
    sig1 = float(sig1)
    return lambda x: (A0*np.exp(-(x-x0)**2/(2*sig0**2)) + offset) * (1 - A1*np.exp(-(x-x1)**2/(2*sig1**2)))

def double_gauss1d_moments(x,data):
    total = data.sum()
    if (x==0):
        x = np.arange(data.size)
    x1 = (x*data).sum() / total
    sig1 = np.sqrt(np.sum((x-x1)**2*data)/np.sum(data))
    A0 = np.amax(data)
    A1 = -A0
    return A1, x1, sig1, A0
    
def fit_double_gauss1d(x,data):
    params = moments1d(x,data)
    if (x==0):
        errorfunction = lambda p: double_gaussian1d(*p)(*np.indices(data.shape))\
                                - data
    else:
        errorfunction = lambda p: double_gaussian1d(*p)(x)\
                                - data
    p, success = spleastsq(errorfunction, params, full_output=0)
    return p


def load_and_average_fluo_no_backg(dname,raw_image,number_of_points_per_cycle,
                                    file_ext=".bmp"):
    ## Load all file names in folder
    files = []
    for file in os.listdir(dname):
        if file.endswith(file_ext):
            files = np.append(files,file)
    files.sort()
    print ('Found %d files' %len(files))
    
    # Do a fit to probe using a reference image. Also correct for jitter in
    # first row from Chameleon
    ref = np.array(plt.imread(dname+files[raw_image]),dtype=np.float64)
    ref = ref[1:]

    #centre = np.unravel_index(ref.argmax(), ref.shape)
    param = fitgaussian2d(ref)
    centre = (np.int(param[1]),np.int(param[2]))

    # Choose a size for the images good to perform division
    # Choose either the waist of the probe or the FWHM
    # The FWHM:
    dx = np.int(param[4]) *2
    dy = np.int(param[3]) *2
    
    # Start the count, assign filename pattern and set the containers
    # for data: res for loading each iteration as transmission
    # and Data to load each data point.
    print ("dy",dy)
    print ("dx",dx)
    res = np.zeros((2*dy, 2*dx))
    Data = np.array([res],dtype=np.float64)
    Data_ind = np.array([res],dtype=np.float64)

    # Start main loop
    for i in range(0,len(files),number_of_points_per_cycle):
        # print pattern2
        for j in range(0,number_of_points_per_cycle,1):
            # Load images for each iteration
            signal = np.array(plt.imread(dname+files[j+i]),dtype=np.float64)
            # Correct for some weird first row jitter from Chameleon    
            signal = signal[1:]
            # Crop around probe beam
            signal = signal[centre[0]-dy:centre[0]+dy,centre[1]-dx:centre[1]+dx]
            res += signal / number_of_points_per_cycle
            Data_ind = np.append(Data_ind,[signal],axis=0)
            print ("Signal image", j)
        # Load final matrix and reset res
        Data = np.append(Data,[res],axis=0)
        res = np.zeros((2*dy, 2*dx))
        print ("Phase", i , "complete")
            
    # Clean the Data points matrix of zeros in the first slot 
    Data = np.delete(Data,0,axis=0)
    Data_ind = np.delete(Data_ind,0,axis=0)
    return Data, Data_ind


def load_and_average_abs_single_ref_no_backg(dname,dname_ref,cycles,fraction,
                                            file_ext=".bmp"):
    ## Load all file names in folders
    ref = np.array([])
    count = 0
    for file in os.listdir(dname_ref):
        if file.endswith(file_ext):
            files_ref = np.array(plt.imread(dname_ref+file),dtype=np.float64)
            files_ref = files_ref[1:]
            try:
                ref += files_ref
                count += 1.
            except:
                ref = np.zeros(files_ref.shape)
                ref += files_ref
                count += 1.
    ref = ref / count
    #print count

    files = []
    for file in os.listdir(dname):
        if file.endswith(file_ext):
            files = np.append(files,file)
    files.sort()
    print ('Found %d files' %len(files))

    # Choose a size for the images good to perform division
    paramy = fitgaussian1d(np.nan,np.average(ref,axis=1))
    paramx = fitgaussian1d(np.nan,np.average(ref,axis=0))
#    print paramx
#    param = fitgaussian2d(ref)
    centre = (np.int(np.round(paramy[1])),np.int(np.round(paramx[1])))
    dx = np.abs(np.int(np.round(paramx[2]*fraction)))
    dy = np.abs(np.int(np.round(paramy[2]*fraction)))
#    ref_adjust = np.sum(ref[:centre[0]-dy,:]) + np.sum(ref[centre[0]+dy:,:]) +\
#                np.sum(ref[centre[0]-dy:centre[0]+dy,:centre[1]-dx]) + \
#                np.sum(ref[centre[0]-dy:centre[0]+dy,centre[1]+dx:])
    ref = ref[np.int(centre[0]-dy):np.int(centre[0]+dy),np.int(centre[1]-dx):np.int(centre[1]+dx)]
     
    # Set the containers
    # for data: res for loading each iteration as transmission
    # and Data to load each data point.
    res_avg = np.zeros((np.int(2*dy), np.int(2*dx)))
    Data_avg = np.array([res_avg],dtype=np.float64)
    Data_ind = np.array([res_avg],dtype=np.float64)
   
    # Start main loop
    for i in range(0,len(files),cycles):
        for j in range(0,cycles):
            # Load images for each iteration
            signal = np.array(plt.imread(dname+files[j+i]),dtype=np.float64)
            # Correct for some weird first row jitter from Chameleon    
            signal = signal[1:]
            # create ref from signal
            signal_x_ref = np.average(signal,axis=0)
            signal_y_ref = np.average(signal,axis=1)
            x_ref = np.append(np.arange(0,np.int(centre[1]-dx)),np.arange(np.int(centre[1]+dx),len(signal_x_ref)))
            y_ref = np.append(np.arange(0,np.int(centre[0]-dy)),np.arange(np.int(centre[0]+dy),len(signal_y_ref)))
            signal_x_ref = np.append(signal_x_ref[:np.int(centre[1]-dx)],signal_x_ref[np.int(centre[1]+dx):])
            signal_y_ref = np.append(signal_y_ref[:np.int(centre[0]-dy)],signal_y_ref[np.int(centre[0]+dy):])
            
            Ax1,x1,sig_x1,offset_x1 = fitgaussian1d(x_ref,signal_x_ref)
            Ay1,y1,sig_y1,offset_y1 = fitgaussian1d(y_ref,signal_y_ref)
#            A1, y1, x1, sigy1, sigx1, offset = fitgaussian2d([y_ref,x_ref],signal[y_ref,x_ref])
            # Crop around probe beam
            signal = signal[np.int(centre[0]-dy):np.int(centre[0]+dy),np.int(centre[1]-dx):np.int(centre[1]+dx)]
            
            # Adjust reference to signal intensity using the gaussian tails
            # from the probe
#            adjust = ref_adjust / (np.sum(signal[:centre[0]-dy,:]) +\
#                                np.sum(signal[centre[0]+dy:,:]) +\
#                                np.sum(signal[centre[0]-dy:centre[0]+dy,:centre[1]-dx]) +\
#                                np.sum(signal[centre[0]-dy:centre[0]+dy,centre[1]+dx:]))
#            print adjust
            
#            res = normalize_by_division(signal,ref)
#            global x0
#            x0 = paramx2[1]
#            global sig0
#            sig0 = dx2
#            global offset
#            offset = paramx2[3]
#            A1, x1, sig1, A0 = fit_double_gauss1d(0,np.average(res,axis=0))
            adjust_y = (Ay1) / paramy[0]
            adjust_x = (Ax1) / paramx[0]
            print (adjust_x)
#            adjust = 1
            # Sum and average
            res = normalize_by_division(signal, (ref-paramx[3])*adjust_x+offset_x1)
            res_avg += res / cycles
            Data_ind = np.append(Data_ind,[res],axis=0)
#            res = signal.sum(axis=1) / (gaussian1d(A0,x0,sig0,paramy[3])(np.arange(len(signal.sum(axis=1)))))
#            res_avg += res / (cycles)
#            Data_ind = np.append(Data_ind,res)
            print ("Signal image", j)
        # Load final matrix and reset res
        Data_avg = np.append(Data_avg,[res_avg],axis=0)
        res_avg = np.zeros((2*dy, 2*dx))
        print ("Phase", i/(cycles)+1 , "complete")
    
    # Clean the Data points matrix of zeros in the first slot 
    Data_avg = np.delete(Data_avg,0,axis=0)
    Data_ind = np.delete(Data_ind,0,axis=0)
    
    return Data_avg, Data_ind, ref


def optical_thickness(dname,freq,Data,trans_fraction,pixel_size,unsorted_cloud_size_point):
    # Create area to look at by using the Data point more suitable for size
    p_binning_y = fitgaussian1d(np.nan,Data[unsorted_cloud_size_point].sum(axis=1))
    p_binning_x = fitgaussian1d(np.nan,Data[unsorted_cloud_size_point].sum(axis=0))
    #p = fitgaussian2d(Data[unsorted_cloud_size_point])
    py = np.abs(p_binning_y[2])*trans_fraction
    p1 = p_binning_y[1]
    #py = np.abs(p[3])*trans_fraction
    px = np.abs(p_binning_x[2])*trans_fraction    
    p2 = p_binning_x[1]
    #px = np.abs(p[4])*trans_fraction
    T = np.array([],dtype=np.float64)
    # Extract transmission through MOT
    for i in range(0, len(Data)):
#        p = fitgaussian2d(Data[i])
#        p[1] = np.abs(p[1])
#        p[2] = np.abs(p[2])
#        Data_roi = Data[i][p[1]-py:p[1]+py,p[2]-px:p[2]+px]
        Data_roi = Data[i][np.int(p1-py):np.int(p1+py),np.int(p2-px):np.int(p2+px)]
        T = np.append(T,np.sum(Data_roi) / (4*px*py))
        print ("Data point", i)
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
    print ("b0 = %.1f, delta = %.1f $\Gamma$"%(b0,nu0))
    sigma_y = np.abs(p_binning_y[2])*pixel_size
    sigma_x = np.abs(p_binning_x[2])*pixel_size
    print ("sigma_x = %.f um"%sigma_x)
    print ("sigma_y = %.f um"%sigma_y)
    lambda_ = 0.78024
    sig_res = 7./15. * (3*lambda_**2) / (2.*np.pi)
    N = 2*np.pi * sigma_y * sigma_x * b0 / sig_res
    print ("N = %.1e"%N)
    fit = extinction_lorentz(b0,nu0)
    plt.text(-0.5,0.35,'b0 =%.1f'%b0)
    plt.scatter(freqs,Ts)
    x = np.linspace(freqs[0],freqs[len(freqs)-1],len(freqs)*10)
    plt.plot(x,fit(x))
    plt.xlabel('$\Delta$ ($\Gamma$)')
    plt.ylabel('T')
    
     # The final path to save to
    savepath = os.path.join(dname, 'b0.pdf') 
    plt.savefig(savepath)
    return b0, freqs, Ts



def optical_thickness2(dname,freq,Data_avg,Data_ind,cycles,pixel_size,size_point,fractionsum):
    # Create area to look at by using the on-resonance data point
    #p1 = fitgaussian2d(Data[size_point])
    
    yavg_data=np.array([])
    yavg=np.float
    
    for i in range(0,len(Data_avg[size_point][1,:])):
        yavg=-np.log(Data_avg[size_point][:,i])
        yavg_data=np.append(yavg_data,np.sum(yavg)/ len(Data_avg[size_point][:,1])) #two different ways of summing were used just because
        
    fit_x=fitgaussian1d(np.nan,yavg_data)
    
    xavg_data=np.zeros(len(Data_avg[size_point][:,1]))
    
    for i in range(0,len(Data_avg[size_point][:,1])):
        xavg_data[i]=np.sum(-np.log(Data_avg[size_point][i,:]))/ len(Data_avg[size_point][1,:]) 
        
    fit_y=fitgaussian1d(np.nan,xavg_data)
    
    center_x = np.int(fit_x[1])
    center_y = np.int(fit_y[1])
    px = np.int(np.abs(fractionsum*fit_x[2]))
    py = np.int(np.abs(fractionsum*fit_y[2]))
    sigx = np.int(np.abs(fit_x[2]))
    sigy = np.int(np.abs(fit_y[2]))
    T_avg = np.array([],dtype=np.float64)
    T_ind = np.array([],dtype=np.float64)  
    T_fixed_avg = np.array([],dtype=np.float64)  
    T_fixed_ind = np.array([],dtype=np.float64)  
    
    # Extract transmission through MOT
    for i in range(0, len(Data_avg)):
#        p = fitgaussian2d(Data[i])
#        p[1] = np.abs(p[1])
#        p[2] = np.abs(p[2])
#        Data_roi = Data[i][p[1]-py:p[1]+py,p[2]-px:p[2]+px]
        Data_roi_avg = Data_avg[i][np.int(center_y-py):np.int(center_y+py),np.int(center_x-px):np.int(center_x+px)]
        T_avg = np.append(T_avg,np.sum(Data_roi_avg) / (4*px*py))
        print ("Data_avg point", i)
#        plt.imshow(Data_roi)
#        savefig(str(i)+'jpg')
#        plt.colorbar()

    for i in range(0, len(Data_avg)):
#calc an average transmission for integration along 0.4x0.4 mm^2 (size of pump beam)    
        Data_roi_fix = Data_avg[i][np.int(center_y-23):np.int(center_y+23),np.int(center_x-23):np.int(center_x+23)]
        T_fixed_avg = np.append(T_fixed_avg,np.sum(Data_roi_fix) / (4*23*23))
        print ("Data_avg_fixed point", i)
        
    for i in range(0, len(Data_ind)):
#        p = fitgaussian2d(Data[i])
#        p[1] = np.abs(p[1])
#        p[2] = np.abs(p[2])
#        Data_roi = Data[i][p[1]-py:p[1]+py,p[2]-px:p[2]+px]
        Data_roi_ind = Data_ind[i][center_y-py:center_y+py,center_x-px:center_x+px]
        T_ind = np.append(T_ind,np.sum(Data_roi_ind) / (4*px*py))
        print ("Data_ind point", i)

    for i in range(0, len(Data_ind)):
#        p = fitgaussian2d(Data[i])
#        p[1] = np.abs(p[1])
#        p[2] = np.abs(p[2])
#        Data_roi = Data[i][p[1]-py:p[1]+py,p[2]-px:p[2]+px]
        Data_roi_fixed = Data_ind[i][center_y-23:center_y+23,center_x-23:center_x+23]
        T_fixed_ind = np.append(T_fixed_ind,np.sum(Data_roi_fixed) / (4*23*23))
        print ("Data__ind_fixed point", i)
        
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
    plt.figure()
    b0, nu0 = fit_extinction_lorentz(freqs,Ts_avg)
    print ("b0 = ", b0)
    sigma_z = sigy*pixel_size
    sigma_x = sigx*pixel_size
    print ("sigma_x = %.f um"%sigma_x)
    print ("sigma_z = %.f um"%sigma_z)
    lambda_ = 0.78024
    sig_res = 7./15. * (3*lambda_**2) / (2.*np.pi)
    N = 2*np.pi * sigma_z * sigma_x * b0 / sig_res
    n0=b0/(sig_res*1e-08)
    print ("N = %.1e"%N)
    print ("n0 = %.1e 1/cm^2"%n0)
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
    #size_point
     # The final path to save to
    savepath = os.path.join(dname, 'b0.pdf') 
    plt.savefig(savepath)
    
    # Perform fit do T versus detuning curve using T = exp(-b0/(1+4*nu**2/Gammasize_point**2)) for a fixed integration of cloud
    plt.figure()
    b0_fixed, nu0_fixed = fit_extinction_lorentz(freqs,Ts_fixed_avg)
    print ("b0 fixed = ", b0_fixed)
    sigma_z = sigy*pixel_size
    sigma_x = sigx*pixel_size
    print ("sigma_x = %.f um"%sigma_x)
    print ("sigma_z = %.f um"%sigma_z)
    lambda_ = 0.78024
    sig_res = 7./15. * (3*lambda_**2) / (2.*np.pi)
    N_fixed = 2*np.pi * sigma_z * sigma_x * b0_fixed / sig_res
    n0_fixed=b0_fixed/(sig_res*1e-08)
    print ("N = %.1e"%N_fixed)
    print ("n0 = %.1e 1/cm^2"%n0_fixed)
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
    #size_point
     # The final path to save to
    savepath = os.path.join(dname, 'b0_fixed.pdf') 
    plt.savefig(savepath)

    return b0, N, n0, Ts_avg





def tof2(time,Data,Data_ind,atomic_mass,cycles,pixel_size,dname):  
    plt.clf()
    tof_x = np.array([])
    tof_y = np.array([])
    cm_x = np.array([])
    cm_y = np.array([])
    yavg_data=np.array([])
    yavg=np.float
    
 #binning and plotting
    for j in range(0, len(Data)): 
        
        x=np.linspace(0,len(Data[j][1,:]),10000)
        y=np.linspace(0,len(Data[j][:,1]),10000)

        
        for i in range(0,len(Data[j][1,:])):
            yavg=Data[j][:,i]
            yavg_data=np.append(yavg_data,np.sum(yavg)/ len(Data[j][:,1])) #two different ways of summing were used just because
        
        fit_x=np.round(fitgaussian1d(np.nan,yavg_data)).astype(int)
        
        xavg_data=np.zeros(len(Data[j][:,1]))
    
        for k in range(0,len(Data[j][:,1])):
            xavg_data[k]=np.sum(Data[j][k,:])/ len(Data[j][1,:]) 
        
        fit_y=np.round(fitgaussian1d(np.nan,xavg_data)).astype(int)
 
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

        print ("Data point", j)

    tofx_ind=np.array([])
    tofy_ind=np.array([])
    yind_data = np.array([])
    
    for j in range(0, len(Data_ind)): 
        
        for i in range(0,len(Data_ind[j][1,:])):
            yind=Data_ind[j][:,i]
            yind_data=np.append(yind_data,np.sum(yind)/ len(Data_ind[j][:,1])) #two different ways of summing were used just because
        
        fitx_ind=fitgaussian1d(np.nan,yind_data)
        yind_data = np.array([])
        
        xind_data=np.zeros(len(Data_ind[j][:,1]))
    
        for k in range(0,len(Data_ind[j][:,1])):
            xind_data[k]=np.sum(Data_ind[j][k,:])/ len(Data_ind[j][1,:]) 
            
        fity_ind=fitgaussian1d(np.nan,xind_data)  
        
        tofy_ind = np.append(tofy_ind,np.abs(fity_ind[2]))
        tofx_ind = np.append(tofx_ind,np.abs(fitx_ind[2]))

    tofx_error=np.array([])
    tofy_error=np.array([])
    tofx_test=np.array([])
    tofy_test=np.array([])

    for i in range(0,len(tofx_ind),cycles):
        # print pattern2
        for j in range(0,cycles,1):
           tofx_test=np.append(tofx_test,tofx_ind[i+j])
           
        tofx_error = np.append(tofx_error,0.5*(max(tofx_test)-min(tofx_test)))
        tofx_test = np.array([])   
    
    for i in range(0,len(tofy_ind),cycles):
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
    print ("Temp_x = %d uK " %Temp_x)
    print ("Temp_z = %d uK"%Temp_y)
    savepath = os.path.join(dname, 'tof.pdf') 
    plt.savefig(savepath)    
    return Temp_x, sigma0x, Temp_y, sigma0y, times, tof_xs, tof_ys, cm_xs, cm_ys
