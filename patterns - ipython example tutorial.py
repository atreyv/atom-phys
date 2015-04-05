
# coding: utf-8

# ####Import all necessary libraries. Libpatterns imports libphys as well

# In[ ]:

get_ipython().magic(u'matplotlib inline')
from libpatternsworkflow import *
from __future__ import division
mpl.rc('font', size=14)
mpl.rcParams['figure.figsize'] = (16.0, 8.0)


# ####Set the directories where the orthogonal and/or parallel polarisation images are. Also set the folder where the logged intensity is

# In[ ]:

dname_ortho = '/home/pedro/LAB/DATA/2014/Dec/PF/12_12_14/pf02/cam448/'
dname_parallel = '/home/pedro/LAB/DATA/2014/Dec/PF/12_12_14/pf02/cam113/'
iname = '/home/pedro/LAB/DATA/2014/Dec/PF/12_12_14/pf02/'
files_ortho = load_files(dname_ortho,".tif")
files_parallel = load_files(dname_parallel,".tif")
files_int = load_files_prefix(iname, prefix='probe', sufix='.dat')


# ####Set the parameters with the information about which images contain patterns (fstart and fstop) and the position in filename to extract the relevant x-axis variable (start and stop)

# In[ ]:

fstart = 0
fstop = 280
start = 5
stop = 8
(files_parallel[280][start:stop])


# ###Set some of experimental parameters.
# 
# * frac is the fraction of gaussian sigma got from 2D fit to ref, to set the size in FFT

# In[ ]:

pixel_size = 3.75 # micrometers
magnification = 2.
raw_image = 280
# Set the fraction of the sigma from the gaussian fit to a reference. This will set the integrable area for
# energy measurements.
frac = 1.
#half_size = 35
symmetry_order = 6
azimuthal_profile_smoothness = 20
scaling_angle_ortho = 12
scaling_angle_parallel = -12
compression_y_over_x = 1./0.84
# After rotation and compression, some cropping is needed to regain the square boundaries for the near field image
image_crop_factor = 0.6
# The pixel bin radius from which the polynomial fit should be applied
start_pos_for_noise_correction_in_fspace = 30
energy_ratio_condition = 0.0


# ###Use a reference image to locate the centre of pump/probe for ortho/parallel pol

# For the orthogonal polarisation; notice that a pure reference is a flat plane intensity profile for linear polarisation at the Brewster angle, so a typical pattern is selected

# In[ ]:

fft_size_ortho, param_ortho, refft_ortho = locate_roi_from_ref(dname_ortho, files_ortho, scaling_angle_ortho,
                                                               0, frac, compression_y_over_x,
                                                               image_crop_factor)


# For the parallel polarisation; a pure reference should be used

# In[ ]:

fft_size_parallel, param_parallel, refft_parallel = locate_roi_from_ref(dname_parallel, files_parallel,
                                                                        scaling_angle_parallel, raw_image, frac,
                                                                        compression_y_over_x, image_crop_factor)


# One single fft_size should used for proper transmission calculations. The smallest size between parallel and orthogonal is typical chosen, although is as good metrics and any other.

# In[ ]:

param_ortho[3:5] = param_parallel[3:5]
fft_size = fft_size_parallel

#param_parallel[3:5] = param_ortho[3:5]
#fft_size = fft_size_ortho
fft_size,param_ortho,param_parallel


# Now pure references should be measured for proper total power sent through the atomic cloud and for correct transmission calculations. Typically this will involve recalculating the fft of a pure reference from the orthogonal pol when the pump pol is linear.

# In[ ]:

def calibrate_intensity(i):
    ref_ortho = read_file_to_ndarray(dname_ortho+files_ortho[i])
    ref_ortho = scale_image(ref_ortho,scaling_angle_ortho,compression_y_over_x,interpolation_value=0)
    ref_ortho = image_crop(ref_ortho,image_crop_factor)
    ref_ortho = prepare_for_fft_crop(ref_ortho,fft_size)
    ref_parallel = read_file_to_ndarray(dname_parallel+files_parallel[i])
    ref_parallel = scale_image(ref_parallel,scaling_angle_parallel,compression_y_over_x,interpolation_value=0)
    ref_parallel = image_crop(ref_parallel,image_crop_factor)
    ref_parallel = prepare_for_fft_crop(ref_parallel,fft_size)
    P0 = ref_parallel.sum() + ref_ortho.sum()
    return P0

P0_ccd = np.array([calibrate_intensity(raw_image)])
for i in xrange(raw_image+1,raw_image+20):
    P0_ccd = np.append(P0_ccd, calibrate_intensity(i))
    
I0_probe_pd = get_probe_intensity_profile_from_txt(iname+'ref.dat',0.08,intensity_plateau_n_points=300,
                                                   smoothness_points=1200,plot_all=False,check_plot=False,
                                                   plot_find_peaks=False)

I0_probe_cal = np.average(I0_probe_pd / P0_ccd)
I0_probe = np.average(P0_ccd * I0_probe_cal)

I0_pump_pd = get_pump_intensity_profile_from_txt(iname+files_int[0],0.08,300,1200,True,True,True)
for i in files_int[1:]:
    I0_pump_pd = np.append(I0_pump_pd,get_pump_intensity_profile_from_txt(iname+i,0.08,300,1200,True,True,True))

#P0, fft_size, refft_ortho.shape, refft_parallel.shape
#ref_ortho.shape, refft_parallel.shape, ref_ortho[0,0]

#%time circle_line_integration(refft_ortho,0)[0]


# In[ ]:

plt.plot(I0_pump_pd)
np.average(I0_pump_pd), np.std(I0_pump_pd)


# In[ ]:

plt.figure(figsize=(12,8))
plt.plot((P0_ccd-np.amin(P0_ccd))/(np.amax(P0_ccd)-np.amin(P0_ccd))*100, marker='o', linewidth=1)
plt.plot((I0_probe_pd-np.amin(I0_probe_pd))/(np.amax(I0_probe_pd)-np.amin(I0_probe_pd))*100,
         marker='x', linewidth=1)
plt.xlabel('# acquisition')
plt.ylabel('% $\Delta$')
#plt.savefig(iname+'intensity_green-image_power_blue-relation.pdf')


# ###Set the main function to all data metrics

# Write the main function for looping. It might be used with joblib parallelisation. In this case the plots don't show in the end, although the call for plot wastes processor time.

# In[ ]:

def get_data_metrics(i, energy_ratio_condition, param1, scaling_angle1, dname1, files1, I_cal=0, plots=False,
                     peak_plot=False, azimuthal_plots=False, voronoi_plots=False,
                     param2=0, scaling_angle2=0, dname2=0, files2=0):
    
    # Read the pattern image from the filename and folder
    image = read_file_to_ndarray(dname1+files1[i])
    # Use the parameters from reference to select relevant area and correct the astigmatism that might be
    # present in optical system
    image = prepare_for_fft_full_image(image,param1,frac)
    image = scale_image(image,scaling_angle1,compression_y_over_x,interpolation_value=0)
    image1 = image_crop(image,image_crop_factor)
    
    # Get the Fourier Trasnf image
    resft1 = do_fft(image1)
    
    # If there is a second image (polarisation):
    if scaling_angle2 != 0:
        image = read_file_to_ndarray(dname2+files2[i])
        image = prepare_for_fft_full_image(image,param2,frac)
        image = scale_image(image,scaling_angle2,compression_y_over_x,interpolation_value=0)
        image2 = image_crop(image,image_crop_factor)
        resft2 = do_fft(image2)
    
    if plots:
        # Set the subplots for visual confirmation of peaks, etc
        fig = plt.figure()
        plot1 = fig.add_subplot(121)
        plot2 = fig.add_subplot(122)
        imshowfft(plot1,resft1,0.2,True)
        plot1.text(5,-5,'t1=%d,i=%d'%(int(files1[i][start:stop]),i),fontsize=16,color='black')

        if scaling_angle2 != 0:
            fig2 = plt.figure()
            plot3 = fig2.add_subplot(121)
            plot4 = fig2.add_subplot(122)
            imshowfft(plot3,resft2,0.2,True)
            plot3.text(5,-5,'t2=%d,i=%d'%(int(files2[i][start:stop]),i),fontsize=16,color='black')
    
    # Set the container that will pick up the azimuthal integral for each radius. Starts at radius = 1,
    # so careful calling the array (that will start at zero).
    radial_plot1 = np.array([])
    for j in range(0,int(len(resft1)/2)):
        radial_plot1 = np.append(radial_plot1,circle_line_integration(resft1,j)[0])
    # Remove the integrated noise. Select the initial readius where one shouldn't get anymore peaks (start_pos)
    correction = gets_integration_noise_on_fourier_space(radial_plot1,
                                                         start_pos=start_pos_for_noise_correction_in_fspace)
    radial_plot1 -= correction(np.arange(0,len(radial_plot1)))
#    if plots:
#        plot2.plot(radial_plot1)

    if scaling_angle2 != 0:
        radial_plot2 = np.array([])
        for j in range(0,int(len(resft2)/2)):
            radial_plot2 = np.append(radial_plot2,circle_line_integration(resft2,j)[0])
        correction = gets_integration_noise_on_fourier_space(radial_plot2,
                                                             start_pos=start_pos_for_noise_correction_in_fspace)
        radial_plot2 -= correction(np.arange(0,len(radial_plot2)))
#        if plots:
#            plot4.plot(radial_plot2)

    # Get the peaks from radial_plot. The method doesn't always work well, so play with ratio of
    # smoothness/interpolation_points. If none is found, nothing will be returned by the main function!
    peaks1_temp = find_peaks(radial_plot1,interpolation_points=1000,peak_finding_smoothness=5,
                             plot=peak_plot, plot_new_fig=True)
    if (peaks1_temp != 0):# and peaks_temp.y_raw[0] > 0):
        
        # Get the the first ring in the radial coordinate
        if plots:
            p1 = fit_ft_peak(wavevector_order=1, radial_spread=1, radial_plot=radial_plot1[:],
                             peaks_temp=peaks1_temp, fit='no_offset', plots=True, subplot=plot2)
        p1 = fit_ft_peak(wavevector_order=1, radial_spread=1, radial_plot=radial_plot1[:],
                         peaks_temp=peaks1_temp, fit='no_offset', plots=False)
        
        if scaling_angle2 != 0:
            peaks2_temp = find_peaks(radial_plot2,interpolation_points=1000,peak_finding_smoothness=5,
                                     plot=peak_plot, plot_new_fig=True)
            if plots:
                    p2 = fit_ft_peak(wavevector_order=1, radial_spread=2, radial_plot=radial_plot2[:],
                                     peaks_temp=peaks2_temp, fit='offset', plots=True, subplot=plot4)
            p2 = fit_ft_peak(wavevector_order=1, radial_spread=2, radial_plot=radial_plot2[:],
                             peaks_temp=peaks2_temp, fit='offset', plots=False)
            
        E1 = energy_ratio_wavevector_ring(1.,p=p1)
        if scaling_angle2 != 0:
            E2 = energy_ratio_wavevector_ring(1.,p=p2)
        
        # The first condition will be the first order sideband (considering there is one) having enough weigth,
        # and the peak detection method not mess around the profile during normalisation.
        if (E1 > energy_ratio_condition):
            
            # Collect the real position of the peak
            Lambda1 = 1. / (p1[1] / (pixel_size/magnification*fft_size))
            # Collect the total intensity transmitted through the atoms
            trans_pow1 = radial_plot1[0] * I_cal
            # Collect all from the ring
            ring1_area = E1
            ring1_amplitude = p1[0]
            ring1_width = 2 * p1[2]
            symmetry_order_measured1,            azimuthal_peak1 = get_data_azimuthal_metrics(resft1, p1[1], which_sideband=0,
                                                         radial_epsilon=2,
                                                         interpolated_points=1000,
                                                         azimuthal_profile_smoothness=20,
                                                         plots=azimuthal_plots)
            count1, position_x1, position_y1, side_ave1,            side_std1 = get_data_voronoi_metrics(image1, symmetry_order_measured1, gaussian_sigma=10,
                                                 pixel_size=pixel_size, magnification=magnification,
                                                 plots=voronoi_plots)
            if scaling_angle2 != 0:
                Lambda2 = 1. / (p2[1] / (pixel_size/magnification*fft_size))
                trans_pow2 = radial_plot2[0] * I_cal
                ring2_area = E2
                ring2_amplitude = p2[0]
                ring2_width = 2 * p2[2]
                #azimuthal_peaks2 = get_data_azimuthal_metrics(resft2, p2[1])

            # Collect the x-axis coordinate
            t = float(files1[i][start:stop])
            if plots:
                plt.show()
            if scaling_angle2 != 0:
                return np.array([t, Lambda1, trans_pow1, ring1_area, ring1_amplitude, ring1_width,
                                Lambda2, trans_pow2, ring2_area, ring2_amplitude, ring2_width,
                                symmetry_order_measured1, azimuthal_peak1])
            else:
                return np.array([t, Lambda1, trans_pow1, ring1_area, ring1_amplitude, ring1_width,
                                symmetry_order_measured1, azimuthal_peak1])
        elif scaling_angle2 != 0:
            return np.ones(13)*(-1)
        else:
            return np.ones(8)*(-1)
    elif scaling_angle2 != 0:
        return np.ones(13)*(-1)
    else:
        return np.ones(8)*(-1)


# def main():
#     data = [get_data_metrics(fstart, energy_ratio_condition, param_ortho, scaling_angle_ortho,
#                              dname_ortho, files_ortho, I0_probe_cal, True, False, False, True,
#                              param_parallel, scaling_angle_parallel, dname_parallel, files_parallel)]
#     for i in xrange(fstart+1,fstop):
#         data = np.append(data, [get_data_metrics(i, energy_ratio_condition, param_ortho,
#                                                  scaling_angle_ortho, dname_ortho, files_ortho, I0_probe_cal,
#                                                  True, False, False, True, param_parallel, scaling_angle_parallel,
#                                                  dname_parallel, files_parallel)],axis=0)
#     return data
# %time data = main()

# In[ ]:

get_ipython().magic(u'time data = np.asarray(Parallel(n_jobs=4)(delayed(get_data_metrics)(i, energy_ratio_condition, param_ortho,                                                                     scaling_angle_ortho, dname_ortho,                                                                     files_ortho,                                                                     I0_probe_cal, False, False, False, False,                                                                     param_parallel,                                                                     scaling_angle_parallel, dname_parallel,                                                                     files_parallel)                                           for i in range(fstart,fstop)))')


# ###Start the main for-loop that will collect all data from each image and append it to the created container
# 
# The returned tuple is composed in order by:
# 
# 0. * t represents the x-axis coordinate relevant from the present analysed data: might be time, mirror distance, etc.
# 1. * Lambda1 is the polarisation1 lengthscale.
# 2. * trans_pow1 is its transmitted power, should be divided by total power
# 3. * ring1_area is its selected peak order are contained in 2*sigma
# 4. * ring1_amplitude is its slected peak amplitude
# 5. * ring1_width is the 2*sigma width
# 6. * Lambda2 is the polarisation2 lengthscale.
# 7. * trans_pow2 is its transmitted power, should be divided by total power
# 8. * ring2_area
# 9. * ring2_amplitude
# 10. * ring2_width
# 11. * intensities is the container for logged intensity at the photodiode.
# 12. * voronoi_maps contains the relevant data for checking translational symmetry breaking, distinguish positive from negative hexagons

# In[ ]:

def truth_intensities(imin,imax):
    truth_int1 = I0_pump_pd[:] > imin
    truth_int2 = I0_pump_pd[:] < imax
    truth_int = truth_int1 * truth_int2
    return truth_int

def average_std_data(x_var, data_var, accept_interval_in_std, truth_int):
    value1_mean = np.array([])
    value1_errors = np.array([])
    t_value = np.array([])
    values_final_count = np.array([])
    j = np.inf
    for i in x_var:
        if i != j and i != -1:
            truth = x_var==i
            t_value = np.append(t_value,i)
            value1_std = np.std(data_var[truth*truth_int])
            value1_average = np.average(data_var[truth*truth_int])
            value1_truth1 = data_var < value1_average + accept_interval_in_std * value1_std
            value1_truth2 = data_var > value1_average - accept_interval_in_std * value1_std
            value1_truth = value1_truth1 * value1_truth2
            value1_mean = np.append(value1_mean, np.average(data_var[truth*truth_int*value1_truth]))
            value1_errors = np.append(value1_errors, np.std(data_var[truth*truth_int*value1_truth]))
            j = i
            values_final_count = np.append(values_final_count,(truth_int * truth * value1_truth).sum())
        else:
            pass
    return t_value, value1_mean, value1_errors, values_final_count

accept_interval_in_std = 2
truth_int = truth_intensities(17,2100)

t_value, value1_mean_ortho, value1_errors_ortho,count1_ortho = average_std_data(data[:,0], data[:,2], accept_interval_in_std, truth_int)

t_value, value2_mean_ortho, value2_errors_ortho,count2_ortho = average_std_data(data[:,0], data[:,4], accept_interval_in_std, truth_int)

t_value, value1_mean_parallel, value1_errors_parallel,count1_parallel = average_std_data(data[:,0], data[:,7]+data[:,2], accept_interval_in_std, truth_int)

t_value, value2_mean_parallel, value2_errors_parallel,count2_parallel = average_std_data(data[:,0], data[:,9], accept_interval_in_std, truth_int)

value2_mean_ortho = value2_mean_ortho / (value1_mean_ortho / I0_probe_cal)
value2_errors_ortho = value2_errors_ortho / (value1_mean_ortho / I0_probe_cal)
value2_mean_parallel = value2_mean_parallel / (value1_mean_parallel / I0_probe_cal)
value2_errors_parallel = value2_errors_parallel / (value1_mean_parallel / I0_probe_cal)


# In[ ]:

fig = plt.figure(figsize=(8,6))
plot1 = fig.add_subplot(111)
plot1.plot(t_value,value1_mean_ortho,linewidth=0,marker='o')
plot1.errorbar(t_value,value1_mean_ortho,yerr=value1_errors_ortho,ls='none')
plot1.set_xlabel(r't ($\mathrm{\mu s}$)')
#plot1.set_ylabel(r'$\Lambda$ ($\mathrm{\mu m}$)')
plot1.set_ylabel(r'$I_{trans}$ ($mW/cm^2$)')
plot1.set_xlim(xmin=0,xmax=160)
plot1.set_ylim(ymin=0)


# In[ ]:

fig = plt.figure(figsize=(8,6))
plot2 = fig.add_subplot(111)
plot2.plot(t_value,value2_mean_ortho,linewidth=0,marker='o')
plot2.errorbar(t_value,value2_mean_ortho,yerr=value2_errors_ortho,ls='none')
plot2.set_xlabel(r't ($\mathrm{\mu s}$)')
plot2.set_ylabel(r'$E_{q_1}$/$E_{0}$')
plot2.set_xlim(xmin=0,xmax=160)
plot2.set_ylim(ymin=0)


# In[ ]:

fig = plt.figure(figsize=(8,6))
plot3 = fig.add_subplot(111)
plot3.plot(t_value,value1_mean_parallel,linewidth=0,marker='o')
plot3.errorbar(t_value,value1_mean_parallel,yerr=value1_errors_parallel,ls='none')
plot3.set_xlabel(r't ($\mathrm{\mu s}$)')
#plot3.set_ylabel(r'$\Lambda$ ($\mathrm{\mu m}$)')
plot3.set_ylabel(r'$I_{trans}$ ($mW/cm^2$)')
plot3.set_xlim(xmin=0,xmax=160)
plot3.set_ylim(ymin=0)


# In[ ]:

fig = plt.figure(figsize=(8,6))
plot4 = fig.add_subplot(111)
plot4.plot(t_value,value2_mean_parallel,linewidth=0,marker='o')
plot4.errorbar(t_value,value2_mean_parallel,yerr=value2_errors_parallel,ls='none')
#plot4.plot(t_value,intense/np.average(intense))
plot4.set_xlabel(r't ($\mathrm{\mu s}$)')
plot4.set_ylabel(r'$E_{q_1}$/$E_{0}$')
plot4.set_xlim(xmin=0,xmax=160)
plot4.set_ylim(ymin=0)
#fig.savefig(dname+'ortho-postmirror-pol.pdf')


# In[ ]:

count1_ortho, count2_ortho, count1_parallel, count2_parallel


# In[ ]:

for i in xrange(0,20):
    resft1 = refft_ortho
    # Read the pattern image from the filename and folder
    image = read_file_to_ndarray(dname_ortho+files_ortho[i])
    # Use the parameters from reference to select relevant area and correct the astigmatism that might be
    # present in optical system
    image = prepare_for_fft_full_image(image,param_ortho,frac)
    image = scale_image(image,scaling_angle_ortho,compression_y_over_x,interpolation_value=0)
    image1 = image_crop(image,image_crop_factor)
    try:
        data_voronoi = np.append(data_voronoi, [get_data_voronoi_metrics(image1,6,10,True)],axis=0)
    except:
        data_voronoi = [get_data_voronoi_metrics(image1,6,10,True)]


# In[ ]:

data_voronoi[:,1:3]


# In[ ]:



