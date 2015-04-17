# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 16:50:55 2015

@author: pedro

This library is not self-suficient. It requires initialisation
of global variables for proper call of methods. 

"""

#from libphys import *
from libpatterns import *
from joblib import Parallel, delayed

def locate_roi_from_ref(dname, files, scaling_angle, raw_image, frac,
                        compression_y_over_x, image_crop_factor, plots=True):
    image = read_file_to_ndarray(dname+files[raw_image])
    param = use_ref_to_locate_centre(image)
    # Crop the image for FFT using the gaussian fit; the returned image is a square.
    ref = prepare_for_fft_full_image(image,param,frac)
    # Correct the astigmatism from the optical system
    ref = scale_image(ref,scaling_angle,compression_y_over_x,interpolation_value=0)
    ref = image_crop(ref,image_crop_factor)
    fft_size = np.shape(ref)[0]
    refft = do_fft(ref)
    if plots:
        fig = plt.figure()
        plot1 = fig.add_subplot(121)
        plot2 = fig.add_subplot(122)
        temp = imshowfft(plot1,ref,1,logscale=False)
        temp = imshowfft(plot2,refft,.2,logscale=True)
    return fft_size, param, refft

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

def get_data_azimuthal_metrics(resft1, pos, which_sideband, radial_epsilon, interpolated_points=1000,
                               azimuthal_profile_smoothness=20, plots=False):
    azimuth_interp,\
    azimuthal = get_azimuthal_profile_from_ft_integrated_along_radius(resft1,int(np.round(pos)),
                                                                      radial_epsilon,interpolated_points)
    peaks_azimuthal = find_peaks(azimuthal,interpolated_points,azimuthal_profile_smoothness,plots)
    if (len(peaks_azimuthal.peaks['peaks'][0]) % 2 != 0):
        azimuthal = np.append(azimuthal,azimuthal[:int(30/360.*interpolated_points)])
        peaks_azimuthal = find_peaks(azimuthal,interpolated_points*390/360,azimuthal_profile_smoothness,plots)
    symmetry_order_measured = len(peaks_azimuthal.peaks['peaks'][0])
    sorted2 = np.sort(peaks_azimuthal.peaks['peaks'][0])
    if sorted2[symmetry_order_measured-1] >= interpolated_points:
        sorted2[symmetry_order_measured-1] -= interpolated_points
        sorted2 = np.sort(sorted2)
    azimuthal_peaks = azimuth_interp[sorted2.astype(np.int32)]
    return  np.array([symmetry_order_measured, azimuthal_peaks[which_sideband]])

def get_data_voronoi_metrics(image1, symmetry_order_measured, gaussian_sigma=10,
                             pixel_size=3.75, magnification=2, 
                             plots=False):
    resft_smoothed = gaussian_filter(image1, sigma=gaussian_sigma)
    detected = detect_peaks(resft_smoothed)
    vor = np.where(detected)
    vor = np.transpose(np.array([vor[1],vor[0]]))
    vor = Voronoi(vor)
    if plots:
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')
        ax.imshow(resft_smoothed,origin='upper')
        voronoi = voronoi_plot_2d(vor,ax)
        plt.axis('off')
        ax.invert_yaxis()
    count = 0
    for i in vor.regions:
        if len(i)==symmetry_order_measured and i.count(-1)==0:
            count += 1.
    side = np.array([])
    for j in vor.regions:
        if len(j)==symmetry_order_measured and j.count(-1)==0:
            for i in xrange(0,len(j)):
                side = np.append(side,np.linalg.norm(vor.vertices[j[i]]-vor.vertices[j[i-1]]))
    side = side * pixel_size / magnification
    position = np.array([np.average(vor.points[:,0]), np.average(vor.points[:,1])]) * pixel_size / magnification
    return np.array([count, position[0], position[1], np.average(side), np.std(side)])

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

def truth_intensities(I0_pump_pd, imin, imax):
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




# This the main function. Should be always last, except if some other
# function depends on it.

def get_data_metrics(i,
                     start,
                     stop,
                     pixel_size,
                     magnification,
                     raw_image,
                     frac,
                     symmetry_order,
                     azimuthal_profile_smoothness,
                     compression_y_over_x,
                     image_crop_factor,                    
                     start_pos_for_noise_corr_in_fspace,
                     energy_ratio_condition,
                     intensity_plateau_n_points,
                     int_smoothness,
                     fft_size,
                     param1, scaling_angle1, dname1,
                     files1, I_cal=0, plots=False,
                     peak_plot=False, azimuthal_plots=False,
                     voronoi_plots=False,
                     param2=0, scaling_angle2=0, dname2=0, files2=0):
    """Main function for pattern metrics.
    The returned tuple is composed in order by:

0. * t represents the x-axis coordinate relevant from the present analysed data: might be time, mirror distance, etc.
1. * Lambda1 is the polarisation1 lengthscale.
2. * trans_pow1 is its transmitted power, should be divided by total power
3. * ring1_area is its selected peak order are contained in 2*sigma
4. * ring1_amplitude is its slected peak amplitude
5. * ring1_width is the 2*sigma width
6. * Lambda2 is the polarisation2 lengthscale.
7. * trans_pow2 is its transmitted power, should be divided by total power
8. * ring2_area
9. * ring2_amplitude
10. * ring2_width
11. * symmetry_order_measured1 is the counted number of sideband peaks
12. * azimuthal_peak1 is the azimuthal angle of the first peak from quadrant 1
13. * count1 is the number of voronoi regions with the counted symmetry order
14. * position_x1 is the c.o.m. in the x-coordinate of the voronoi point
15. * position_y1 is the c.o.m. in the y-coordinate of the voronoi point
16. * side_ave1 is the average of all sides of all voronoi regions with the highest counted symmetry order
17. * side_std1 is standard deviation of all sides of all voronoi regions with the highest counted symmetry order

voronoi metrics contains the relevant data for checking translational symmetry breaking, distinguish positive from negative hexagons
"""
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
                                                         start_pos=start_pos_for_noise_corr_in_fspace)
    radial_plot1 -= correction(np.arange(0,len(radial_plot1)))
#    if plots:
#        plot2.plot(radial_plot1)

    if scaling_angle2 != 0:
        radial_plot2 = np.array([])
        for j in range(0,int(len(resft2)/2)):
            radial_plot2 = np.append(radial_plot2,circle_line_integration(resft2,j)[0])
        correction = gets_integration_noise_on_fourier_space(radial_plot2,
                                                             start_pos=start_pos_for_noise_corr_in_fspace)
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
            p1 = fit_ft_peak(wavevector_order=1, radial_spread=2, radial_plot=radial_plot1[:],
                             peaks_temp=peaks1_temp, fit='no_offset', plots=True, subplot=plot2)
        p1 = fit_ft_peak(wavevector_order=1, radial_spread=2, radial_plot=radial_plot1[:],
                         peaks_temp=peaks1_temp, fit='no_offset', plots=False)
        
        if scaling_angle2 != 0:
            peaks2_temp = find_peaks(radial_plot2,interpolation_points=1000,peak_finding_smoothness=5,
                                     plot=peak_plot, plot_new_fig=True)
            if plots:
                p2 = fit_ft_peak(wavevector_order=1, radial_spread=2, radial_plot=radial_plot2[:],
                                 peaks_temp=peaks2_temp, fit='no_offset', plots=True, subplot=plot4)
            p2 = fit_ft_peak(wavevector_order=1, radial_spread=2, radial_plot=radial_plot2[:],
                             peaks_temp=peaks2_temp, fit='no_offset', plots=False)
        
        if (np.size(p1) > 1):
            E1 = energy_ratio_wavevector_ring(1.,p=p1)
            if (scaling_angle2 != 0 and np.size(p2) > 1):
                E2 = energy_ratio_wavevector_ring(1.,p=p2)

            # The first condition will be the first order sideband (considering there is one) having
            # enough weigth, and the peak detection method not mess around the profile during normalisation.
            if (E1 > energy_ratio_condition and np.size(p1) > 1):
                
                # Collect the x-axis coordinate
                t = float(files1[i][start:stop])
                
                # Collect the real position of the peak
                Lambda1 = 1. / (p1[1] / (pixel_size/magnification*fft_size))
                # Collect the total intensity transmitted through the atoms
                trans_pow1 = radial_plot1[0] * I_cal
                # Collect all from the ring
                ring1_area = E1
                ring1_amplitude = p1[0]
                ring1_width = 2 * p1[2]
                symmetry_order_measured1,\
                azimuthal_peak1 = get_data_azimuthal_metrics(resft1, p1[1], which_sideband=0,
                                                             radial_epsilon=2,
                                                             interpolated_points=1000,
                                                             azimuthal_profile_smoothness=20,
                                                             plots=azimuthal_plots)
                count1, position_x1, position_y1, side_ave1,\
                side_std1 = get_data_voronoi_metrics(image1, symmetry_order_measured1, gaussian_sigma=10,
                                                     pixel_size=pixel_size, magnification=magnification,
                                                     plots=voronoi_plots)
                if (scaling_angle2 != 0):
                    if (np.size(p2) > 1):
                        Lambda2 = 1. / (p2[1] / (pixel_size/magnification*fft_size))
                        trans_pow2 = radial_plot2[0] * I_cal
                        ring2_area = E2
                        ring2_amplitude = p2[0]
                        ring2_width = 2 * p2[2]
                        return np.array([t, Lambda1, trans_pow1, ring1_area, ring1_amplitude, ring1_width,
                                         Lambda2, trans_pow2, ring2_area, ring2_amplitude, ring2_width,
                                         symmetry_order_measured1, azimuthal_peak1, count1, position_x1,
                                         position_y1, side_ave1, side_std1])
                    else:
                        return np.ones(18)*(-1)
                else:
                    return np.array([t, Lambda1, trans_pow1, ring1_area, ring1_amplitude, ring1_width,
                                     symmetry_order_measured1, azimuthal_peak1, count1, position_x1,
                                     position_y1, side_ave1, side_std1])
            elif scaling_angle2 != 0:
                return np.ones(18)*(-1)
            else:
                return np.ones(13)*(-1)
        elif scaling_angle2 != 0:
            return np.ones(18)*(-1)
        else:
            return np.ones(13)*(-1)
    elif scaling_angle2 != 0:
        return np.ones(18)*(-1)
    else:
        return np.ones(13)*(-1)