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
    if sorted2[symmetry_order_measured-1] > interpolated_points:
        sorted2[symmetry_order_measured-1] -= interpolated_points
        sorted2 = np.sort(sorted2)
    azimuthal_peaks = azimuth_interp[sorted2.astype(np.int16)]
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