# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 14:04:08 2015

@author: pedro
"""

from libphys import *
from scipy.ndimage import zoom, rotate
from scipy.spatial import KDTree, Voronoi, voronoi_plot_2d
from scipy.ndimage import gaussian_filter

# scale_image for camera 448 is different from camera 113! because the images
# are mirrored.

def scale_image(image,angle,ratio_compression_y_over_x,interpolation_value=0):
    """Corrects for scaling distortions in image plane, e.g. from astigmatism
    in the optical system. Returns the scaled image with the same rotation."""
    image_rot = rotate(image,angle,order=interpolation_value)
    image_norm = zoom(image_rot,[ratio_compression_y_over_x,1],order=interpolation_value)
    image_final = rotate(image_norm,-angle,order=interpolation_value)
    return image_final


def circle_line(image,radius):
    """Returns the radial profile from an input image for a given radius"""
    if radius==0:
#        return image[len(image)/2,len(image)/2], 1
        return 0, 0
    if radius == 1:
        return image[np.int(len(image)/2)-1:np.int(len(image)/2)+2,np.int(len(image)/2)-1:np.int(len(image)\
            /2)+2].sum(), 9# - image[len(image)/2,len(image)/2], 9
    else:
        lx, ly = np.shape(image)
        x, y = np.ogrid[0:lx,0:ly]
        circle1 = (x-lx/2)**2 + (y-ly/2)**2 <= radius**2+1
        circle2 = (x-lx/2)**2 + (y-ly/2)**2 <= (radius-1)**2+1
#        image[circle1-circle2]=0
        return circle1,circle2


def get_azimuthal_profile_from_ft(ff,pos_actual):
    circle_circumscribed = pos_actual + 1
    circle1, circle2 = circle_line(ff,np.round(pos_actual))
    a = circle1^circle2
    y,x = a.shape
    b1 = a[np.int(y/2)-circle_circumscribed:np.int(y/2),np.int(x/2):np.int(x/2)+circle_circumscribed]
    b2 = a[np.int(y/2)-circle_circumscribed:np.int(y/2),np.int(x/2)-circle_circumscribed:np.int(x/2)]
    b3 = a[np.int(y/2):np.int(y/2)+circle_circumscribed,np.int(x/2)-circle_circumscribed:np.int(x/2)]
    b4 = a[np.int(y/2):np.int(y/2)+circle_circumscribed,np.int(x/2):np.int(x/2)+circle_circumscribed]
    ff_circle1 = (ff*a)[np.int(y/2)-circle_circumscribed:np.int(y/2),np.int(x/2):np.int(x/2)+circle_circumscribed]#.astype(np.int8)
    ff_circle2 = (ff*a)[np.int(y/2)-circle_circumscribed:np.int(y/2),np.int(x/2)-circle_circumscribed:np.int(x/2)]
    ff_circle3 = (ff*a)[np.int(y/2):np.int(y/2)+circle_circumscribed,np.int(x/2)-circle_circumscribed:np.int(x/2)]
    ff_circle4 = (ff*a)[np.int(y/2):np.int(y/2)+circle_circumscribed,np.int(x/2):np.int(x/2)+circle_circumscribed]
    ff_circle1_read = ff_circle1[::-1,::-1]
    ff_circle2_read = ff_circle2[::1,::-1]
    ff_circle3_read = ff_circle3
    ff_circle4_read = ff_circle4[::-1,::1]
    ff_profile = np.array([])
    ff_profile=np.append(ff_profile,ff_circle1_read[ff_circle1_read.astype(np.bool)])
    ff_profile=np.append(ff_profile,ff_circle2_read[ff_circle2_read.astype(np.bool)])
    ff_profile=np.append(ff_profile,ff_circle3_read[ff_circle3_read.astype(np.bool)])
    ff_profile=np.append(ff_profile,ff_circle4_read[ff_circle4_read.astype(np.bool)])
    angle = np.linspace(0,360,len(ff_profile))
    return angle, ff_profile

def get_azimuthal_profile_from_ft_integrated_along_radius(ff,pos_actual,radial_epsilon,interpolated_points):
    for pos_actual in range(pos_actual+radial_epsilon+1,pos_actual-radial_epsilon,-1):
        angle, ff_profile = get_azimuthal_profile_from_ft(ff,pos_actual)
        angle_interp = np.linspace(0,np.amax(angle),interpolated_points)
        ff_interp = np.interp(angle_interp,angle,ff_profile)
        try:
            ff_final += ff_interp
        except:
            ff_final = np.zeros(np.shape(ff_interp))
            ff_final += ff_interp
    return angle_interp, ff_final
    
def get_real_size(value,pixel_size,magnification,fft_size):
    return 1. / (value / (pixel_size/magnification*fft_size))

def fit_ft_peak(wavevector_order,radial_spread,radial_plot,peaks_temp,
                fit='no_offset',plots=False,subplot=plt, centre_tol=1.2):
    valleys_sorting = np.argsort(peaks_temp.peaks['valleys'][0])
    peaks_sorting = np.argsort(peaks_temp.peaks['peaks'][0])
    valley1 = peaks_temp.peaks['valleys'][0][valleys_sorting[2*wavevector_order-2]]
    valley2 = peaks_temp.peaks['valleys'][0][valleys_sorting[2*wavevector_order-1]]
#    valley1_actual = np.argmin(radial_plot[valley1-radial_spread:valley1+radial_spread+1]) + valley1-radial_spread# + 1
#    valley2_actual = np.argmin(radial_plot[valley2-radial_spread:valley2+radial_spread]) + valley2-radial_spread# + 1
    peak_after_valley1 = peaks_temp.peaks['peaks'][0][peaks_sorting[wavevector_order-1]]

    x_peak1 = (peaks_temp.x > np.round(peak_after_valley1) - radial_spread)
    x_peak2 = (peaks_temp.x < np.round(peak_after_valley1) + radial_spread)
    x_peak = x_peak1 * x_peak2
    if fit == 'no_offset' and x_peak.sum() > 3:
        p = fitgaussian1d_no_offset(peaks_temp.x[x_peak],peaks_temp.y_raw[x_peak])
        test = p[:3] < 0
    elif fit == 'offset' and x_peak.sum() > 4:
        x_peak3 = peaks_temp.x == valley1
        x_peak4 = peaks_temp.x == valley2
        x_peak = x_peak + x_peak3 + x_peak4
        p = fitgaussian1d(peaks_temp.x[x_peak],peaks_temp.y_raw[x_peak])
        test = p[:3] < 0
    else:
        print ('fit not understood or possible\n')
        return -1
#    plt.plot(peaks.x[x_peak],gaussian1d_no_offset(*p)(peaks.x[x_peak]))
    if test.any() == True:
        print ('fit output has negative values!... excluding.\n')
        return -1
    if p[1] > centre_tol * peak_after_valley1:
        print ('fit output has invalid centre!... excluding.\n')
        return -1
    if plots:
        v1 = np.round(valley1)
        v2 = np.round(valley2)
        x = np.linspace(v1,v2,(v2-v1)*10)
        if fit == 'no_offset':
            fit = gaussian1d_no_offset(*p)(x)
        else:
            fit = gaussian1d(*p)(x)
        subplot.plot(radial_plot[:np.int(np.round(2*v2-v1))])
        subplot.plot(x,fit)
    return p#, valley1_actual, valley2_actual

def energy_ratio_wavevector_ring(ref, p):
    A = p[0]
    x0 = p[1]
    sig = p[2]
    int_gauss_from_sig_minus_to_sig_plus =  A * np.sqrt(2 * np.pi) * sig * sp.special.erf(np.sqrt(2)/2)
    return int_gauss_from_sig_minus_to_sig_plus / np.sum(ref)

def get_pump_intensity_profile_from_txt(fname,beam_waist,intensity_plateau_n_points,
                                   smoothness_points,plot_all=False,
                                   check_plot=False,plot_find_peaks=False,
                                   averaging=True):
    file1 = np.loadtxt(fname)
    file1 = np.nan_to_num(file1)
    a = file1[:]
    a = 2*a / (np.pi * beam_waist**2)
    a -= np.amin(a)
    if plot_all:
        plt.figure()
        plt.plot(a)
    if averaging:        
        intensities = np.array(np.average(a[np.int(np.round(0.05*intensity_plateau_n_points)):np.int(np.round(intensity_plateau_n_points*0.98))]))
        if check_plot:
            plt.figure()
            plt.plot(a[np.int(np.round(0.05*intensity_plateau_n_points)):np.int(np.round(intensity_plateau_n_points*0.98))])
    else:
        intensities = np.array(np.amax(a[np.int(np.round(0.0*intensity_plateau_n_points)):np.int(np.round(intensity_plateau_n_points*0.98))]))
        if check_plot:
            plt.figure()
            plt.plot(a[np.int(np.round(0.0*intensity_plateau_n_points)):np.int(np.round(intensity_plateau_n_points*0.98))])
    a_extended = np.append(a,a[:smoothness_points])
    peaks = find_peaks_big_array(a_extended[:],len(a_extended[:])*1,smoothness_points,plot_find_peaks)
    peaks_pos = np.sort(peaks.peaks['peaks'][0])
    #print peaks_pos
    for i in peaks_pos:
        if averaging:
            intensities = np.append(intensities,np.average(a[np.int(np.round(i-(9*intensity_plateau_n_points/20))):
                                                         np.int(np.round(i+(12*intensity_plateau_n_points/25)))]))
        else:
            intensities = np.append(intensities,np.amax(a[np.int(np.round(i-(20*intensity_plateau_n_points/20))):
                                                         np.int(np.round(i+(25*intensity_plateau_n_points/25)))]))
    return intensities

def get_probe_intensity_profile_from_txt(fname,beam_waist,intensity_plateau_n_points,
                                   smoothness_points,plot_all=False,
                                   check_plot=False,plot_find_peaks=False):
    file1 = np.loadtxt(fname)
    file1 = np.nan_to_num(file1)
    a = file1[:]
    a = 2*a / (np.pi * beam_waist**2)
    a -= np.amin(a)
    if plot_all:
        plt.figure()
        plt.plot(a)
    if check_plot:
        plt.figure()
        plt.plot(a[np.int(np.round(1.02*intensity_plateau_n_points)):np.int(np.round(intensity_plateau_n_points*2))])
    intensities = np.array(np.amax(a[np.int(np.round(1.02*intensity_plateau_n_points)):np.int(np.round(intensity_plateau_n_points*2))]))
    a_extended = np.append(a,a[:intensity_plateau_n_points])
    peaks = find_peaks_big_array(a_extended[:],len(a_extended[:])*1,smoothness_points,plot_find_peaks)
#    peaks = find_peaks_big_array(a[:],len(a[:])*1,smoothness_points,plot_find_peaks)
    peaks_pos = np.sort(peaks.peaks['peaks'][0])
    #print peaks_pos
    for i in peaks_pos:
#        peaks_probe = find_peaks_big_array(a[i + intensity_plateau_n_points/2.2 :
#                                             i + 1.*intensity_plateau_n_points],
#                                           len(a[i + intensity_plateau_n_points/2.2 :
#                                                 i + 1.*intensity_plateau_n_points]),
#                                           5.,
#                                           plot_find_peaks)
#        peaks_probe_pos = np.sort(peaks_probe.peaks['peaks'][0])
#        intensities = np.append(intensities, np.amax(a[peaks_probe_pos[0]:peaks_probe_pos[0]+smoothness_points]))
                                           
        intensities = np.append(intensities,np.amax(a[np.int(np.round(i + intensity_plateau_n_points/1.95)) :
                                                      np.int(np.round(i + 1.*intensity_plateau_n_points))]))
    return intensities
