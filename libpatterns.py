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
        return image[len(image)/2-1:len(image)/2+2,len(image)/2-1:len(image)\
            /2+2].sum(), 9# - image[len(image)/2,len(image)/2], 9
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
    a = circle1-circle2
    y,x = a.shape
    b1 = a[y/2-circle_circumscribed:y/2,x/2:x/2+circle_circumscribed]
    b2 = a[y/2-circle_circumscribed:y/2,x/2-circle_circumscribed:x/2]
    b3 = a[y/2:y/2+circle_circumscribed,x/2-circle_circumscribed:x/2]
    b4 = a[y/2:y/2+circle_circumscribed,x/2:x/2+circle_circumscribed]
    ff_circle1 = (ff*a)[y/2-circle_circumscribed:y/2,x/2:x/2+circle_circumscribed]#.astype(np.int8)
    ff_circle2 = (ff*a)[y/2-circle_circumscribed:y/2,x/2-circle_circumscribed:x/2]
    ff_circle3 = (ff*a)[y/2:y/2+circle_circumscribed,x/2-circle_circumscribed:x/2]
    ff_circle4 = (ff*a)[y/2:y/2+circle_circumscribed,x/2:x/2+circle_circumscribed]
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
    for pos_actual in range(pos_actual+radial_epsilon,pos_actual-radial_epsilon-1,-1):
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
                fit='no_offset',plots=False,subplot=plt):
    valleys_sorting = np.argsort(peaks_temp.peaks['valleys'][0])
    peaks_sorting = np.argsort(peaks_temp.peaks['peaks'][0])
    valley1 = peaks_temp.peaks['valleys'][0][valleys_sorting[2*wavevector_order-2]]
    #valley2 = peaks_temp.peaks['valleys'][0][valleys_sorting[2*wavevector_order-1]]
    valley1_actual = np.argmin(radial_plot[valley1-radial_spread:valley1+radial_spread]) + valley1-radial_spread# + 1
    #valley2_actual = np.argmin(radial_plot[valley2-radial_spread:valley2+radial_spread]) + valley2-radial_spread# + 1
    peak_after_valley1 = peaks_temp.peaks['peaks'][0][peaks_sorting[wavevector_order-1]]
    valley2_actual = peak_after_valley1 + (peak_after_valley1 - valley1_actual)
    if fit == 'no_offset':
        p = fitgaussian1d_no_offset(None,radial_plot[np.round(valley1_actual):np.round(valley2_actual)])
    elif fit == 'offset':
        p = fitgaussian1d(None,radial_plot[np.round(valley1_actual):np.round(valley2_actual)])
    else:
        print 'fit not understood'
        return -1
    p[1] += np.round(valley1_actual)
    if plots:
        v1 = np.round(valley1_actual)
        v2 = np.round(valley2_actual)
        x = np.linspace(v1,v2,(v2-v1)*10)
        if fit == 'no_offset':
            fit = gaussian1d_no_offset(*p)(x)
        else:
            fit = gaussian1d(*p)(x)
        subplot.plot(radial_plot[:2*v2-v1])
        subplot.plot(x,fit)
    return p#, valley1_actual, valley2_actual

def energy_ratio_wavevector_ring(ref, p):
    A = p[0]
    x0 = p[1]
    sig = p[2]
#    B = p[3]
    int_gauss_from_sig_minus_to_sig_plus =  A * np.sqrt(2 * np.pi) * sig * sp.special.erf(np.sqrt(2)/2 * sig**2)# + 2 * B * sig
    return int_gauss_from_sig_minus_to_sig_plus / np.sum(ref)
