# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 16:50:55 2015

@author: pedro
"""

#from libphys import *
from libpatterns import *
from joblib import Parallel, delayed

def locate_roi_from_ref(dname, files, scaling_angle, raw_image, frac,
                        compression_y_over_x, image_crop_factor, plots=True):
    image = read_file_to_ndarray(dname+files[raw_image])
    param = use_ref_to_locate_centre(image)
    # Crop the image for FFT using the gaussian fit; the returned image is a square.
    ref = prepare_for_fft(image,param,frac)
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

