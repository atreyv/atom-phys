import numpy as np  # NumPy (multidimensional arrays, linear algebra, ...)
import scipy as sp  # SciPy (signal and image processing library)

import matplotlib as mpl         # Matplotlib (2D/3D plotting library)
#mpl.use('Agg')
import matplotlib.pyplot as plt  # Matplotlib's pyplot: MATLAB-like syntax
from pylab import *              # Matplotlib's pylab interface

from PIL import Image
import scipy.fftpack as ft
from scipy.optimize import leastsq as spleastsq
from matplotlib.colors import LogNorm
from scipy.ndimage import correlate as ndcorrelate
from scipy.ndimage import convolve as ndconvolve
from scipy.signal import convolve2d, correlate2d
from scipy.constants import k, u

# Customisations
mpl.rcParams['mathtext.fontset'] = 'stix'

# Turn on Matplotlib's interactive mode - in pylab
ion()


def gaussian_hwhm_to_radius(hwhm):
    return hwhm*np.sqrt(2/np.log(2))

def gaussian_radius_to_hwhm(radius):
    return radius*np.sqrt(np.log(2)/2)
    
def gaussian_radius_to_sig(radius):
    return radius / 2.

def gaussian_sig_to_radius(sig):
    return sig * 2.

def gaussian_hwhm_to_sig(hwhm):
    return hwhm/np.sqrt(2*np.log(2))

def gaussian_sig_to_hwhm(sig):
    return sig * np.sqrt(2*np.log(2))

#def gaussian1d(height, x0, hwhm, offset):
#    """Returns a function of the Gauss distribution for a given parameter set.\n
#    Example:\n
#    >>> g = gaussian1d(1,0,1.2,0.)\
#    x = np.linspace(0,2,100)\
#    plt.plot(x,g(x))"""
#    hwhm = float(hwhm)
#    return lambda x: height*np.exp(-1*((x-x0)/(hwhm))**2*np.log(2))\
#                    + offset

def gaussian1d(height, x0, sig, offset):
    """Returns a function of the Gauss distribution for a given parameter set.\n
    Example:\n
    >>> g = gaussian1d(1,0,1.2,0.)\
    x = np.linspace(0,2,100)\
    plt.plot(x,g(x))"""
    sig = float(sig)
    return lambda x: height*np.exp(-0.5*((x-x0)/sig)**2) + offset

#def residuals_gaussian1d(p,y,x):
#    """DEPRECATED\n
#    Calculates the array of residuals from Gaussian distribution"""
#    gmax, gx0, gfw, goffset = p
#    err = y - gmax*np.exp(-1*pow((x-gx0)/(gfw/2/np.log(2)),2)) - goffset
#    return err

def moments1d(x,data):
    """Returns (height, x0, stdev, offset) the gaussian parameters of 1D
    distribution found by a fit"""
    total = data.sum()
    if (x==None):
        x = np.arange(data.size)
    x0 = (x*data).sum()/total
    stdev = sqrt(((x-x0)**2*data).sum()/data.sum())
    height = np.amax(data)
    offset = np.amin(data)
    return height, x0, stdev, offset

def fitgaussian1d(x,data):
    """Returns (height, centre, sigma, offset)
    the gaussian parameters of a 1D distribution found by a fit"""
    params = moments1d(x,data)
    if (x==None):
        errorfunction = lambda p: gaussian1d(*p)(*np.indices(data.shape))\
                                - data
    else:
        errorfunction = lambda p: gaussian1d(*p)(x)\
                                - data
    p, success = spleastsq(errorfunction, params, full_output=0)
    return p



#def gaussian2d(height, x0, y0, hwhm_x, hwhm_y,offset):
#    """Returns a gaussian function with the given parameters"""
#    hwhm_x = float(hwhm_x)
#    hwhm_y = float(hwhm_y)
#    return lambda x,y: height*np.exp(
#                 -(((x-x0)/hwhm_x)**2+((y-y0)/hwhm_y)**2)*np.log(2))+offset

def gaussian2d(height, x0, y0, sig_x, sig_y,offset):
    """Returns a gaussian function with the given parameters:
    (height, y, x, sig_y, sig_x, offset)"""
    sig_x = float(sig_x)
    sig_y = float(sig_y)
    return lambda x,y: height*np.exp(
                 -0.5*(((x-x0)/sig_x)**2+((y-y0)/sig_y)**2))+offset


def moments2d(data):
    """Receives numpy array data with dim=2 and 
    returns (height, y, x, sig_y, sig_x, offset)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    sig_x = sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    sig_y = sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = np.nanmax(data)
    offset = np.nanmin(data)
    return height, x, y, sig_x, sig_y, offset

def fitgaussian2d(data):
    """Returns (height, y, x, sig_y, sig_x, offset)
    the gaussian parameters of a 2D distribution found by a fit"""
#    data = np.transpose(data)    
    params = moments2d(data)
    errorfunction = lambda p: np.ravel(gaussian2d(*p)(*np.indices(data.shape)) -
                                 data)
    p, success = spleastsq(errorfunction, params, xtol=1e-16,ftol=1e-16)
    return p


def lorentz1d(height, x0, fwhm, offset):
    """Returns a function of the Lorentzian distribution for a given parameter set.\n
    Example:\n
    >>> lor = lorentz1d(1,0,1.2,0.)\
    x = np.linspace(0,2,100)\
    plt.plot(x,lor(x))"""
    fwhm = float(fwhm)    
    return lambda x: height/(1 + (2*(x-x0)/fwhm)**2) + offset

def residuals_lorentz(p,y,x):
    """DEPRECATED\n
    Calculates the array of residuals from Lorentz distribution"""
    lmax, lx0, lfw, loffset = p
    err = y - lmax/(1 + (2*(x-lx0)/lfw)**2) - loffset
    return err

def fitlorentz1d(x,data):
    """Returns (height, centre, fwhm, offset)
    the lorentzian parameters of a 1D distribution found by a fit"""
    params = moments1d(x,data)
    if (x==None):
        errorfunction = lambda p: lorentz1d(*p)(*np.indices(data.shape))\
                                - data
    else:
        errorfunction = lambda p: lorentz1d(*p)(x)\
                                - data        
    p, success = spleastsq(errorfunction, params, full_output=0)
    return p


def extinction_lorentz(b0,nu0,nu=1):
    """Function to calculate the transmission through a medium with a
    specific thickness. Receives b0, optical thickness at resonance,
    nu0, the offset from zero in centre frequency and offset in value
    from zero transmission, normally derived from laser linewidth."""
#    fwhm = float(fwhm)
#    nu = 6.066
    return lambda x: np.exp(-b0/(1 + (2*(x-nu0)/nu)**2))

def moments_extinction_lorentz(x,data):
    """Returns (b0, nu0, offset) the moments of transmission
    distribution for a fit"""
#    total = data.sum()
    if (x==None):
        x = np.arange(data.size)
    nu0 = x[np.argmin(data)]
    b0 = sqrt(((x-nu0)**2*data).sum()/data.sum())
#    offset = np.amin(data)
    return b0, nu0

def fit_extinction_lorentz(x,data):
    """Returns the (b0, nu0, offset) exponential decay with a lorentzian
    argument (b(nu)) parameters of a distribution found by a fit"""
    params = moments_extinction_lorentz(x,data)
#    params = np.array([8,0])
    if (x==None):
        errorfunction = lambda p: extinction_lorentz(*p)(*np.indices(data.shape))\
                                - data
    else:
        errorfunction = lambda p: extinction_lorentz(*p)(x)\
                                - data        
    p, success = spleastsq(func=errorfunction, x0=params)#, xtol=1e-16,ftol=1e-16)
    return p

    
def kinetic_expansion(Temp, sigma0):
    """Returns a function for time-of-flight measurements"""
    return lambda t: np.sqrt(Temp*t**2 + sigma0**2)

def moments_tof(x,data):
    """Calculates the initial parameters for a fit of time-of-flight to data"""
    sigma0 = np.amin(data)
    Temp = np.average((data**2-sigma0**2)/x**2)
    return Temp, sigma0

def fit_tof(x,data):
    """Returns the sigma0 and Temp for a time-of-flight measurements"""
    params = moments_tof(x,data)
    if (x==None):
        errorfunction = lambda p: kinetic_expansion(*p)(*np.indices(data.shape)) - data
    else:
        errorfunction = lambda p: kinetic_expansion(*p)(x) - data
    p, success = spleastsq(func=errorfunction, x0=params)#, xtol=1e-16,ftol=1e-16)
    return p



def low_pass_rfft(curve, low_freqs):
    """Filters the curve by setting to zero the high frequencies"""
    a = ft.rfft(curve)
    for i in range(2*low_freqs, len(curve)):
        a[i]=0    
    return np.array(ft.irfft(a))

def FourierFilter(function, half_interval):
    """Returns a fourier space filtered function by setting to zero all
    frequencies above half-interval and below -half-interval"""
    f_fft=ft.fft(function)
    for j in range(0,len(f_fft)):
        if(j>half_interval and j<len(f_fft)-half_interval):
            f_fft[j] = 0
    return ft.ifft(f_fft)



def prepare_for_fft(input_image,fft_size,image_centre):
    """Returns an image cropped around image_centre with size fft_size.
    Image_centre should be a tuple with image centre coordinates
    or zero if centre should be found"""
    x,y = input_image.shape
    if image_centre == 0:
        centre_x, centre_y = unravel_index(input_image.argmax(), input_image.shape)        
    if image_centre != 0:
        centre_x,centre_y = image_centre
    if (x - centre_x < fft_size/2 or y - centre_y < fft_size/2):
        print "FFT size is bigger than the image itself!"
        return -1
    return input_image[centre_x-fft_size/2:centre_x+fft_size/2,\
        centre_y-fft_size/2:centre_y+fft_size/2]

def circle_line_integration(image,radius):
    """Calculates the integral in a radial perimeter with input radius
    on input image and returns integral and pixels in integral"""
    if radius==0:
#        return image[len(image)/2,len(image)/2], 1
        return 0, 0
    if radius == 1:
        return image[len(image)/2-1:len(image)/2+2,len(image)/2-1:len(image)\
            /2+2].sum() - image[len(image)/2,len(image)/2], 9
    else:
        lx, ly = shape(image)
        x, y = np.ogrid[0:lx,0:ly]
        circle1 = (x-lx/2)**2 + (y-ly/2)**2 <= radius**2+1
        circle2 = (x-lx/2)**2 + (y-ly/2)**2 <= (radius-1)**2+1
#        image[circle1-circle2]=0
        return image[circle1-circle2].sum(), (circle1-circle2).sum()
    
def normalize_by_division(signal_image,ref_image):
    """Receives two images, the first one is the signal and the second
    the reference (e.g. gaussian spatial profile of a pump beam). 
    Divides the first by the second and removes any inf or nan in 
    the resultant matrix"""
    signal = signal_image / ref_image
    for pos in np.nditer(signal,op_flags=['readwrite']):
        if(pos==inf):
            pos[...] = 1.
        if(pos==-inf):
            pos[...] = 0
        if(pos < 0):
            pos[...] = 0
    signal = np.nan_to_num(signal)
    return signal
