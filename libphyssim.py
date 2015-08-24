# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 16:09:13 2015

@author: pedro
"""

from libpatternsworkflow import *
from __future__ import division
from scipy.constants import k, hbar, u
from scipy.stats import maxwell

def alpha0L(b0, Delta):
    return b0 / (2 * (1 + 4 * Delta**2))

def F(alpha0L, Delta, n, w, n_normalize, w_normalize = 1):
    n = n / n_normalize
    w = w / w_normalize
    return np.exp(-alpha0L * (1 - 2j * Delta) * n * w)

def B(F, d, lambda_, R):
    F = ft.fft(F)
    q = np.arange(0,len(F)/2)
    q = np.append(q,np.arange(np.int(-len(F)/2),0))
#    F = np.absolute(F)
    return np.sqrt(R) * F * np.exp(1j * d * lambda_ * np.absolute(q)**2 / (4 * np.pi))

def s(F, B):
    B = ft.ifft(B)
#    B = np.absolute(B)
    return np.absolute(F)**2 + np.absolute(B)**2

def U_dip(s, delta):
    return hbar*delta/2 * np.log(1+s)

def v_dip(x, U_dip, M):
    return np.sqrt(2 / M * (U_dip[x-1] - np.amax(U_dip)))
    
    


def reset_grating(Lambda, number_atoms, x_bins, r):
    average_number_atoms = number_atoms / x_bins
    x = np.arange(1, x_bins + 1)
    n_model =  np.round( (1 + r * np.cos(2*np.pi*x/Lambda + np.pi)) * average_number_atoms)# * gaussian1d(1,800,800,0)(x)
    initial_grating = np.array([])
    for i in xrange(0,len(x)):
        initial_grating = np.append(initial_grating,np.ones(n_model[i])*x[i])
    n, bins, patches = plt.hist(initial_grating, x_bins)
    plt.scatter(x,n)
    return n, bins, initial_grating

def evolve_grating(grating, x_bins, M, T, dt):
    new_grating = grating.copy()
    v = np.sqrt(k*T/M)
    veloc = v * np.random.randn(len(grating))
    new_grating += veloc * dt
#    for i in xrange(0,len(grating)):
        # Maxwell-Boltzmann in 1D
#        v = np.sqrt(k*T/M)
#        veloc = v * np.random.randn(1)
#        new_grating[i] += veloc * dt #np.random.rand(1) * dt

        # Maxwell-Boltzmann in 3D
#        if grating[i] > x_bins/2:
#            new_grating[i] += maxwell.rvs(scale=veloc, size=1) * dt
#        else:
#            new_grating[i] -= maxwell.rvs(scale=veloc, size=1) * dt

    new_x_bins = np.amax(new_grating) - np.amin(new_grating) + 1
    plt.figure()
    n, bins, patches = plt.hist(new_grating, new_x_bins)
#    plt.plot(bins[1:],n)
    border = np.absolute(np.amin(new_grating))
    n = n[border:border+x_bins]
#    plt.figure()
#    plt.plot(n)
    return n, bins, new_grating

def bunching(n, Lambda_, x_bins, q_wavenumber_spread):
    """Define metric function - this case it will be the B = Sum(F{qc}) / F{0}"""
    B = ft.fft(n)
    #B = ft.fftshift(B)
    B = np.absolute(B)
#    plt.figure()
#    p = fit_ft_peak(1,1,b,peak,plots=False,fit='no_offset')
#    if np.size(p) > 1:
#        return 2*p[0]/b[0]
#    else:
#        return 0
#    b1 = B[:q_wavenumber_max]
#    peak = find_peaks_big_array(b1,len(b1)*10,10,True)
#    valleys_sorting = np.argsort(peak.peaks['valleys'][0])
#    peaks_sorting = np.argsort(peak.peaks['peaks'][0])
#    valley1 = peak.peaks['valleys'][0][valleys_sorting[0]]
#    peak_after_valley1 = peak.peaks['peaks'][0][peaks_sorting[0]]
#    p1 = np.amax(b1[valley1 : 2*peak_after_valley1-valley1])
#    b2 = B[len(B)-q_wavenumber_max:]
#    b2 = b2[::-1]
#    plt.plot(b2)
#    peak = find_peaks_big_array(b2,len(b2)*10,5,True)
#    valleys_sorting = np.argsort(peak.peaks['valleys'][0])
#    peaks_sorting = np.argsort(peak.peaks['peaks'][0])
#    valley1 = peak.peaks['valleys'][0][valleys_sorting[0]]
#    peak_after_valley1 = peak.peaks['peaks'][0][peaks_sorting[0]]
#    p2 = np.amax(b2[valley1 : 2*peak_after_valley1-valley1])
    peak_pos = x_bins / Lambda_
    p1 = np.amax(B[peak_pos - q_wavenumber_spread : peak_pos + q_wavenumber_spread])
    plt.figure()
    plt.plot(B[:q_wavenumber_max])
    return 2*(p1) / B[0]

def simulation(grating, x_bins, number_atoms, M, T, b0, Delta_, Lambda_, lambda_, l_, d_t, q_wavenumber_spread):
#    n_normalize = number_atoms / x_bins
    ballistic_motion = np.array([])
    for dt in d_t:
        n, bins, new_grating = evolve_grating(grating, x_bins, M, T, dt)
        n_normalize = np.average(n)
        a0l = alpha0L(b0, Delta_)
        F1 = F(a0l, Delta_, n, 1, n_normalize)
#        F1 = np.absolute(F1)**2
        F2 = B(F1, l_*Lambda_**2/(2*lambda_), lambda_, 1)
        # Cut higher spatial frequencies that q^2 > 3 * q_c^2
#        F2[15:len(F2)-15] = 0
        F2 = np.absolute(ft.ifft(F2))**2
        ballistic_motion = np.append(ballistic_motion, bunching(F2, Lambda_, x_bins, q_wavenumber_spread))
    return ballistic_motion