
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

dsave = '/home/pedro/Dropbox/PhD at Strathclyde/Thesis/fig/'
dname_ortho = '/home/pedro/LAB/DATA/2014/Dec/PF/05_12_14/pf01/cam448/'
dname_parallel = '/home/pedro/LAB/DATA/2014/Dec/PF/05_12_14/pf01/cam113/'
iname = '/home/pedro/LAB/DATA/2014/Dec/PF/05_12_14/pf01/'
files_ortho = load_files(dname_ortho,".bmp")
files_parallel = load_files(dname_parallel,".bmp")
files_int = load_files_prefix(iname, prefix='tpump', sufix='.dat')


# ####Set the parameters with the information about which images contain patterns (fstart and fstop) and the position in filename to extract the relevant x-axis variable (start and stop)

# In[ ]:

fstart = 20
fstop = 880
raw_image = 100
start = 5
stop = 9
(files_parallel[600][start:stop])


# ###Set some of experimental parameters.
# 
# * frac is the fraction of gaussian sigma got from 2D fit to ref, to set the size in FFT

# In[ ]:

pixel_size = 3.75 # micrometers
magnification = 2.
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
start_pos_for_noise_corr_in_fspace = 30
energy_ratio_condition = 0.0
# The next parameters affect the power logged in a text file
intensity_plateau_n_points = 1000
int_smoothness = 300


# ###Use a reference image to locate the centre of pump/probe for ortho/parallel pol

# For the orthogonal polarisation; notice that a pure reference is a flat plane intensity profile for linear polarisation at the Brewster angle, so a typical pattern is selected

# In[ ]:

fft_size_ortho, param_ortho, refft_ortho = locate_roi_from_ref(dname_ortho, files_ortho, scaling_angle_ortho,
                                                               raw_image, frac, compression_y_over_x,
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

P0_ref_ccd = np.array([calibrate_intensity(raw_image)])
for i in xrange(raw_image+1,raw_image+20):
    P0_ref_ccd = np.append(P0_ref_ccd, calibrate_intensity(i))
    
I0_ref_pd = get_pump_intensity_profile_from_txt(iname+'ref.dat',0.08,intensity_plateau_n_points,
                                                 smoothness_points=int_smoothness,plot_all=True,
                                                 check_plot=True, plot_find_peaks=True,averaging=False)

I0_pump_cal = np.average(I0_ref_pd / P0_ref_ccd)
I0_pump = np.average(P0_ref_ccd * I0_pump_cal)

I0_pump_pd = get_pump_intensity_profile_from_txt(iname+files_int[0],0.08,intensity_plateau_n_points,
                                                 int_smoothness,True,True,True,averaging=False)
#I0_probe_pd = get_probe_intensity_profile_from_txt(iname+files_int[0],0.08,intensity_plateau_n_points,
#                                                   int_smoothness,True,True,True)
for i in files_int[1:]:
    I0_pump_pd = np.append(I0_pump_pd,get_pump_intensity_profile_from_txt(iname+i,0.08,
                                                                          intensity_plateau_n_points,
                                                                          int_smoothness,True,True,False,
                                                                          averaging=False))
#    I0_probe_pd = np.append(I0_probe_pd,get_probe_intensity_profile_from_txt(iname+i,0.08,
#                                                                             intensity_plateau_n_points,
#                                                                             int_smoothness,True,True,False))


# In[ ]:

plt.figure(figsize=(12,8))
plt.plot((P0_ref_ccd-np.amin(P0_ref_ccd))/(np.amax(P0_ref_ccd)-np.amin(P0_ref_ccd)) * 100, marker='o', linewidth=1)
plt.plot((I0_ref_pd-np.amin(I0_ref_pd))/(np.amax(I0_ref_pd)-np.amin(I0_ref_pd)) * 100,
         marker='x', linewidth=1)
plt.xlabel('acquisition number')
plt.ylabel('% $\Delta$')
#plt.savefig(iname+'intensity_green-image_power_blue-relation.pdf')


# ###Call the main function that will collect all data from each image and append it to the created container
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
# 11. * symmetry_order_measured1 is the counted number of sideband peaks
# 12. * azimuthal_peak1 is the azimuthal angle of the first peak from quadrant 1
# 13. * count1 is the number of voronoi regions with the counted symmetry order
# 14. * position_x1 is the c.o.m. in the x-coordinate of the voronoi point
# 15. * position_y1 is the c.o.m. in the y-coordinate of the voronoi point
# 16. * side_ave1 is the average of all sides of all voronoi regions with the highest counted symmetry order
# 17. * side_std1 is standard deviation of all sides of all voronoi regions with the highest counted symmetry order
# 
# voronoi metrics contains the relevant data for checking translational symmetry breaking, distinguish positive from negative hexagons

# Write the main function for looping. It might be used with joblib parallelisation. In this case the plots don't show in the end, although the call for plot wastes processor time.

# def main():
#     data = [get_data_metrics(600,
#                              start,
#                              stop,
#                              pixel_size,
#                              magnification,
#                              raw_image,
#                              frac,
#                              symmetry_order,
#                              azimuthal_profile_smoothness,
#                              compression_y_over_x,
#                              image_crop_factor,                    
#                              start_pos_for_noise_corr_in_fspace,
#                              energy_ratio_condition,
#                              intensity_plateau_n_points,
#                              int_smoothness,
#                              fft_size,
#                              param_ortho, scaling_angle_ortho,
#                              dname_ortho, files_ortho, I0_pump_cal, True, False, True, True,
#                              param_parallel, scaling_angle_parallel, dname_parallel, files_parallel)]
#     for i in xrange(600+1,640):
#         data = np.append(data, [get_data_metrics(i,
#                                                  start,
#                                                  stop,
#                                                  pixel_size,
#                                                  magnification,
#                                                  raw_image,
#                                                  frac,
#                                                  symmetry_order,
#                                                  azimuthal_profile_smoothness,
#                                                  compression_y_over_x,
#                                                  image_crop_factor,                    
#                                                  start_pos_for_noise_corr_in_fspace,
#                                                  energy_ratio_condition,
#                                                  intensity_plateau_n_points,
#                                                  int_smoothness,
#                                                  fft_size,
#                                                  param_ortho,
#                                                  scaling_angle_ortho, dname_ortho, files_ortho, I0_pump_cal,
#                                                  True, False, True, True, param_parallel,
#                                                  scaling_angle_parallel,
#                                                  dname_parallel, files_parallel)],axis=0)
#     return data
# %time data = main()

# In[ ]:

get_ipython().magic(u'time data = np.asarray(Parallel(n_jobs=4)(delayed(get_data_metrics)(i,                                                                     start,                                                                     stop,                                                                     pixel_size,                                                                     magnification,                                                                     raw_image,                                                                     frac,                                                                     symmetry_order,                                                                     azimuthal_profile_smoothness,                                                                     compression_y_over_x,                                                                     image_crop_factor,                                                                     start_pos_for_noise_corr_in_fspace,                                                                     energy_ratio_condition,                                                                     intensity_plateau_n_points,                                                                     int_smoothness,                                                                     fft_size,                                                                     param_ortho,                                                                     scaling_angle_ortho, dname_ortho,                                                                     files_ortho,                                                                     I0_pump_cal, False, False, False, False,                                                                     param_parallel,                                                                     scaling_angle_parallel, dname_parallel,                                                                     files_parallel)                                           for i in range(fstart,fstop)))')


# In[ ]:

pump = I0_pump_pd
#probe = I0_probe_pd[20:]
plt.plot(pump)
#plt.plot(probe)
print "Ipump = %.2f mW/cm2 with std = %.2f mW/cm2"%(np.average(pump), np.std(pump))
#print "Iprobe = %.2f mW/cm2 with std = %.2f mW/cm2"%(np.average(probe), np.std(probe))
print I0_pump_pd.shape


# Correct desired array for invalid metrics...

# patched_data9 = np.copy(data[:,9])
# for i in xrange(3):
#     element = np.argmax(patched_data9)
#     patched_data9[element] = patched_data9[element-1]

# In[ ]:

rejected_truth = data[:,0] != -1
accept_interval_in_std = 2
truth_int = rejected_truth#truth_intensities(I0_pump_pd,imin=0,imax=300)

#t_value, total_int_mean, total_int_errors,\
#total_int_count = average_std_data(data[:,0], I0_pump_pd, accept_interval_in_std, truth_int)

t_value, value1_mean_ortho, value1_errors_ortho,count1_ortho = average_std_data(data[:,0], data[:,2], accept_interval_in_std, truth_int)

t_value, value2_mean_ortho, value2_errors_ortho,count2_ortho = average_std_data(data[:,0], data[:,4], accept_interval_in_std, truth_int)

t_value, value1_mean_parallel, value1_errors_parallel,count1_parallel = average_std_data(data[:,0], data[:,7], accept_interval_in_std, truth_int)

t_value, value2_mean_parallel, value2_errors_parallel,count2_parallel = average_std_data(data[:,0], data[:,9], accept_interval_in_std, truth_int)

value2_mean_ortho = value2_mean_ortho / ((value1_mean_ortho+value1_mean_parallel) / I0_pump_cal)
value2_errors_ortho = value2_errors_ortho / ((value1_mean_ortho+value1_mean_parallel) / I0_pump_cal)
value2_mean_parallel = value2_mean_parallel / ((value1_mean_ortho+value1_mean_parallel) / I0_pump_cal)
value2_errors_parallel = value2_errors_parallel / ((value1_mean_ortho+value1_mean_parallel) / I0_pump_cal)


# In[ ]:

count1_ortho, count2_ortho, count1_parallel, count2_parallel


# In[ ]:

fig = plt.figure(figsize=(8,6))
plot1 = fig.add_subplot(111)
plot1.plot(t_value,value1_mean_ortho,linewidth=0,marker='o')
plot1.errorbar(t_value,value1_mean_ortho,yerr=value1_errors_ortho,ls='none')
plot1.set_xlabel(r't ($\mathrm{\mu s}$)')
#plot1.set_ylabel(r'$\Lambda$ ($\mathrm{\mu m}$)')
plot1.set_ylabel(r'$I_{trans}$ ($mW/cm^2$)')
plot1.set_xlim(xmin=0)#,xmax=250)
plot1.set_ylim(ymin=0)
plt.savefig(dsave+'decay-02-12_12_14-trans-total.pdf')


# In[ ]:

fig = plt.figure(figsize=(8,6))
plot2 = fig.add_subplot(111)
plot2.plot(t_value,value2_mean_ortho,linewidth=0,marker='o')
plot2.errorbar(t_value,value2_mean_ortho,yerr=value2_errors_ortho,ls='none')
plot2.set_xlabel(r't ($\mathrm{\mu s}$)')
plot2.set_ylabel(r'$E_{q_1}$/$E_{0}$')
plot2.set_xlim(xmin=0)#,xmax=250)
plot2.set_ylim(ymin=0)#,ymax=1)


# In[ ]:

fig = plt.figure(figsize=(8,6))
plot3 = fig.add_subplot(111)
plot3.plot(t_value,value1_mean_parallel,linewidth=0,marker='o')
plot3.errorbar(t_value,value1_mean_parallel,yerr=value1_errors_parallel,ls='none')
#plot3.plot(t_value,total_int_mean,linewidth=0,marker='d')
#plot3.errorbar(t_value,total_int_mean,yerr=total_int_errors,ls='none')
plot3.set_xlabel(r't ($\mathrm{\mu s}$)')
#plot3.set_ylabel(r'$\Lambda$ ($\mathrm{\mu m}$)')
plot3.set_ylabel(r'$I_{trans}$ ($mW/cm^2$)')
plot3.set_xlim(xmin=0)#,xmax=250)
plot3.set_ylim(ymin=0)


# In[ ]:

fig = plt.figure(figsize=(8,6))
plot4 = fig.add_subplot(111)
plot4.plot(t_value,value2_mean_parallel,linewidth=0,marker='o')
plot4.errorbar(t_value,value2_mean_parallel,yerr=value2_errors_parallel,ls='none')
#plot4.plot(t_value,intense/np.average(intense))
plot4.set_xlabel(r't ($\mathrm{\mu s}$)')
plot4.set_ylabel(r'$E_{q_1}$/$E_{0}$')
plot4.set_xlim(xmin=0)#,xmax=250)
plot4.set_ylim(ymin=0)#,ymax=1)
#fig.savefig(dname+'ortho-postmirror-pol.pdf')


# In[ ]:

bincount = np.bincount(np.round(data[:,12][rejected_truth]).astype(np.int64))
bincount.argmax(), bincount


# In[ ]:

sym = np.bincount(np.round(data[:,11][rejected_truth]).astype(np.int64))
sym_truth = data[:,11] == sym.argmax()
sym.argmax(), sym


# In[ ]:

rot_truth = rejected_truth * sym_truth
rotational = data[:,12][rot_truth]
n, bins, patches = plt.hist(rotational,70,range=(-5,65),histtype='bar')
plt.xlabel('azimuthal angle (deg)')
plt.ylabel('counts')


# In[ ]:

azi_truth1 = data[:,12] >= bincount.argmax() - 1.5
azi_truth2 = data[:,12] < bincount.argmax() + 1.5
trans_truth = rejected_truth * sym_truth * azi_truth1 * azi_truth2

first_pattern = np.sqrt(data[:,14][trans_truth][0]**2 + data[:,15][trans_truth][0]**2)
translational = (np.sqrt(data[:,14][trans_truth]**2 + data[:,15][trans_truth]**2) - first_pattern)                * pixel_size/magnification

fig, ax1 = plt.subplots()
ax1.set_ylabel('angle (deg)')
ax1.set_ylim(ymin=0)#,ymax=60)
ax1.set_xlabel('image number')
for tl in ax1.get_yticklabels():
    tl.set_color('r')
ax2 = ax1.twinx()
ax2.set_ylabel('absolute translation ($\mathrm{\mu m}$)')
for tl in ax2.get_yticklabels():
    tl.set_color('b')
ax1.plot(data[:,12][trans_truth],'r--')
ax2.plot(translational, 'b')


# In[ ]:

n, bins, patches = plt.hist(translational,80,range=(-40,40),histtype='bar')
plt.xlabel('absolute translation ($\mathrm{\mu m}$)')
plt.ylabel('counts')

