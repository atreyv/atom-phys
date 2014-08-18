from libphys import *

plt.ion()
plt.clf()

size_pz_ramp = 9000
single_mode_interval_ini = -0.015
single_mode_interval_end = 0.06

pzfile=scp.loadtxt('PZ-FP-2011.11.21-0001.dat')
pz=scp.transpose(pzfile)
peaksfile=scp.loadtxt('PZ-FP-2011.11.21-0000.dat')
peaks=scp.transpose(peaksfile)

peaks[1] *= 100.

fit = scp.polyfit(pz[0][:size_pz_ramp],pz[1][:size_pz_ramp],1)

fsr_peak_pos, peaks_temp = [], []
for i in range(0, size_pz_ramp):
    if (peaks[0][i] > single_mode_interval_ini and peaks[0][i] < single_mode_interval_end):
        if (peaks[1][i] > fit[0]*peaks[0][i]+fit[1]):
            peaks_temp += [peaks[1][i]]
#            print peaks[1][i], fit[0]*peaks[0][i]+fit[1]
        if (peaks_temp != [] and peaks[1][i] < fit[0]*peaks[0][i]+fit[1]):
            a = scp.argmax(peaks_temp)
            fsr_peak_pos += [i-len(peaks_temp)+a]
            peaks_temp = []

fsr_order = [0,1,2,3,4]
pz_volt=[]
for i in range(0, len(fsr_peak_pos)):
    pz_volt += [fit[0]*peaks[0][fsr_peak_pos[i]]+fit[1]]
pz_volt=scp.array(pz_volt)*100.
plt.scatter(fsr_order,pz_volt,label='PZ voltage Data')
plt.xlabel('Laser freq tuning (GHz)')
plt.ylabel('PZ Voltage (V)')

fit1=scp.polyfit(fsr_order,pz_volt,1,full=True)
fit2=scp.polyfit(fsr_order,pz_volt,2,full=True)
pz_volt1,pz_volt2=[],[]
for i in range(0,5):
	pz_volt1+=[fit1[0][0]*i+fit1[0][1]]
	pz_volt2+=[fit2[0][0]*i*i+fit2[0][1]*i+fit2[0][2]]
plt.plot(fsr_order,pz_volt1,label=r'1st order fit with $\chi^2 = $'+'%4.2f'%fit1[1][0])
plt.plot(fsr_order,pz_volt2,label=r'2nd order fit with $\chi^2 = $'+'%4.2f'%fit2[1][0]+' and $a_2 = $'+'%4.2f'%fit2[0][0])
#plt.text(2,52,r'$\chi^2$',fontsize=10)
plt.legend(loc=4)
