from libphys import *

plt.ion()
plt.clf()

i1file = scp.loadtxt('TDS-FP-20111214-i1', comments='#')
i2file = scp.loadtxt('TDS-FP-20111214-i2', comments='#')
i3file = scp.loadtxt('TDS-FP-20111214-i3', comments='#')

i1=scp.transpose(i1file)
i2=scp.transpose(i2file)
i3=scp.transpose(i3file)

i1[1]/=1e9
i2[1]/=1e9
i3[1]/=1e9

fit1 = scp.polyfit(i1[0],i1[1], 1)
fit2 = scp.polyfit(i2[0],i2[1], 1)
fit3 = scp.polyfit(i3[0],i3[1], 1)

print fit1, fit2, fit3

plt.xlabel('I (mA)')
plt.ylabel('Freq detuning (GHz)')

freq1, freq2, freq3 = [], [], []
 
for i in range(0,len(i1[0])):
    freq1 += [fit1[0]*i1[0][i]+fit1[1]]
for i in range(0,len(i2[0])):
    freq2 += [fit2[0]*i2[0][i]+fit2[1]]
for i in range(0,len(i3[0])):
    freq3 += [fit3[0]*i3[0][i]+fit3[1]]

plt.scatter(i1[0],i1[1],label='_nolegend_')
plt.plot(i1[0],freq1,label='slope = '+'%4.2f'%fit1[0]+r' $GHz/mA$')
plt.scatter(i2[0],i2[1],label='_nolegend_')
plt.plot(i2[0],freq2,label='slope = '+'%4.2f'%fit2[0]+r' $GHz/mA$')
plt.scatter(i3[0],i3[1],label='_nolegend_')
plt.plot(i3[0],freq3,label='slope = '+'%4.2f'%fit3[0]+r' $GHz/mA$')
plt.legend(loc=4)

