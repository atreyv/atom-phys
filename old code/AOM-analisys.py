from libphys import *

plt.ion()

laser_power_lm_0301 = 9.97
laser_power_pm = 10.40
laser_power_un = 9.62
laser_power_lm_0305 = 12.30


#Load data
lm0 = scp.transpose(scp.loadtxt('LM110_20120301_00.txt'))
lm1 = scp.transpose(scp.loadtxt('LM110_20120301_01.txt'))
pm0 = scp.transpose(scp.loadtxt('PM080_20120301_00.txt'))
pm1 = scp.transpose(scp.loadtxt('PM080_20120301_01.txt'))
un0 = scp.transpose(scp.loadtxt('unknown65_20120301_00.txt'))
un1 = scp.transpose(scp.loadtxt('unknown65_20120301_01.txt'))
lm050 = scp.transpose(scp.loadtxt('LM110_20120305_00.txt'))
lm051 = scp.transpose(scp.loadtxt('LM110_20120305_01.txt'))
lm052 = scp.transpose(scp.loadtxt('LM110_20120305_02.txt'))
lm053 = scp.transpose(scp.loadtxt('LM110_20120305_03.txt'))


#Process Data
lm0[1] /= laser_power_lm_0301
pm0[1] /= laser_power_pm
un0[1] /= laser_power_un
lm051[1] /= laser_power_lm_0305
pm1[1] /= laser_power_pm
un1[1] /= laser_power_un
lm050[1] /= laser_power_lm_0301
lm052[1] /= laser_power_lm_0305
lm053[1] /= laser_power_lm_0305
lm1[1] /= laser_power_lm_0301


#Fixed RF freq, Diffrated beam eff vs RF power
plt.clf()
plt.xlabel('RF Power (W)')
plt.ylabel('Diffracted beam efficiency')

plt.scatter(lm0[0],lm0[1],label='LM110 @110MHz')
plt.scatter(pm0[0],pm0[1],c='g',label='PM080 @80MHz')
plt.scatter(un0[0],un0[1],c='r',label='unknown @80MHz')
plt.legend(loc=4)

plt.savefig('Fixed RF freq, Diffracted beam eff vs RF power for 3 Isle Optics AOM.eps')


#Fixed RF power, Diffrated beam eff vs RF freq
plt.clf()
plt.xlabel('RF freq (MHz)')
plt.ylabel('Diffracted beam efficiency')

plt.scatter(lm051[0],lm051[1],label='LM110 @ 0.5W')
plt.scatter(pm1[0],pm1[1],c='g',label='PM080 @ 0.5W')
plt.scatter(un1[0],un1[1],c='r',label='unknown @ 1W')
plt.legend(loc=1)

plt.savefig('Fixed RF power, Diffracted beam eff vs RF freq for 3 Isle Optics AOM.eps')


#Fixed RF freq, Diffracted beam eff vs RF power for LM110
plt.clf()
plt.xlabel('RF Power (W)')
plt.ylabel('Diffracted beam efficiency')

plt.scatter(lm050[0],lm050[1],label='LM110 @110MHz')
plt.scatter(lm052[0],lm052[1],c='g',label='LM110 @90MHz')
plt.scatter(lm053[0],lm053[1],c='r',label='LM110 @95MHz')
plt.legend(loc=4)

plt.savefig('3 diff fixed RF freq for LM110, Diffracted beam eff vs RF power.eps')


#Fixed RF power, Diffracted beam eff vs RF freq for LM110
plt.clf()
plt.xlabel('RF freq (MHz)')
plt.ylabel('Diffracted beam efficiency')

plt.scatter(lm051[0],lm051[1],label='LM110 @ 0.5W')
plt.scatter(lm1[0],lm1[1],c='g',label='LM110 @ 1W')
plt.legend(loc=1)

plt.savefig('2 diff fixed RF power for LM110, Diffracted beam eff vs RF freq.eps')
