# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 12:06:56 2014

@author: pedro
"""

from libphys import *
mpl.rcParams['figure.figsize'] = (16.0, 8.0)

#dname1 = '/home/pedro/LAB/DATA/2014/Dec/PF/05_12_14/pf03/cam113/'
#dname2 = '/home/pedro/LAB/DATA/2014/Dec/PF/05_12_14/pf04/cam113/'
#dname3 = '/home/pedro/LAB/DATA/2014/Dec/PF/05_12_14/pf01/cam113/'
#f1 = 'pattern-parallel_pump_dyn-circpol-no_molasses.dat'
#f2 = 'pattern-parallel_pump_dyn-circpol-molasses.dat'
#f3 = 'pattern-parallel_pump_dyn-linpol_no-molasses.dat'
#nomol = np.loadtxt(dname1+f1, skiprows=1)
#mol = np.loadtxt(dname2+f2, skiprows=1)
#lin = np.loadtxt(dname3+f3, skiprows=1)
#nomol = np.transpose(nomol)
#mol = np.transpose(mol)
#lin = np.transpose(lin)

time_nomol = np.array([0,0,0,0,1,1,1,1,2,2,2,2,\
                    4,4,5,6,6,7,8,8,9,10,\
                    ]) 
b0_nomol = np.array([34.3,34.0,33.3,33.5,31.5,32.3,30.8,31.2,29.1,28.3,29.0,27.3,\
                    19.1,21.4,15.2,12.5,14.3,10.9,9.0,9.0,7.7,6.6,\
                    ])
time_mol = np.array([4.1,4.5,5,5.5,6,7,8,9,10])
b0_mol = np.array([20.1,18.7,19.3,16.9,15.2,15.5,14.6,14.8,14])

tnomol_value = np.array([])
b0nomol_mean = np.array([])
b0nomol_errors = np.array([])
tmol_value = np.array([])
b0mol_mean = np.array([])
b0mol_errors = np.array([])

j = np.inf
for i in time_nomol:
    if i != j:
        temp = time_nomol==i
        tnomol_value = np.append(tnomol_value,i)
        b0nomol_mean = np.append(b0nomol_mean, np.average(b0_nomol[temp==True]))
        b0nomol_errors = np.append(b0nomol_errors, np.std(b0_nomol[temp==True]))
        j = i
    else:
        pass

j = np.inf
for i in time_mol:
    if i != j:
        temp = time_mol==i
        tmol_value = np.append(tmol_value,i)
        b0mol_mean = np.append(b0mol_mean, np.average(b0_mol[temp==True]))
        b0mol_errors = np.append(b0mol_errors, np.std(b0_mol[temp==True]))
        j = i
    else:
        pass

fig = plt.figure()
#fig.subplots_adjust(left=0.1,right=0.9)
ax1 = fig.add_subplot(111)
ax1.set_ylabel(r'$b_0$')

#ax2 = fig.add_subplot(212)
ax1.set_xlabel('t (ms)')


ax1.plot(tnomol_value,b0nomol_mean,linewidth=0.5,marker='o',label='b0 - no molasses')
ax1.errorbar(tnomol_value,b0nomol_mean,yerr=b0nomol_errors,ls='none')

ax1.plot(tmol_value,b0mol_mean,linewidth=0.5,marker='o',label='b0 - molasses')
ax1.errorbar(tmol_value,b0mol_mean,yerr=b0mol_errors,ls='none')

plt.legend()

plt.savefig('/home/pedro/LAB/DATA/2014/Dec/MOT/comparison-b0-molasses.pdf')