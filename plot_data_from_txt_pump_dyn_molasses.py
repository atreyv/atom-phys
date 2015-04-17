# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 12:06:56 2014

@author: pedro
"""

from libphys import *
mpl.rcParams['figure.figsize'] = (16.0, 8.0)
mpl.rc('font', size=14)

dname1 = '/home/pedro/LAB/DATA/2014/Dec/PF/05_12_14/pf03/cam113/'
dname2 = '/home/pedro/LAB/DATA/2014/Dec/PF/05_12_14/pf04/cam113/'
dname3 = '/home/pedro/LAB/DATA/2014/Dec/PF/05_12_14/pf01/cam113/'
f1 = 'pattern-parallel_pump_dyn-circpol-no_molasses.dat'
f2 = 'pattern-parallel_pump_dyn-circpol-molasses.dat'
f3 = 'pattern-parallel_pump_dyn-linpol_no-molasses.dat'
nomol = np.loadtxt(dname1+f1, skiprows=1)
mol = np.loadtxt(dname2+f2, skiprows=1)
lin = np.loadtxt(dname3+f3, skiprows=1)
nomol = np.transpose(nomol)
mol = np.transpose(mol)
lin = np.transpose(lin)

fig = plt.figure()
#fig.subplots_adjust(left=0.1,right=0.9)
ax1 = fig.add_subplot(111)
ax1.set_ylabel('a.u.')

#ax2 = fig.add_subplot(212)
ax1.set_xlabel(r't ($\mathrm{\mu s}$)')


ax1.plot(nomol[0],nomol[1],linewidth=0.5,marker='o',label='circ pol - no molasses')
ax1.errorbar(nomol[0],nomol[1],yerr=nomol[2],ls='none',color='b')

ax1.plot(mol[0],mol[1],linewidth=0.5,marker='x',label='circ pol - molasses')
ax1.errorbar(mol[0],mol[1],yerr=mol[2],ls='none',color='g')

ax1.plot(lin[0],lin[1],linewidth=0.5,label='lin pol - no molasses')
ax1.errorbar(lin[0],lin[1],yerr=lin[2],ls='none',color='r')

plt.legend()

#fig.savefig('/home/pedro/LAB/DATA/2014/Dec/PF/05_12_14/'+'comparison-molasses-circ_pol.pdf')