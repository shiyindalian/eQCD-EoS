#!/usr/bin/env python
# -*- coding: utf-8 -*-
# sphinx_gallery_thumbnail_number = 3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib.ticker as ticker
import matplotlib as mpl

mpl.style.use('classic')


# Data for plotting
T=np.loadtxt('Deltals-data/muB0-ms0200/TMeV.dat')
T2p1=np.loadtxt('2p1data/mu0/TMEV.DAT')
Delta_ls_200=np.loadtxt('Deltals-data/muB0-ms0200/Delta_ls.dat')
Delta_ls_nor_200=Delta_ls_200/Delta_ls_200[1]

Delta_ls_120=np.loadtxt('Deltals-data/muB0-ms0120/Delta_ls.dat')
Delta_ls_nor_120=Delta_ls_120/Delta_ls_120[1]

Delta_ls_140=np.loadtxt('Deltals-data/muB0-ms0140/Delta_ls.dat')
Delta_ls_nor_140=Delta_ls_140/Delta_ls_140[1]

Delta_ls_150=np.loadtxt('Deltals-data/muB0-ms0150/Delta_ls.dat')
Delta_ls_nor_150=Delta_ls_150/Delta_ls_150[1]

Delta_ls_155=np.loadtxt('Deltals-data/muB0-ms0155/Delta_ls.dat')
Delta_ls_nor_155=Delta_ls_155/Delta_ls_155[1]

Delta_ls_160=np.loadtxt('Deltals-data/muB0-ms0160/Delta_ls.dat')
Delta_ls_nor_160=Delta_ls_160/Delta_ls_160[1]

Delta_ls_180=np.loadtxt('Deltals-data/muB0-ms0180/Delta_ls.dat')
Delta_ls_nor_180=Delta_ls_180/Delta_ls_180[1]

Delta_ls_2p1=np.loadtxt('2p1data/mu0/Ndelta_ls.dat')

# Wuppertal-Budapest
Deltal_ls_WB=np.loadtxt('DeltalR-WB/Deltals.dat')


# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)

line_Delta_ls,=ax1.plot(T,Delta_ls_nor_120,'b:',linewidth=1,markersize=5,label=r'$\bar m_s-\bar m_l=120\,\mathrm{MeV}$')
ax1.plot(T,Delta_ls_nor_150,'k-',linewidth=1,markersize=5,label=r'$\bar m_s-\bar m_l=150\,\mathrm{MeV}$')
ax1.plot(T2p1,Delta_ls_2p1,'c-',linewidth=1,markersize=5,label=r'2p1')
ax1.plot(T,Delta_ls_nor_155,'r--',linewidth=1,markersize=5,label=r'$\bar m_s-\bar m_l=155\,\mathrm{MeV}$')
ax1.plot(T,Delta_ls_nor_160,'g-.',linewidth=1,markersize=5,label=r'$\bar m_s-\bar m_l=160\,\mathrm{MeV}$')
#ax1.text(120, 0.3, r'$\mu_B=0$',fontsize=10, color='k')
ax1.text(230, 0.5, r'$\mu_B=0$',fontsize=10, color='k')
#ax1.axis([100,220,0,0.5])
ax1.axis([0,300,0,1])

points_WB=ax1.errorbar(Deltal_ls_WB[:,0],Deltal_ls_WB[:,1],yerr=Deltal_ls_WB[:,2],fmt='bs',ecolor='b',markersize=5,alpha=0.8,label=r'Lattice: WB Continuum')


ax1.set_xlabel('$T\,[\mathrm{MeV}]$', fontsize=14, color='black')
ax1.set_ylabel(r'$\Delta_{l,s}$', fontsize=14, color='black')

ax1.legend(loc=0,fontsize='xx-small',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)


#for ax in fig.axes:
#    ax.grid(True)

# Format the minor tick labels of the y-axis into empty strings with
# `NullFormatter`, to avoid cumbering the axis with too many labels.
#plt.gca().yaxis.set_minor_formatter(NullFormatter())
# Adjust the subplot layout, because the logit one may take more space
# than usual, due to y-tick labels like "1 - 10^{-3}"
fig.subplots_adjust(top=0.9, bottom=0.15, left=0.15, right=0.95, hspace=0.35,
                    wspace=0.35)


fig.savefig("Deltals_Range.pdf")

#plt.show()
