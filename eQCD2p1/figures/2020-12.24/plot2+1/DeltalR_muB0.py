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
T=np.loadtxt('muB-Tto300/mu0/TMeV.dat')
T_2p1=np.loadtxt('2p1data/mu0/TMEV.DAT')

fpi0=np.loadtxt('muB-Tto300/mu0/fpi.dat')
c0=np.loadtxt('muB-Tto300/mu0/c.dat')

fpi0_2p1=np.loadtxt('2p1data/mu0/FPI.DAT')
c0_2p1=np.loadtxt('2p1data/mu0/CLCS.DAT')

DeltalR_rescale=0.0414192*0.975
D2p1_rescale=0.0414192*0.80

DeltalR0_1=-c0*fpi0
DeltalR0_2=DeltalR0_1-DeltalR0_1[0]
DeltalR0=DeltalR0_2*DeltalR_rescale

D2p1_1=-c0_2p1[:,0]*fpi0_2p1[:,0]
D2p1_2=D2p1_1-D2p1_1[0]
Del_2p1=D2p1_2*D2p1_rescale

# Wuppertal-Budapest
DeltalR_WB=np.loadtxt('DeltalR-WB/DeltalR.dat')


# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)

line_DeltalR0,=ax1.plot(T,DeltalR0,'k-',linewidth=2,markersize=5,label=r'FRG $N_f=2$')
line_DeltalR0,=ax1.plot(T_2p1,Del_2p1,'c-',linewidth=2,markersize=5,label=r'FRG $N_f=2+1$')
#ax1.text(120, 0.3, r'$\mu_B=0$',fontsize=10, color='k')
ax1.text(50, 0.3, r'$\mu_B=0$',fontsize=10, color='k')
#ax1.axis([100,220,0,0.5])
ax1.axis([0,300,0,0.5])

points_WB=ax1.errorbar(DeltalR_WB[:,0],DeltalR_WB[:,1],yerr=DeltalR_WB[:,2],fmt='bs',ecolor='b',markersize=5,alpha=0.8,label=r'Lattice: WB Continuum')


ax1.set_xlabel('$T\,[\mathrm{MeV}]$', fontsize=14, color='black')
ax1.set_ylabel(r'$\Delta_{l,R}$', fontsize=14, color='black')

ax1.legend(loc=0,fontsize='x-small',frameon=False,shadow=False,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1)

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


fig.savefig("DeltalR_muB0.pdf")

#plt.show()
