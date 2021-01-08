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
T=np.loadtxt('TMEV.DAT')
MF=np.loadtxt('MF.DAT')
#MF2=np.loadtxt('MF2.DAT')
MBO=np.loadtxt('MBO_PHY.DAT')

# Create figure
fig=plt.figure(figsize=(9., 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(121)

ax1.plot(T,MF[:,0],'k-',linewidth=2,markersize=5,label=r'$M_l$')
ax1.plot(T,MF[:,1],'r-',linewidth=2,markersize=5,label=r'$M_s$')
#ax1.plot(T,MF2,'g.-',linewidth=2,markersize=5,label=r'$M_s$')
ax1.axis([0,250,0.,550])

ax1.set_xlabel('$T [\mathrm{MeV}]$', fontsize=14, color='black')
ax1.set_ylabel('$M_f [\mathrm{MeV}]$', fontsize=14, color='black')

ax1.legend(loc=0,fontsize='x-small',frameon=False,shadow=False,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

ax2=fig.add_subplot(122)

f0,=ax2.plot(T,MBO[:,0],'r:',linewidth=2,markersize=5,label=r'$f_0$')
a0,=ax2.plot(T,MBO[:,1],'b:',linewidth=2,markersize=5,label=r'$a_0$')
kappa,=ax2.plot(T,MBO[:,2],'g:',linewidth=2,markersize=5,label=r'$\kappa$')
sigma,=ax2.plot(T,MBO[:,3],'c:',linewidth=2,markersize=5,label=r'$\sigma$')
eta,=ax2.plot(T,MBO[:,4],'m-',linewidth=2,markersize=5,label='$\eta$\'')
pion,=ax2.plot(T,MBO[:,5],'k-',linewidth=2,markersize=5,label=r'$\pi$')
kaon,=ax2.plot(T,MBO[:,6],'r-',linewidth=2,markersize=5,label=r'$K$')
etap,=ax2.plot(T,MBO[:,7],'y-',linewidth=2,markersize=5,label=r'$\eta$')
ax2.axis([0,250,0,1800])

ax2.set_xlabel('$T [\mathrm{MeV}]$', fontsize=14, color='black')
ax2.set_ylabel('$M_{meson} [\mathrm{MeV}]$', fontsize=14, color='black')


legend1=ax2.legend(handles=[eta,pion,kaon,etap],loc=2,fontsize='x-small',frameon=False,shadow=False,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1)
legend0=ax2.legend(handles=[f0,a0,kappa,sigma],loc=9,fontsize='x-small',frameon=False,shadow=False,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1)

plt.gca().add_artist(legend1)

for label in ax2.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax2.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.17, left=0.10, right=0.95, hspace=0.35,
                    wspace=0.35)

fig.savefig("MFMBO.pdf")
