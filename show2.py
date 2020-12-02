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
T=np.loadtxt('2p1225057/TMeV.dat')
chi2_250057=np.loadtxt('data/chi2_250.dat')
chi4_250057=np.loadtxt('data/chi4_250.dat')
chi6_250057=np.loadtxt('data/chi6_250.dat')
chi2_2p1=np.loadtxt('2p1225057/chi2.dat')
chi4_2p1=np.loadtxt('2p1225057/chi4.dat')
chi6_2p1=np.loadtxt('2p1225057/chi6.dat')
hotQCDR42=np.loadtxt('data/hotQCD_R42.dat')
WBT=np.loadtxt('../data1/WB_chix.dat')
WBR42=np.loadtxt('../data1/WB_R42.dat')
WBchi2=np.loadtxt('../data1/WB_chi2.dat')
WBchi6=np.loadtxt('../data1/WB_chi6.dat')
WBchi6err=np.loadtxt('../data1/WB_chi6erro.dat')
WBR62=WBchi6/WBchi2
hotQCDR62=np.loadtxt('data/hotQCD_R62.dat')
hotQCDb=np.loadtxt('data/hotQCDR62_b.dat')
hotQCDg=np.loadtxt('data/hotQCDR62_g.dat')

R42_2p1=chi4_2p1/chi2_2p1
R62_2p1=chi6_2p1/chi2_2p1
R42_250057=chi4_250057/chi2_250057
R62_250057=chi6_250057/chi2_250057
#T1=T/177
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)

#ax1.plot(T,R42_250057,'k-',linewidth=2,markersize=5,label=r'$FRG,Tc=250,\alpha=0.57$')
#ax1.plot(T,R42Z,'r-',linewidth=2,markersize=5,label=r'$FRG,Tc=250,\alpha=0.57$')
ax1.plot(T,R42_2p1,'k-',linewidth=2,markersize=5,label=r'$FRG2p1,Tc=225,\alpha=0.57$')
#ax1.plot(T,R42_2p1270,'k-',linewidth=2,markersize=5,label=r'$FRG2p1,Tc=270,\alpha=0.57$')
ax1.fill_between(hotQCDR42[:,0],hotQCDR42[:,1]+hotQCDR42[:,2],hotQCDR42[:,1]-hotQCDR42[:,2],alpha=0.25,facecolor='green',edgecolor='',label=r'HotQCD')
ax1.errorbar(WBR42[:,0],WBR42[:,1],yerr=WBR42[:,2],color='blue',marker='o',linestyle='',linewidth=2,markersize=5,fillstyle='none',alpha=1,label=r'Wuppertal-Budaspest')
ax1.axis([100,200,-0.1,1.2])

ax1.set_xlabel('$T [\mathrm{MeV}]$', fontsize=14, color='black')
ax1.set_ylabel('$\chi^B_4/\chi^B_2$', fontsize=14, color='black')

ax1.legend(loc=0,fontsize='x-small',frameon=True,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1,numpoints=1)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)

fig.subplots_adjust(top=0.9, bottom=0.17, left=0.15, right=0.95, hspace=0.35,
                    wspace=0.35)

fig.savefig("R42R62.pdf")
