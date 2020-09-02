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
T=np.loadtxt('./fun_ck/TMEV.DAT')
chi2=np.loadtxt('./fun_ck/chi2.dat')
chi4=np.loadtxt('./fun_ck/chi4.dat')
R42=chi4/chi2
# Create figure
fig=plt.figure(figsize=(4.5, 3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)

ax1.plot(T,R42,'k-',linewidth=2,markersize=5)
ax1.axis([100,200,0.,1.2])

ax1.set_xlabel('$T [\mathrm{MeV}]$', fontsize=14, color='black')
ax1.set_ylabel('$R_{42}$', fontsize=14, color='black')

fig.subplots_adjust(top=0.9, bottom=0.17, left=0.17, right=0.95, hspace=0.35,
                    wspace=0.35)

fig.savefig("R42.pdf")
