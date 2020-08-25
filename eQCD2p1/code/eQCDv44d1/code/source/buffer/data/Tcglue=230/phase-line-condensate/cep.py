#!/usr/bin/env python
# -*- coding: utf-8 -*-
# sphinx_gallery_thumbnail_number = 3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib.ticker as ticker


# Data for plotting

#dmfdT_qti=np.loadtxt('muB/muB0/dmfdT.dat')
#TMeV=np.loadtxt('muB/muB0/TMeV.dat')
list1 =  np.zeros((81,3))

for ii in range(81):
    fns = str(ii)
    filename_fpi = '../data-T300/mu'+fns+'/fpi.dat'
    filename_c = '../data-T300/mu'+fns+'/c.dat'
    filename_TMeV = '../data-T300/mu'+fns+'/TMeV.dat'

    fpi= np.loadtxt(filename_fpi)
    c= np.loadtxt(filename_c)
    Delta=fpi*c
#mf=np.loadtxt('muB/muB0/mf.dat')

    T=np.loadtxt(filename_TMeV)

#muB=np.loadtxt('CEP/muB.dat')
#T_low=np.loadtxt('CEP/Tlow.dat')
#T_high=np.loadtxt('CEP/Thigh.dat')

#Tc=np.loadtxt('CEP/Tc.dat')
#Tc1=np.loadtxt('CEP/Tc1.dat')


    y_data=Delta
    n_data=len(y_data)

    dydx_data=range(n_data)
    for i in range(n_data):
        if i==0:
            dydx_data[0]=(y_data[1]-y_data[0])
        elif i==n_data-1:
            dydx_data[i]=(y_data[i]-y_data[i-1])
        else:
            dydx_data[i]=(y_data[i+1]-y_data[i-1])/2.
    dydx_data=dydx_data/np.amin(dydx_data)

    y_cut=0.8

    for i in range(n_data-1):
        if dydx_data[i]<=y_cut and dydx_data[i+1]>y_cut :
            T_cut_low=(y_cut-dydx_data[i])/(dydx_data[i+1]-dydx_data[i])+T[i]
        elif dydx_data[i]>y_cut and dydx_data[i+1]<=y_cut :
            T_cut_high=(y_cut-dydx_data[i])/(dydx_data[i+1]-dydx_data[i])+T[i]



    list1[ii,0] = ii * 10
    list1[ii,1] = T_cut_low
    list1[ii,2] = T_cut_high
#print('T_cut_low=',T_cut_low)
#print('T_cut_high=',T_cut_high)
np.savetxt('./pb-Nf2p1.dat',list1)

dmfdT=dydx_data


