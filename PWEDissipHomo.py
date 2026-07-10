#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 11:19:58 2026

@author: michael
"""


import sys
sys.path.append("WAPOD/")

import numpy as np
import matplotlib.pyplot as plt

import funcPWEDissip

tau=1/30
R=1200
c=2800
E=R*c**2
#############
M0=20000
delM=-0.9
C0=1e-9/2.45
delC=0.9
h=10
GammaC=0e-8
GammaM=5e4
##%%%%%
valR=R+M0/h
valE=1/(1/E+C0/h)
deltaE=delC*C0*E/(h+C0*E)
deltaR=delM*M0/(h*R+M0)
#############
valGD=GammaC/h
valGM=GammaM/h

deltaGD=0.0
deltaGM=0.0
fctModR,fctModE,fctModGD,fctModGM="np.sin","np.sin","np.sin","np.sin"
NE=8

NG=49
kgap1m=[]
kgap1p=[]
kgap2m=[]
kgap2p=[]
kgap3m=[]
kgap3p=[]
kgap4m=[]
kgap4p=[]
valGammaC=np.linspace(1e-8,1e-7,NG)
valGammaM=np.linspace(1e4,1e5,NG)
for iG in range(NG):
    valGD=valGammaC[iG]/h
    valGM=valGammaM[iG]/h
    omega,VP=funcPWEDissip.fillMatrix(tau,valR,valE,valGD,valGM,deltaR,deltaE,deltaGD,deltaGM,fctModR,fctModE,fctModGD,fctModGM,N=10,NE=8)
    kgap1m.append(max(VP[0]))
    kgap1p.append(min(VP[1]))
    kgap2m.append(max(VP[1]))
    kgap2p.append(min(VP[2]))
    kgap3m.append(max(VP[2]))
    kgap3p.append(min(VP[3]))
    kgap4m.append(max(VP[3]))
    kgap4p.append(min(VP[4]))
    for i in range(NE):
        plt.plot((VP[i,:])/2/np.pi,omega/2/np.pi,linewidth=3)
plt.xlabel(r'$k$ (rad/m)',fontsize=16)
plt.ylabel(r'$\omega$ (rad/s)',fontsize=16)
# plt.legend(fontsize=16)
plt.title(r'$N=16$',fontsize=16)
plt.tight_layout()

plt.figure()
plt.plot(valGammaM/valR/(2*np.pi/tau),kgap1m,label=r"$k^+_1$")
plt.plot(valGammaM/valR/(2*np.pi/tau),kgap1p,label=r"$k^-_2$")
plt.plot(valGammaM/valR/(2*np.pi/tau),kgap2m,label=r"$k^+_2$")
plt.plot(valGammaM/valR/(2*np.pi/tau),kgap2p,label=r"$k^-_3$")
plt.plot(valGammaM/valR/(2*np.pi/tau),kgap3m,label=r"$k^+_3$")
plt.plot(valGammaM/valR/(2*np.pi/tau),kgap3p,label=r"$k^-_4$")
plt.plot(valGammaM/valR/(2*np.pi/tau),kgap4m,label=r"$k^+_4$")
plt.plot(valGammaM/valR/(2*np.pi/tau),kgap4p,label=r"$k^-_5$")
# plt.axvline()
plt.xlabel(r"$\gamma_M$",fontsize=16)
plt.ylabel(r"$k$",fontsize=16)
plt.legend(fontsize=14)
plt.tight_layout()

plt.figure()
for i in range(NE):
    plt.plot((VP[i,:])/2/np.pi,omega/2/np.pi,'red',linewidth=3)
plt.xlabel(r'$k$ (rad/m)',fontsize=16)
plt.ylabel(r'$\omega$ (rad/s)',fontsize=16)
# plt.legend(fontsize=16)
plt.title(r'$N=16$',fontsize=16)
plt.tight_layout()
#    ...: #plt.pcolor(k_map,Om_plot,np.log10(np.abs(DD_micro)))

# deltan09=[0.182065, 0.408018, 0.967324, 1.387383, 2.360391, 2.920547]
# Omegan=2*np.pi/tau
# knth09=Omegan/c*np.sqrt(deltan09)
# knth09=np.array([0.003048,0.004563,0.007025,0.008413,0.010974,0.012207,0.014904,0.015978,0.018826,0.019737,0.022746,0.023489 ])
# gamma09=valGM/Omegan/valR

# deltaE=0.9
# omega,VP=funcPWEDissip.fillMatrix(tau,valR,valE,valGD,valGM,deltaR,deltaE,deltaGD,deltaGM,fctModR,fctModE,fctModGD,fctModGM,N=32)

# plt.figure()
# #    ...: #plt.pcolor(k_map,Om_plot,np.log10(np.abs(DD_micro)))
# for i in range(NE):
#     plt.plot(np.real(VP[i,:]),omega,'red',linewidth=3)
#     plt.plot(np.imag(VP[i,:]),omega,'green',linewidth=3)
# for i in range(knth09.shape[0]):
#     plt.axvline(knth09[i])
# plt.xlabel(r'$k$ (rad/m)',fontsize=16)
# plt.ylabel(r'$\omega$ (rad/s)',fontsize=16)
# # plt.legend(fontsize=16)
# plt.title(r'$N=32$',fontsize=16)
# plt.tight_layout()

# # deltan07=[0.182065, 0.408018, 0.967324, 1.387383, 2.360391, 2.920547]
# # Omegan=2*np.pi/tau
# # knth09=Omegan/c*np.sqrt(deltan09)

# # gamma07=valGM/Omegan/valR

# # deltaE=0.7
# # omega,VP=funcPWEDissip.fillMatrix(tau,valR,valE,valGD,valGM,deltaR,deltaE,deltaGD,deltaGM,fctModR,fctModE,fctModGD,fctModGM)

# # plt.figure()
# # #    ...: #plt.pcolor(k_map,Om_plot,np.log10(np.abs(DD_micro)))
# # for i in range(NE):
# #     plt.plot((VP[i,:]),omega,'red',linewidth=3)
# # for i in range(knth09.shape[0]):
# #     plt.axvline(knth09[i])
# # plt.xlabel(r'$k$ (rad/m)',fontsize=16)
# # plt.ylabel(r'$\omega$ (rad/s)',fontsize=16)
# # # plt.legend(fontsize=16)
# # plt.title(r'$N=16$',fontsize=16)
# # plt.tight_layout()

