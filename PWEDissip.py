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

tau=2*np.pi/20
valR=1200
c=2800
valE=valR*c**2
valGD=0
valGM=0#2400
deltaE=0.25
deltaR=0.0
deltaGD=0.0
deltaGM=0.0
fctModR,fctModE,fctModGD,fctModGM="np.cos","np.cos","np.cos","np.cos"
NE=12


Omega=2*np.pi/tau

mu=valGM/valR/Omega

delta1=1/(4-2*deltaE*np.sqrt(1-mu**2))#     1/(4*deltaE+np.sqrt(1-mu**2))
delta2=1/(4+2*deltaE*np.sqrt(1-mu**2))
delta3=2/(2-deltaE*np.sqrt(1-mu**2))#     1/(4*deltaE+np.sqrt(1-mu**2))
delta4=2/(2+deltaE*np.sqrt(1-mu**2))

k1=Omega/c*np.sqrt(delta1)
k2=Omega/c*np.sqrt(delta2)
k3=Omega/c*np.sqrt(delta3)
k4=Omega/c*np.sqrt(delta4)



# deltan05=[0.223, 0.304]
# # deltan05=[0.199, 0.328, 0.980, 1.112, 2.292, 2.388]
# # Omegan=2*np.pi/tau
# knth05=Omegan/c*np.sqrt(deltan05)

# gamma05=valGM/Omegan/valR


omega,VP=funcPWEDissip.fillMatrix(tau,valR,valE,valGD,valGM,deltaR,deltaE,deltaGD,deltaGM,fctModR,fctModE,fctModGD,fctModGM,N=16)


plt.figure()
#    ...: #plt.pcolor(k_map,Om_plot,np.log10(np.abs(DD_micro)))
for i in range(NE):
    plt.plot((VP[i,:]),omega,'red',linewidth=3)
# for i in range(knth05.shape[0]):
#     plt.axvline(knth05[i])
plt.axvline(k1)
plt.axvline(k2)
plt.axvline(k3)
plt.axvline(k4)
plt.xlabel(r'$k$ (rad/m)',fontsize=16)
plt.ylabel(r'$\omega$ (rad/s)',fontsize=16)
# plt.legend(fontsize=16)
plt.title(r'$N=16$',fontsize=16)
plt.tight_layout()

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

