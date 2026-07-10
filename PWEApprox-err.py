#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 10:37:30 2026

@author: michael
"""


import sys
sys.path.append("WAPOD/")

import numpy as np
import matplotlib.pyplot as plt

import funcPWEDissip

plt.rc('text',usetex=True)
plt.rc('font',family='serif')


tau=2*np.pi/20
valR=1200
c=2800
valE=valR*c**2
valGD=0
valGM=2400
deltaE=0.5
deltaR=0.0
deltaGD=0.0
deltaGM=0.0
fctModR,fctModE,fctModGD,fctModGM="np.cos","np.cos","np.cos","np.cos"
NE=12

Nu=np.linspace(0,1,21)
kgap1m=[]
kgap1p=[]
kgap2m=[]
kgap2p=[]
for inu in range(21):
    deltaE=Nu[inu]
    omega,VP=funcPWEDissip.fillMatrix(tau,valR,valE,0,0,deltaR,deltaE,deltaGD,deltaGM,fctModR,fctModE,fctModGD,fctModGM,N=16,NE=8)
    kgap1m.append(max(VP[0]))
    kgap1p.append(min(VP[1]))
    kgap2m.append(max(VP[1]))
    kgap2p.append(min(VP[2]))

def FuncK1m(x):
    y =2*np.pi/tau/(2*c*np.sqrt(1+x/2))
    return y

def FuncK1p(x):
    y =2*np.pi/tau/(2*c*np.sqrt(1-x/2))
    return y

def FuncK2m(x):
    alpha=6/x/x*(-1+np.sqrt(1+1/3*x**2))
    y =2*np.pi/tau/c*np.sqrt(alpha)
    return y

def FuncK2p(x):
    Delta = 6 / 5 / (x*x) * (1 - np.sqrt(1 - 5* x**2/ 3))
    y =2*np.pi/tau/c*np.sqrt(Delta)
    return y

def FuncK2mb(x):
    y=2*np.pi/tau/c*(1+x*x/24)
    return y

def FuncK2pb(x):
    y=2*np.pi/tau/c*(1+5*x*x/24)
    return y

def FuncK1mq(x):
    alpha=-(1+x/2)+np.sqrt((1+x/2)**2+x**2/8)
    y=2*2*np.pi/tau/c/x*np.sqrt(alpha)
    #y =2*np.pi/tau/(2*c)*(1-x/4)
    return y

def FuncK1pq(x):
    alpha=-(1-x/2)+np.sqrt((1-x/2)**2+x**2/8)
    y=2*2*np.pi/tau/c/x*np.sqrt(alpha)
    #y =2*np.pi/tau/(2*c)*(1+x/4)
    return y

def FuncK1mq2(x):
    alpha=-(1+x/2)+np.sqrt((1+x/2)**2+x**2/8)
    y=2*2*np.pi/tau/c/x*np.sqrt(alpha)
    y =2*np.pi/tau/(2*c)*(1-x/4)
    return y

def FuncK1pq2(x):
    alpha=-(1-x/2)+np.sqrt((1-x/2)**2+x**2/8)
    y=2*2*np.pi/tau/c/x*np.sqrt(alpha)
    y =2*np.pi/tau/(2*c)*(1+x/4)
    return y


branche1=FuncK1m(Nu)
branche2=FuncK1p(Nu)
branche1q=FuncK1mq(Nu)
branche2q=FuncK1pq(Nu)
branche1q2=FuncK1mq2(Nu)
branche2q2=FuncK1pq2(Nu)
branche3=FuncK2m(Nu)
branche3b=FuncK2mb(Nu)
branche4=FuncK2p(Nu)
branche4b=FuncK2pb(Nu)

plt.figure()
plt.plot(branche1, Nu, linewidth=2.5, color='k',label='Lin')
plt.plot(branche2, Nu, linewidth=2.5, color='k')
plt.plot(branche1q, Nu, linewidth=2.5, color='lightblue',label='Quad')
plt.plot(branche2q, Nu, linewidth=2.5, color='lightblue')
# plt.plot(branche1q2, Nu, linewidth=2.5, color='lightblue',label='Quad')
# plt.plot(branche2q2, Nu, linewidth=2.5, color='lightblue')
plt.plot(branche3,Nu, linewidth=2.5, color='k')
plt.plot(branche4, Nu, linewidth=2.5, color='k')
# plt.plot(branche3b, Nu,'--', linewidth=2.5, color='green')
plt.plot(branche4b, Nu,'--', linewidth=2.5, color='green',label="Lin Approx$^2$")
plt.plot(kgap1m,Nu,'*',linestyle="None", color='red',label='PWE')
plt.plot(kgap1p,Nu,'*',linestyle="None", color='red')
plt.plot(kgap2m,Nu,'*',linestyle="None", color='red')
plt.plot(kgap2p,Nu,'*',linestyle="None", color='red')
plt.xlabel(r'$k$',fontsize=16)
plt.xlim([0.002,0.012])
plt.ylabel(r'$\nu$',fontsize=16)
plt.ylim([0,1])
plt.legend(fontsize=16)
plt.tight_layout()

err1=(branche1-kgap1m)**2#/kgap1m**2
err2=(branche2-kgap1p)**2#/kgap1p**2
err1q=(branche1q-kgap1m)**2#/kgap1m**2
err2q=(branche2q-kgap1p)**2#/kgap1p**2
err3=(branche3-kgap2m)**2#/kgap2m**2
err4=(branche4-kgap2p)**2#/kgap2m**2
err3b=(branche3b-kgap2m)**2#/kgap2m**2
err4b=(branche4b-kgap2p)**2#/kgap2m**2



plt.figure()
# plt.loglog(Nu, Nu**3/10000000000,"--")
plt.loglog(Nu, err1,label="1$^{st}$ branch")
plt.loglog(Nu, err2,label="2$^{nd}$ branch")
plt.loglog(Nu, err1q,label="1$^{st}$ branch quad")
plt.loglog(Nu, err2q,label="2$^{nd}$ branch quad")
# plt.loglog(Nu, err3,label="3$^{rd}$ branch")
# # plt.loglog(Nu, err4,label="4$^{th}$ branch")
# plt.loglog(Nu, err3,label="3$^{rd}$ branch")
# plt.loglog(Nu, err4,label="4$^{rd}$ branch")
# plt.loglog(Nu, err4b,label="4$^{th}$ branch approx$^2$")
plt.loglog(Nu, Nu**3.87/550000000,":",label=r"$\propto\nu^{3.87}$")
plt.loglog(Nu, Nu**4.175/165000000,"--",label=r"$\propto\nu^{4.175}$")
plt.loglog(Nu, Nu**3.90/5700000000000,":",label=r"$\propto\nu^{3.9}$")
plt.loglog(Nu, Nu**4.75/900000000000,"--",label=r"$\propto\nu^{4.75}$")
plt.loglog(Nu, Nu**4.18/1900000000000,"-.",label=r"$\propto\nu^{4.18}$")
plt.loglog(Nu, Nu**5.78/50000000000,"-.",label=r"$\propto\nu^{5.8}$")
plt.xlabel(r'$\nu$',fontsize=16)
plt.xlabel(r'$\nu$',fontsize=16)
plt.xlim([0.048,1])
plt.ylabel(r'$\varepsilon_k$',fontsize=16)
plt.ylim([1.5e-18,4e-7])
plt.legend(fontsize=12)
plt.tight_layout()

plt.figure()
# plt.loglog(Nu, Nu**3/10000000000,"--")
# plt.loglog(Nu, err1,label="1$^{st}$ branch")
# plt.loglog(Nu, err2,label="2$^{nd}$ branch")
# plt.loglog(Nu, err1q,label="1$^{st}$ branch quad")
# plt.loglog(Nu, err2q,label="2$^{nd}$ branch quad")
# plt.loglog(Nu, err3,label="3$^{rd}$ branch")
# # plt.loglog(Nu, err4,label="4$^{th}$ branch")
plt.loglog(Nu, err3,label="3$^{rd}$ branch")
plt.loglog(Nu, err4,label="4$^{rd}$ branch")
plt.loglog(Nu, err4b,label="4$^{th}$ branch approx$^2$")
plt.loglog(Nu, Nu**8.175/950000,"--",label=r"$\propto\nu^{8.175}$")
plt.loglog(Nu, Nu**8.08/11000000,":",label=r"$\propto\nu^{8.08}$")
plt.loglog(Nu, Nu**7.30/80000000000,"-.",label=r"$\propto\nu^{7.3}$")
plt.loglog(Nu, Nu**5.30/5000000000000,":",label=r"$\propto\nu^{5.3}$")
plt.xlabel(r'$\nu$',fontsize=16)
plt.xlim([0.048,1])
plt.ylabel(r'$\varepsilon_k$',fontsize=16)
plt.ylim([1.5e-21,1e-6])
plt.legend(fontsize=12)
plt.tight_layout()