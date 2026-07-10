#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 09:24:58 2023

@author: michael
"""
import numpy as np
################## Source function #################     
#% JKPS C6 spatial   
# def jkps_x(x,xm,xM,lo):
#     A0=1
#     A1=-21/32
#     A2=63/768
#     A3=-1/512
#     k=2*np.pi/lo
#     if x>=xm and x<=xM:
#         fct=A0*np.sin(2**0*k*x)+A1*np.sin(2**1*k*x)+A2*np.sin(2**2*k*x)+A3*np.sin(2**3*k*x)
#     else:
#         fct=0
#     return fct
#% JKPS C6 temporal
def jkps_t(t,f):
    A0=1
    A1=-21/32
    A2=63/768
    A3=-1/512
    omega=2*np.pi*f
    if t>=0 and t<=1/f :
        fct=A0*np.sin(2**0*omega*t)+A1*np.sin(2**1*omega*t)+A2*np.sin(2**2*omega*t)+A3*np.sin(2**3*omega*t)
    else:
        fct=0
    return fct
#% Ricker spatial   
# def ricker_x(x,xm,xM,lo):
#     A0=1
#     A1=-21/32
#     A2=63/768
#     A3=-1/512
#     k=2*np.pi/lo
#     if x>=xm and x<=xM:
#         fct=A0*np.sin(2**0*k*x)+A1*np.sin(2**1*k*x)+A2*np.sin(2**2*k*x)+A3*np.sin(2**3*k*x)
#     else:
#         fct=0
#     return fct
#% Ricker temporal
def ricker_t(t,f):
    tau=1.5/f
    T=t-tau
    fct=(1-2*(np.pi*f*T)**2)*np.exp(-(np.pi*f*T)**2)
    return fct   
#% infinite harmonic
def harmonic_t(t,f):
    omega=2*np.pi*f
    if t<=0:
        fct=0
    else:
        fct=np.sin(omega*t)
    
    return fct
#% finite harmonic
def harmonic_n(t,f,n):
    omega=2*np.pi*f
    if t>=0 and t<=1/f*n :
        fct=np.sin(omega*t)#A0*np.sin(2**0*omega*t)+A1*np.sin(2**1*omega*t)+A2*np.sin(2**2*omega*t)+A3*np.sin(2**3*omega*t)
    else:
        fct=0
    return fct

def expplus(t,f):
    omega=2*np.pi*f
    fct=np.exp(-1j*omega*t)#A0*np.sin(2**0*omega*t)+A1*np.sin(2**1*omega*t)+A2*np.sin(2**2*omega*t)+A3*np.sin(2**3*omega*t)
    return fct

def fileSource(file):
    sce=np.readtxt(file)
    
    return 0
    
temporalFct= {"JKPS C6":jkps_t,
              "Ricker":ricker_t,
              "SINUS":harmonic_t,
              "N-Harmonic":harmonic_n,
              "Exp+":expplus,
              "File":fileSource}

def choice_timefct(sourceparameters, t):
    f=float(sourceparameters["Frequence"])
    amp=float(sourceparameters["Force"])
    if sourceparameters["Temporel"] in temporalFct:
        g = temporalFct[sourceparameters["Temporel"]]
        if sourceparameters["Temporel"]=="N-Harmonic":
            n=float(sourceparameters["N-Arches"])
            fct = amp*g(t,f,n)
        elif sourceparameters["Temporel"]=="File":
            fct=0
        else:
            fct = amp*g(t,f)
        return fct
    else:
        return "Source function not defined"

# import matplotlib.pyplot as plt
def TFsce(sce,f,Nf,fc,Kf,Tmax):
    df=Kf*fc/Nf 
    t=np.linspace(0,Tmax,Nf)
    dt=t[1]
    TF=0
    # fsce=np.zeros(Nf+1)
    # for k in range(Nf):
    #     fsce[k]=sce(t[k],fc)
    # plt.plot(t,fsce)
    for k in range(int(t.shape[0])):
        if k==0 or k==int(t.shape[0]):
            TF=TF+1/2*sce(t[k],fc)*np.exp(-2j*np.pi*f*df*t[k])*dt
        else:
            TF=TF+sce(t[k],fc)*np.exp(-2j*np.pi*f*df*t[k])*dt
    return TF


# Tmax=0.06
# NFourier=1024
# Kf=10
# fc=50
# TF=np.zeros(NFourier)
# fftF=np.zeros(NFourier)
# freq=np.zeros(NFourier)
# df=Kf*fc/NFourier
# sce=jkps_t

# for f in range(NFourier):
#     freq[f]=f*df
#     TF[f]=TFsce(sce,f,NFourier,fc,Kf,Tmax)

# Ak=np.zeros(4)
# Ak[0]= 1.0
# Ak[1]= - 21 / 32
# Ak[2]=  63 / 768
# Ak[3]= - 1 / 512      
# Eps = 1e-7 
# FF0=np.zeros(NFourier)
# OmegaC= 2 * np.pi * fc
# for f in range(NFourier):
#     Omega =2 * np.pi * f*df
#     for K in range(4):
#         BetaK = 2**(K)
#         if abs(Omega - BetaK * OmegaC) < Eps * BetaK * OmegaC:
#             FF0[f]=FF0[f]-1j* Ak[K] * np.pi / (OmegaC)
#         else:
#             Xint= 2 * np.pi * Omega / OmegaC
#             Denom=(Omega)**2 - (BetaK * OmegaC)**2
#             Facteur= Ak[K] * BetaK * OmegaC
#             FF0[f]=FF0[f] + Facteur * ((np.cos(Xint) - 1)-1j*np.sin(Xint)) / Denom
# # freq2=np.fft.fftfreq(Nf,d=t[1])
# # Ak=np.zeros(4)
# # Ak[0]= 1.0
# # Ak[1]= - 21 / 32
# # Ak[2]=  63 / 768
# # Ak[3]= - 1 / 512      
# # Eps = 1e-7 
# # FF0=np.zeros(Nf)
# # Tampon= 0.0
# # OmegaC= 2 * np.pi * fc
# # for f in range(Nf):
# #     Omega =2 * np.pi * f*df
# #     for K in range(4):
# #         BetaK = 2**(K-1)
# #         if abs(Omega - BetaK * OmegaC) < Eps * BetaK * OmegaC:
# #             FF0[f]=FF0[f]+1j* Ak[K] * np.pi / (OmegaC)
# #         else:
# #             X= 2 * np.pi * Omega / OmegaC
# #             Denom=np.sqrt(Omega) - np.sqrt(BetaK * OmegaC)
# #             Facteur= Ak[K] * BetaK * OmegaC
# #             FF0[f]=FF0[f] + Facteur * ((np.cos(X) - 1)-1j*np.sin(X)) / Denom

# plt.figure()
# plt.plot(freq,TF)
# plt.plot(freq,FF0)  
        
