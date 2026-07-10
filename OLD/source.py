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
    fct=(1-2*(np.pi*f*t)**2)*np.exp(-(np.pi*f*t)**2)
    return fct   
#% infinite harmonic
def harmonic_t(t,f):
    omega=2*np.pi*f
    fct=np.sin(omega*t)
    return fct
#% finite harmonic
def harmonic_n(t,f):
    omega=2*np.pi*f
    n=input("number of periods")
    if t>=0 and t<=1/f*n :
        fct=np.sin(omega*t)#A0*np.sin(2**0*omega*t)+A1*np.sin(2**1*omega*t)+A2*np.sin(2**2*omega*t)+A3*np.sin(2**3*omega*t)
    else:
        fct=0
    return fct
    
    
temporalFct= {"JKPS C6":jkps_t,
              "Ricker":ricker_t,
              "SINUS":harmonic_t,
              "N-Harmonic":harmonic_n}

def choice_timefct(sourceparameters, t):
    f=float(sourceparameters["Frequence"])
    amp=float(sourceparameters["Force"])
    if sourceparameters["Temporel"] in temporalFct:
        g = temporalFct[sourceparameters["Temporel"]]
        fct = amp*g(t,f)
        return fct
    else:
        return "Source function not defined"


    


        
        

        
        
