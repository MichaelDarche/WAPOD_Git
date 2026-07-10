#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 13:44:10 2024

@author: michael
"""

import numpy as np
import sympy as sp

def propMatHomo(t,x,rho,cm,mat,modulation):
    contacti="Contact_"+str(1)
    if mat[contacti]=="solide module":
        ampR=float(modulation["DeltaR"])
        ampC=float(modulation["DeltaC"])
        if modulation["Synchro"]=="Yes":
            FreqR=float(modulation["Freq"])
            FreqC=float(modulation["Freq"])
            celR=float(modulation["ModC"])
            celC=float(modulation["ModC"])
        else:
            FreqR=float(modulation["FreqR"])
            FreqC=float(modulation["FreqC"])
            celR=float(modulation["ModR"])
            celC=float(modulation["ModC"])
        if celR==0:
            phiRm=0
        else:
            phiRm=FreqR*2*np.pi/celR
        if celC==0:
            phiCm=0
        else:
            phiCm=FreqC*2*np.pi/celC
        if modulation["TypeMod"]=="Sinus":
            fctR=np.sin(2*np.pi*FreqR*t+phiRm*x)
            fctC=np.sin(2*np.pi*FreqC*t+phiCm*x)
        elif modulation["TypeMod"]=="Square":
            fctR=np.sign(np.sin(2*np.pi*FreqR*t+phiRm*x))
            fctC=np.sign(np.sin(2*np.pi*FreqC*t+phiCm*x))
        elif modulation["TypeMod"]=="Square NS":
            ratio=float(modulation["Ratio"])
            fctR=-np.sign(2*np.pi*FreqR*t-int(2*np.pi*FreqR*t)-ratio)
            fctC=-np.sign(2*np.pi*FreqC*t-int(2*np.pi*FreqC*t)-ratio)
        nameRhoi="Rho_"+str(1)
        nameCi="Cel_"+str(1)
        R=float(mat[nameRhoi])
        C=float(mat[nameCi])
        Cmod=(1+ampC*fctC)*C
        Rmod=(1+ampR*fctR)*R
    return Rmod,Cmod

def propMatHomoForm(mat,modulation):
    t,x=sp.symbols("t x")
    if mat["Nat"]=="solide module":
        ampC=float(modulation["DeltaC"])
        ampR=float(modulation["DeltaR"])
        if modulation["Synchro"]=="Yes":
            FreqR=float(modulation["Freq"])
            FreqC=float(modulation["Freq"])
            celR=float(modulation["ModC"])
            celC=float(modulation["ModC"])
        else:
            FreqR=float(modulation["FreqR"])
            FreqC=float(modulation["FreqC"])
            celR=float(modulation["ModR"])
            celC=float(modulation["ModC"])
        if celR==0:
            phiRm=0
        else:
            phiRm=FreqR*2*np.pi/celR
        if celC==0:
            phiCm=0
        else:
            phiCm=FreqC*2*np.pi/celC
        if modulation["TypeMod"]=="Sinus":
            fctR=sp.sin(2*np.pi*FreqR*t+phiRm*x)
            fctC=sp.sin(2*np.pi*FreqC*t+phiCm*x)
        elif modulation["TypeMod"]=="Square":
            fctR=sp.sign(sp.sin((2*np.pi*FreqR*t+phiRm*x)))
            fctC=sp.sign(sp.sin((2*np.pi*FreqC*t+phiCm*x)))
        elif modulation["TypeMod"]=="Square NS":
            ratio=float(modulation["Ratio"])
            fctR=-sp.sign(2*np.pi*FreqR*t-sp.floor(2*np.pi*FreqR*t)-ratio)
            fctC=-sp.sign(2*np.pi*FreqC*t-sp.floor(2*np.pi*FreqC*t)-ratio)
        facC=(1+ampC*fctC)
        facR=(1+ampR*fctR)
    return facR,facC