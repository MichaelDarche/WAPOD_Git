#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 13:44:10 2024

@author: michael
"""

import numpy as np
import sympy as sp

def propMatHomo(t,rho,cm,mat,modulation):
    if mat["Nat_-1"]=="solides modules":
        ampR=float(modulation["DeltaR"])
        ampC=float(modulation["DeltaC"])
        if modulation["Synchro"]=="Yes":
            FreqR=float(modulation["Freq"])
            FreqC=float(modulation["Freq"])
        else:
            FreqR=float(modulation["FreqR"])
            FreqC=float(modulation["FreqC"])
        if modulation["TypeMod"]=="Sinus":
            fctR=np.sin(2*np.pi*FreqR*t)
            fctC=np.sin(2*np.pi*FreqC*t)
        elif modulation["TypeMod"]=="Square":
            fctR=np.sign(np.sin(2*np.pi*FreqR*t))
            fctC=np.sign(np.sin(2*np.pi*FreqC*t))
        elif modulation["TypeMod"]=="Square NS":
            ratio=float(modulation["Ratio"])
            fctR=-np.sign(FreqR*t-int(FreqR*t)-ratio)
            fctC=-np.sign(FreqC*t-int(FreqC*t)-ratio)
        nameRhoi="Rho_"+str(1)
        nameCi="Cel_"+str(1)
        R=float(mat[nameRhoi])
        C=float(mat[nameCi])
        Cmod=(1+ampC*fctC)*C
        Rmod=(1+ampR*fctR)*R
    if mat["Nat_-1"]=="solides modules ER":
        ampR=float(modulation["DeltaR"])
        ampE=float(modulation["DeltaE"])
        if modulation["Synchro"]=="Yes":
            FreqR=float(modulation["Freq"])
            FreqE=float(modulation["Freq"])
        else:
            FreqR=float(modulation["FreqR"])
            FreqE=float(modulation["FreqE"])
        if modulation["TypeMod"]=="Sinus":
            fctR=np.sin(2*np.pi*FreqR*t)
            fctE=np.sin(2*np.pi*FreqE*t)
        elif modulation["TypeMod"]=="Square":
            fctR=np.sign(np.sin(2*np.pi*FreqR*t))
            fctE=np.sign(np.sin(2*np.pi*FreqE*t))
        elif modulation["TypeMod"]=="SquareSym":
            fctR=np.sign(np.cos(2*np.pi*FreqR*t))
            fctE=np.sign(np.cos(2*np.pi*FreqE*t))
        elif modulation["TypeMod"]=="Square NS":
            ratio=float(modulation["Ratio"])
            fctR=-np.sign(FreqR*t-int(FreqR*t)-ratio)
            fctE=-np.sign(FreqE*t-int(FreqE*t)-ratio)
        nameRhoi="Rho_"+str(0)
        nameCi="Cel_"+str(0)
        R=float(mat[nameRhoi])
        C=float(mat[nameCi])
        Rmod=(1+ampR*fctR)*R
        Emod=R*C*C/(1+ampE*fctE)
        Cmod=np.sqrt(Emod/Rmod)
        print(1/(1+ampR*fctR)/(1+ampE*fctE))
    return Rmod,Cmod

def propMatHomoForm(mat,modulation):
    t,x=sp.symbols("t x")
    if mat["Nat_-1"]=="solides modules CR":
        ampC=float(modulation["DeltaC"])
        ampR=float(modulation["DeltaR"])
        if modulation["Synchro"]=="Yes":
            FreqR=float(modulation["Freq"])
            FreqC=float(modulation["Freq"])
        else:
            FreqR=float(modulation["FreqR"])
            FreqC=float(modulation["FreqC"])
        if modulation["TypeMod"]=="Sinus":
            fctR=sp.sin(2*np.pi*FreqR*t)
            fctC=sp.sin(2*np.pi*FreqC*t)
        elif modulation["TypeMod"]=="Square":
            fctR=sp.sign(sp.sin((2*np.pi*FreqR*t)))
            fctC=sp.sign(0.001+sp.sin((2*np.pi*FreqC*t)))
        elif modulation["TypeMod"]=="Square NS":
            ratio=float(modulation["Ratio"])
            fctR=-sp.sign(FreqR*t-sp.floor(FreqR*t)-ratio)
            fctC=-sp.sign(0.001+FreqC*t-sp.floor(FreqC*t)-ratio)
        facC=(1+ampC*fctC)
        facR=(1+ampR*fctR)
    if mat["Nat_-1"]=="solides modules ER":
        ampE=float(modulation["DeltaE"])
        ampR=float(modulation["DeltaR"])
        if modulation["Synchro"]=="Yes":
            FreqR=float(modulation["Freq"])
            FreqE=float(modulation["Freq"])
        else:
            FreqR=float(modulation["FreqR"])
            FreqE=float(modulation["FreqE"])
        if modulation["TypeMod"]=="Sinus":
            fctR=sp.sin(2*np.pi*FreqR*t)
            fctE=sp.sin(2*np.pi*FreqE*t)
        elif modulation["TypeMod"]=="Cosinus":
            fctR=sp.cos(2*np.pi*FreqR*t)
            fctE=sp.cos(2*np.pi*FreqE*t)
        elif modulation["TypeMod"]=="Square":
            fctR=sp.sign(sp.sin((2*np.pi*FreqR*t)))
            fctE=sp.sign(0.0000001+sp.sin((2*np.pi*FreqE*t)))
        elif modulation["TypeMod"]=="SquareSym":
            fctR=sp.sign(sp.cos((2*np.pi*FreqR*t)))
            fctE=sp.sign(sp.cos((2*np.pi*FreqE*t)))
        elif modulation["TypeMod"]=="Square NS":
            ratio=float(modulation["Ratio"])
            fctR=-sp.sign(FreqR*t-sp.floor(FreqR*t)-ratio)
            fctE=-sp.sign(0.001+FreqE*t-sp.floor(FreqE*t)-ratio)
        facE=(1+ampE*fctE)
        facR=(1+ampR*fctR)
        facC=sp.sqrt(1/facE/facR)
    if mat["Nat_-1"]=="solides modules dephases":
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
            fctR=-sp.sign(FreqR*t-sp.floor(FreqR*t)-ratio)
            fctC=-sp.sign(FreqC*t-sp.floor(FreqC*t)-ratio)
        facC=(1+ampC*fctC)
        facR=(1+ampR*fctR)
    return facR,facC

def updateCFL(Nt,dt,t,rho,cm,mat,modulation):
    fCFL=1
    if mat["Nat_-1"]=="solides":
        rhoeff,Ceff=propMatHomo(0,rho[0],cm[0],mat,modulation)
        fCFL=Ceff/cm[0]
    for ts in range(Nt):
        if mat["Nat_-1"]=="solides modules":
            rhoeff,Ceff=propMatHomo(ts*dt,rho[0],cm[0],mat,modulation)
            fCFL=max(Ceff/cm[0], fCFL)
        if mat["Nat_-1"]=="solides modules ER":
            rhoeff,Ceff=propMatHomo(ts*dt,rho[0],cm[0],mat,modulation)
            fCFL=max(Ceff/cm[0], fCFL)
        if mat["Nat_-1"]=="masse-ressort module dephase":
            rhoeff,Ceff=propMatHomo(ts*dt,rho[0],cm[0],mat,modulation)
            fCFL=max(Ceff/cm[0], fCFL)
    return fCFL

def timeProp(Nt,dt,t,rho,cm,mat,modulation):
    rhoeffT=np.zeros(Nt)
    CeffT=np.zeros(Nt)
    if mat["Nat_-1"]=="solides":
        rhoeffT,CeffT=propMatHomo(0,rho[0],cm[0],mat,modulation)
    for ts in range(Nt):
        if mat["Nat_-1"]=="solides modules":
            rhoeffT[ts],CeffT[ts]=propMatHomo(ts*dt,rho[0],cm[0],mat,modulation)
        if mat["Nat_-1"]=="solides modules ER":
            rhoeffT[ts],CeffT[ts]=propMatHomo(ts*dt,rho[0],cm[0],mat,modulation)
        if mat["Nat_-1"]=="solides modules dephases":
            rhoeffT[ts],CeffT[ts]=propMatHomo(ts*dt,rho[0],cm[0],mat,modulation)
    return rhoeffT,CeffT