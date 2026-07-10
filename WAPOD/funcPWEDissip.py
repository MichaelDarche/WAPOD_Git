#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 17:48:31 2025

@author: michael
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as spl
import scipy.sparse.linalg as ssl

def fctheavyside(t,tau,phi):
    Tphi=tau*phi
    res=t.copy()
    for i in range(t.shape[0]):
        if t[i]<Tphi:
            res[i]=+1
        else:
            res[i]=-1
    return res    


def rhot(t,tau,rho,Delta,fctM):
    if fctM=="np.sin":
        rhot=rho*(1+Delta*np.sin((2*np.pi*t/tau)))
    if fctM=="np.cos":
        rhot=rho*(1+Delta*np.cos((2*np.pi*t/tau)))
    elif fctM=="np.sign":
        rhot0=np.sin(2*np.pi*t/tau)
        rhot=rho*(1+Delta*np.sign(rhot0))
    elif fctM=="heavy":
        rhot=rho*(1+Delta*fctheavyside(t,tau))
    return rhot

def Et(t,tau,E,delta,fctM):
    if fctM=="np.sin":
        Et=(E/(1+delta*np.sin((2*np.pi*t/tau))))
    if fctM=="np.cos":
        Et=E*(1+delta*np.cos((2*np.pi*t/tau)))
    elif fctM=="np.sign":
        e0=np.sin((2*np.pi*t/tau))
        Et=(E/(1+delta*np.sign(e0)))
    elif fctM=="heavy":
        Et=(E/(1+delta*fctheavyside(t,tau)))
    return Et

def GammaMt(t,tau,GM,DGM,fctGM):
    if fctGM=="np.sin":
        GMt=GM*(1+DGM*np.sin((2*np.pi*t/tau)))
    if fctGM=="np.cos":
        GMt=GM*(1+DGM*np.cos((2*np.pi*t/tau)))
    elif fctGM=="np.sign":
        rhot0=np.sin(2*np.pi*t/tau)
        GMt=GM*(1+DGM*np.sign(rhot0))
    elif fctGM=="heavy":
        GMt=GM*(1+DGM*fctheavyside(t,tau))
    return GMt

def GammaDt(t,tau,GD,dGD,fctGD):
    if fctGD=="np.sin":
        GDt=GD*(1+dGD*np.sin((2*np.pi*t/tau)))
    if fctGD=="np.cos":
        GDt=GD*(1+dGD*np.cos((2*np.pi*t/tau)))
    elif fctGD=="np.sign":
        e0=np.sin((2*np.pi*t/tau))
        GDt=GD*(1+dGD*np.sign(e0))
    elif fctGD=="heavy":
        GDt=GD*(1+dGD*fctheavyside(t,tau))
    return GDt

def rhotD1(t,tau,rho,Delta,fctM):
    if fctM=="np.sin":
        rhotD=rho*(2*np.pi/tau*Delta*np.cos((2*np.pi*t/tau)))
    if fctM=="np.cos":
        rhotD=rho*(-1*2*np.pi/tau*Delta*np.sin((2*np.pi*t/tau)))
    elif fctM=="np.sign":
        rhotD=0
    elif fctM=="heavy":
        rhotD=0
    return rhotD

def rhotD2(t,tau,rho,Delta,fctM):
    if fctM=="np.sin":
        rhotD2=rho*(-2*np.pi/tau*Delta*2*np.pi/tau*np.sin((2*np.pi*t/tau)))
    if fctM=="np.cos":
        rhotD2=rho*(-1*2*np.pi/tau*2*np.pi/tau*Delta*np.cos((2*np.pi*t/tau)))
    elif fctM=="np.sign":
        rhotD2=0
    elif fctM=="heavy":
        rhotD2=0
    return rhotD2

def DtD(t,tau,E,delta,fctM):
    if fctM=="np.sin":
        DtD=1/E*(2*np.pi/tau*delta*np.cos((2*np.pi*t/tau)))
    if fctM=="np.cos":
        DtD=-1/E*(-delta*2*np.pi/tau*np.sin(2*np.pi*t/tau))/(1+delta*np.cos(2*np.pi*t/tau))**2
        #1/E*(-2*np.pi/tau*delta*np.sin((2*np.pi*t/tau)))
    elif fctM=="np.sign":
        DtD=0
    elif fctM=="heavy":
        DtD=0
    return DtD

def GammaMtD(t,tau,GM,DGM,fctGM):
    if fctGM=="np.sin":
        GMtD=GM*(2*np.pi/tau*DGM*np.cos((2*np.pi*t/tau)))
    if fctGM=="np.cos":
        GMtD=GM*(-2*np.pi/tau*DGM*np.sin((2*np.pi*t/tau)))
    elif fctGM=="np.sign":
        GMtD=0
    elif fctGM=="heavy":
        GMtD=0
    return GMtD

def GammaDtD(t,tau,GD,dGD,fctGD):
    if fctGD=="np.sin":
        GDtD=GD*(2*np.pi/tau*dGD*np.cos((2*np.pi*t/tau)))
    if fctGD=="np.cos":
        GDtD=GD*(-2*np.pi/tau*dGD*np.sin((2*np.pi*t/tau)))
    elif fctGD=="np.sign":
        GDtD=0
    elif fctGD=="heavy":
        GDtD=0
    return GDtD

def xit(t,tau,rho,E,GD,GM,Delta,delta,dGD,DGM,fctR,fctE,fctGD,fctGM):
    xi1=rhotD1(t,tau,rho,Delta,fctR)
    xi2=rhot(t,tau,rho,Delta,fctR)*Et(t,tau,E,delta,fctE)*DtD(t,tau,E,delta,fctE)
    xi3=GammaMt(t,tau,GM,DGM,fctGM)
    xi4=rhot(t,tau,rho,Delta,fctR)*Et(t,tau,E,delta,fctE)*GammaDt(t, tau, GD, dGD, fctGD)
    xi=xi1+xi2+xi3+xi4
    return  xi

def thetat(t,tau,rho,E,GD,GM,Delta,delta,dGD,DGM,fctR,fctE,fctGD,fctGM):
    theta1=rhotD2(t,tau,rho,Delta,fctR)
    theta2=rhotD1(t,tau,rho,Delta,fctR)*Et(t,tau,E,delta,fctE)*DtD(t,tau,E,delta,fctE)
    theta3=Et(t,tau,E,delta,fctE)*DtD(t,tau,E,delta,fctE)*GammaMt(t,tau,GM,DGM,fctGM)
    theta4=GammaMtD(t,tau,GM,DGM,fctGM)
    theta5=rhotD1(t,tau,rho,Delta,fctR)*Et(t,tau,E,delta,fctE)*GammaDt(t, tau, GD, dGD, fctGD)
    theta6=Et(t,tau,E,delta,fctE)*GammaMt(t,tau,GM,DGM,fctGM)*GammaDt(t, tau, GD, dGD, fctGD)
    theta=theta1+theta2+theta3+theta4+theta5+theta6
    return theta


def SF(fct,n,tau,val,delta,fctMod,M):
    t=np.linspace(0,tau,M)
    dt=t[1]-t[0]
    an=dt/tau*np.sum((fct(t[1:],tau,val,delta,fctMod)+fct(t[:-1],tau,val,delta,fctMod))/2*np.exp(-2j*n*np.pi/tau*(t[1:]+t[:-1])/2))   #(fct(t[:-1])*np.exp(-2j*n*np.pi/tau*(t[:-1]+dt/2)))
    return an

def SFB(fct,n,tau,val1,val2,val3,val4,delta1,delta2,delta3,delta4,fctMod1,fctMod2,fctMod3,fctMod4,M):
    t=np.linspace(0,tau,M)
    dt=t[1]-t[0]
    an=dt/tau*np.sum((fct(t[1:],tau,val1,val2,val3,val4,delta1,delta2,delta3,delta4,fctMod1,fctMod2,fctMod3,fctMod4)+fct(t[:-1],tau,val1,val2,val3,val4,delta1,delta2,delta3,delta4,fctMod1,fctMod2,fctMod3,fctMod4))/2*np.exp(-2j*n*np.pi/tau*(t[1:]+t[:-1])/2))   #(fct(t[:-1])*np.exp(-2j*n*np.pi/tau*(t[:-1]+dt/2)))
    return an

def fillMatrix(tau,valR,valE,valGD,valGM,deltaR,deltaE,deltaGD,deltaGM,fctModR,fctModE,fctModGD,fctModGM,Nom=50,N=16,Ni=500,NE=24):
    MatrixA=np.zeros([2*N+1,2*N+1],dtype=complex)
    MatrixB=np.zeros([2*N+1,2*N+1],dtype=complex)
    VP=np.zeros([NE,Nom+1],dtype=complex)
   # VPI=np.zeros([2*N+1,Nom+1],dtype=complex)
    om=np.linspace(-np.pi/tau,np.pi/tau,Nom+1)
    for omeg in range(Nom+1):
        print(omeg)
        for n in range(-N,N+1):
            nn=n+N
            for p in range(-N,N+1):
                pp=p+N
                anE=SF(Et,p-n,tau,valE,deltaE,fctModE,Ni)
                anR=SF(rhot,p-n,tau,valR,deltaR,fctModR,Ni)
                anX=SFB(xit,p-n,tau,valR,valE,valGD,valGM,deltaR,deltaE,deltaGD,deltaGM,fctModR,fctModE,fctModGD,fctModGM,Ni)
                anT=SFB(thetat,p-n,tau,valR,valE,valGD,valGM,deltaR,deltaE,deltaGD,deltaGM,fctModR,fctModE,fctModGD,fctModGM,Ni)
                MatrixA[pp,nn]=(2*np.pi*n/tau+om[omeg])*(2*np.pi*p/tau+om[omeg])*anR-1j*(2*np.pi*n/tau+om[omeg])*anX-anT
                MatrixB[pp,nn]=anE
        EigV=ssl.eigs(MatrixA,k=NE,M=MatrixB,return_eigenvectors=False,which='SM')
        #EigVI=spl.eigvals(-MatrixA,k=NE,M=MatrixB,which='SM',return_eigenvectors=False)
        VP[:,omeg]=np.sqrt(EigV)
        #VPI[:,omeg]=np.sqrt(EigVI)
    VP=np.nan_to_num(VP)
    VP.sort(axis=0)
    #VPI.sort(axis=0)
    print("tau: ",tau)
    print("Omega: ",2*np.pi/tau)
    print("gamma: ",valGM/valR/(2*np.pi/tau))
    return om,VP
