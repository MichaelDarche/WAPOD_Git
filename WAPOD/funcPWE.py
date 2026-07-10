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
        Et=1/(1/E*(1+delta*np.sin((2*np.pi*t/tau))))
    if fctM=="np.cos":
        Et=1/(1/E*(1+delta*np.cos((2*np.pi*t/tau))))
    elif fctM=="np.sign":
        e0=np.sin((2*np.pi*t/tau))
        Et=1/(1/E*(1+delta*np.sign(e0)))
    elif fctM=="heavy":
        Et=1/(1/E*(1+delta*fctheavyside(t,tau)))
    return Et

def SF(fct,n,tau,val,delta,fctMod,M):
    t=np.linspace(0,tau,M)
    dt=t[1]-t[0]
    an=dt/tau*np.sum((fct(t[1:],tau,val,delta,fctMod)+fct(t[:-1],tau,val,delta,fctMod))/2*np.exp(-2j*n*np.pi/tau*(t[1:]+t[:-1])/2))   #(fct(t[:-1])*np.exp(-2j*n*np.pi/tau*(t[:-1]+dt/2)))
    return an

def fillMatrix(tau,valE,valR,deltaE,deltaR,fctModE,fctModR,NE,Nom=50,N=16,Ni=500):
    MatrixA=np.zeros([2*N+1,2*N+1],dtype=complex)
    MatrixB=np.zeros([2*N+1,2*N+1],dtype=complex)
    VP=np.zeros([NE,Nom+1],dtype=complex)
    VPI=np.zeros([NE,Nom+1],dtype=complex)
    om=np.linspace(-np.pi/tau,np.pi/tau,Nom+1)
    for omeg in range(Nom+1):
        for n in range(-N,N+1):
            nn=n+N
            for p in range(-N,N+1):
                pp=p+N
                anE=SF(Et,p-n,tau,valE,deltaE,fctModE,Ni)
                anR=SF(rhot,p-n,tau,valR,deltaR,fctModR,Ni)
                MatrixA[pp,nn]=(2*np.pi*n/tau-om[omeg])*(2*np.pi*p/tau-om[omeg])*anR
                MatrixB[pp,nn]=anE
        EigV=ssl.eigsh(MatrixA,k=NE,M=MatrixB,which='SM',return_eigenvectors=False)
        EigVI=ssl.eigs(-MatrixA,k=NE,M=MatrixB,which='SM',return_eigenvectors=False)
        VP[:,omeg]=np.sqrt(EigV)
        VPI[:,omeg]=np.sqrt(EigVI)
    VP=np.nan_to_num(VP) 
    VP.sort(axis=0)
    VPI.sort(axis=0)
    return om,VP

def fillMatrixS(lambdaMod,valE,valR,deltaE,deltaR,fctModE,fctModR,NE,Nk=50,N=16,Ni=500):
    MatrixA=np.zeros([2*N+1,2*N+1],dtype=complex)
    MatrixB=np.zeros([2*N+1,2*N+1],dtype=complex)
    VP=np.zeros([NE,Nk+1],dtype=complex)
    VPI=np.zeros([NE,Nk+1],dtype=complex)
    k=np.linspace(-np.pi/lambdaMod,np.pi/lambdaMod,Nk+1)
    for ki in range(Nk+1):
        for n in range(-N,N+1):
            nn=n+N
            for p in range(-N,N+1):
                pp=p+N
                anE=SF(Et,p-n,lambdaMod,valE,deltaE,fctModE,Ni)
                anR=SF(rhot,p-n,lambdaMod,valR,deltaR,fctModR,Ni)
                MatrixA[pp,nn]=(2*np.pi*n/lambdaMod-k[ki])*(2*np.pi*p/lambdaMod-k[ki])*anR
                MatrixB[pp,nn]=anE
        EigV=ssl.eigsh(MatrixA,k=NE,M=MatrixB,which='SM',return_eigenvectors=False)
        EigVI=ssl.eigs(-MatrixA,k=NE,M=MatrixB,which='SM',return_eigenvectors=False)
        VP[:,ki]=np.sqrt(EigV)
        VPI[:,ki]=np.sqrt(EigVI)
    VP=np.nan_to_num(VP) 
    VP.sort(axis=0)
    VPI.sort(axis=0)
    return k,VP