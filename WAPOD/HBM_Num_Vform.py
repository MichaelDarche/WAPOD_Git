#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 14:14:57 2025

@author: michael
"""


import numpy as np

def sysPsi(n,dK,rho,c,K,f,F,Qc,dB):
    if K!=0:
        Km1=1/K
    else:
        Km1=0
    sK=rho*c*Km1
    beta=Qc*rho*c
    APsi=np.zeros([2*n+1,2*n+1],dtype='complex')
    BPsi=np.zeros(2*n+1,dtype='complex')
    for i in range(-n,n+1):
        k=i+n
        ok=2*np.pi*(f+i*F)
        APsi[k,k]=1+beta/2+1j*ok/2*sK
        if i==-n:
            APsi[k,k+1]=-(sK*dK*ok-1j*beta*dB)/4
        elif i==n:
            APsi[k,k-1]=(sK*dK*ok-1j*beta*dB)/4
        else:
            APsi[k,k+1]=-(sK*dK*ok-1j*beta*dB)/4
            APsi[k,k-1]=(sK*dK*ok-1j*beta*dB)/4
        if i==-1:
            BPsi[k]=(sK*dK*ok-1j*beta*dB)/4
        elif i==0:
            BPsi[k]=1-beta/2-1j*sK/2*ok
        elif i==1:
            BPsi[k]=-(sK*dK*ok-1j*beta*dB)/4
    return APsi,BPsi

def sysPhi(n,dM,rho,c,M,f,F,Qm,dZ):
    sM=M/rho/c
    zeta=Qm/rho/c
    APhi=np.zeros([2*n+1,2*n+1],dtype='complex')
    BPhi=np.zeros(2*n+1,dtype='complex')
    for i in range(-n,n+1):
        k=i+n
        ok=2*np.pi*(f+i*F)
        APhi[k,k]=1+zeta/2+1j*ok/2*sM
        if i==-n:
            APhi[k,k+1]=-(sM*dM*ok-1j*dZ*zeta)/4
        elif i==n:
            APhi[k,k-1]=(sM*dM*ok-1j*dZ*zeta)/4
        else:
            APhi[k,k+1]=-(sM*dM*ok-1j*dZ*zeta)/4
            APhi[k,k-1]=(sM*dM*ok-1j*dZ*zeta)/4
        if i==-1:
            BPhi[k]=(sM*dM*ok-1j*dZ*zeta)/4
        elif i==0:
            BPhi[k]=1-zeta/2-1j*sM/2*ok
        elif i==1:
            BPhi[k]=-(sM*dM*ok-1j*dZ*zeta)/4
    return APhi,BPhi

def RTcoefs(n,dM,dK,rho,c,M,K,f,F,Qc=0,dB=0,Qm=0,dZ=0):
    APsi,BPsi=sysPsi(n,dK,rho,c,K,f,F,Qc,dB)
    APhi,BPhi=sysPhi(n,dM,rho,c,M,f,F,Qm,dZ)
    NPsi,NBPsi=normalizeMatrix(APsi,BPsi)
    NPhi,NBPhi=normalizeMatrix(APhi,BPhi)
    Psi=np.linalg.solve(NPsi,NBPsi)
    Phi=np.linalg.solve(NPhi,NBPhi)
    R=(Phi-Psi)/2
    T=(Phi+Psi)/2
    freq=vecFreq(n,f,F)
    return freq,R,T

def vecFreq(n,f,F):
    freqs=np.zeros(2*n+1)
    for k in range(-n,n+1):
        freqs[k+n]=f+k*F
    return freqs

def normV(V,n):
    VN=np.zeros(2*n+1)
    for i in range(-n,n+1):
        VN[i+n]=V[i+n]/V[n]
    return VN

def normalizeMatrix(Mat0,Mat1):
    maxi=0.
    Mat0N=np.zeros(Mat0.shape)
    Mat1N=np.zeros(Mat1.shape)
    maxi0=np.abs(Mat0.max())
    maxi1=np.abs(Mat1.max())
    mini0=np.abs(Mat0.min())
    mini1=np.abs(Mat1.min())
    maxi=max(maxi0,maxi1,mini0,mini1)
    # for i in range(Mat0.shape[0]):
        # for j in range(Mat0.shape[0]):
        #     if np.abs(Mat0[i,j])>maxi:
        #         maxi=np.abs(Mat0[i,j])
        #     if np.abs(Mat1[i,j])>maxi:
        #         maxi=np.abs(Mat1[i,j])
    if maxi !=0.:
        Mat0N=Mat0/maxi
        Mat1N=Mat1/maxi
    return Mat0N, Mat1N