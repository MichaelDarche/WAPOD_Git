#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 21:02:41 2023

@author: michael
"""
import numpy as np 

def transfertM(rho,cm,omega,l):
    k=omega/cm
    Z=rho*cm
    Y=Z*omega
    T=np.zeros([2,2])
    T[0,0]=np.cos(k*l)
    T[1,1]=np.cos(k*l)
    T[0,1]=1/Y*np.sin(k*l)
    T[1,0]=-Y*np.sin(k*l)
    return T
    
def matT(omega,Ncouches,li,rho,cm):
    P=np.identity(2)
    for i in range(Ncouches):
        T=transfertM(rho[i+1], cm[i+1], omega, li[i])
        P=np.dot(P, T)
    return P

def matSM(K,M,omega):
    omega0=2*np.sqrt(K/M)
    Jmoins=np.array([[1,1/2/K],[-M/2*omega**2,1]])
    Jplusm1=1/(1+(omega/omega0)**2)*Jmoins
    SM=np.dot(Jplusm1,Jmoins)
    return SM

def matTJ(omega,Ncouches,li,rho,cm,K,M):
    P=np.identity(2)
    J=matSM(K, M, omega)
    for i in range(Ncouches):
        T=transfertM(rho[i+1], cm[i+1], omega, li[i])
        P=np.dot(P, T)
    return np.dot(P,J)

def blochVP(P,L):
    eV1,eV2=np.linalg.eigvals(P)
    q1=np.log(eV1)/L/1j
    q2=np.log(eV2)/L/1j
    return q1,q2