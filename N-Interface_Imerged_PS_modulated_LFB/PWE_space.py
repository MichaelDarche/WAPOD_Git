#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 14:01:07 2024

@author: michael
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as spl

#Parametres materiaux
rho=1250
Delta=0.2
E=1714285714
delta=-0.714285714285714
#Periode spatiale
h=10


def rhoh(x,h=h,rho=rho,Delta=Delta):
    rhoh=rho+Delta*np.sign(np.sin(2*np.pi*x/h))
    
    return rhoh
#
def Eh(x,h=h,E=E,delta=delta):
    Eh=E+delta*np.sign(np.sin(2*np.pi*x/h))
    return Eh


# Series de Fourier
def SF(fct,n,h,M):
    x=np.linspace(-h/2,h/2,M)
    dx=x[1]-x[0]
    plt.plot(x,fct(x))
    an=dx/h*np.sum(fct(x)*np.exp(-2j*n*np.pi/h*x))
    return an



# Bloch vector
Nk=100
k=np.linspace(-np.pi/h,np.pi/h,Nk)

#Taille de la SF
N=2**10
#Integration pour les coefs
M=500
#Nombre de VP
NV=12

# Initialisation des matrices
MatrixA=np.zeros([2*N,2*N])
MatrixB=np.zeros([2*N,2*N])
VP=np.zeros([2*N,Nk])



# Remplissage des matrices
for kk in range(Nk):
    for p in range(-N,N+1):
        for n in range(-N,N+1):
            MatrixA[p,n]=(2*np.pi*n/h+k[kk])*(-2*np.pi*p/h+k[kk])*SF(Eh,-p-n,h,M)
            MatrixB[p,n]=SF(rhoh,-p-n,h,M)
    EigV=spl.eigvals(MatrixA,MatrixB)
    VP[:,kk]=EigV
    


VP.sort(axis=0)

plt.figure()
for i in range(2*N):
    plt.plot(k,np.sqrt(VP[i,:]))

