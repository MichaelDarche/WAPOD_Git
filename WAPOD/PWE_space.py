#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 14:01:07 2024

@author: michael
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as spl
import scipy.sparse.linalg as ssl
#Parametres materiaux
rho=1000
Delta=0.0
c=2800
E=rho*c**2
E=1.6e9
delta=0.6
#Periode spatiale
h=60
phi=1/3


def fctheavyside(x,h=h,phi=phi):
    Xphi=h*phi+x[0]
    res=x.copy()
    for i in range(x.shape[0]):
        if x[i]<Xphi:
            res[i]=-1
        else:
            res[i]=+1
    return res    


fctM="heavy"
#fct=
#fct=fctheavyside


def rhoh(x,h=h,rho=rho,Delta=Delta):
    if fctM=="np.sinplus":
        rhoh=rho*(1+Delta*np.sin((2*np.pi*x/h)))
    elif fctM=="np.sin":
        rhoh=rho*(1+Delta*np.sin((2*np.pi*x/h)))
    elif fctM=="np.cos":
        rhoh=rho*(1+Delta*np.cos((2*np.pi*x/h)))
    elif fctM=="np.sinsign":
        fctcomp=np.sin((2*np.pi*x/h))
        for i in range(3):
            fctcomp=fctcomp+np.sin((2*i+1)*(2*np.pi*x/h))
        rhoh=rho*(1+Delta*fctcomp)
    elif fctM=="np.sign":
        rhoh0=np.sin(2*np.pi*x/h)
        rhoh=rho*(1+Delta*np.sign(rhoh0))
    elif fctM=="heavy":
        rhoh=rho*(1+Delta*fctheavyside(x,h))
    return rhoh
#


def Eh(x,h=h,E=E,delta=delta):
    if fctM=="np.sin":
        Eh=1/(1/E*(1+delta*np.sin((2*np.pi*x/h))))
    elif fctM=="np.cos":
        Eh=1/(1/E*(1+delta*np.cos((2*np.pi*x/h))))
    elif fctM=="np.sinplus":
        Eh=E*(1+delta*np.sin((2*np.pi*x/h)))
    elif fctM=="np.sinsign":
        fctcomp=np.sin((2*np.pi*x/h))
        for i in range(3):
            fctcomp=fctcomp+np.sin((2*i+1)*(2*np.pi*x/h))
        Eh=1/(1/E*(1+delta*fctcomp))
    elif fctM=="np.sign":
        e0=np.sin((2*np.pi*x/h))
        Eh=1/(1/E*(1+delta*np.sign(e0)))
    elif fctM=="heavy":
        Eh=1/(1/E*(1+delta*fctheavyside(x,h)))
    return Eh


# Series de Fourier
def SF(fct,n,h,M):
    x=np.linspace(-h/2,h/2,M)
    dx=x[1]-x[0]
    an=dx/h*np.sum((fct(x[1:])+fct(x[:-1]))/2*np.exp(-2j*n*np.pi/h*(x[1:]+x[:-1])/2))
    print(an)
    return an



# Bloch vector
Nk=30
k=np.linspace(-np.pi/h,0,Nk+1)

#Taille de la SF
N=8#2**4
#Integration pour les coefs
Ni=600
#Nbre de modes
NE=8

# Initialisation des matrices
MatrixA=np.zeros([2*N+1,2*N+1],dtype=complex)
MatrixB=np.zeros([2*N+1,2*N+1],dtype=complex)
VP=np.zeros([NE,Nk+1])

# Remplissage des matrices

for kk in range(Nk+1):
    print(kk)
    for n in range(0,2*N+1):
        nn=n-N
        for p in range(0,2*N+1):
            pp=p-N
            anE=SF(Eh,pp-nn,h,Ni)
            anR=SF(rhoh,pp-nn,h,Ni)
            MatrixA[p,n]=(2*np.pi*nn/h-k[kk])*(2*np.pi*pp/h-k[kk])*anE
            MatrixB[p,n]=anR
    # print(spl.ishermitian(MatrixA),spl.ishermitian(MatrixA))
    # Bm1=spl.inv(MatrixB)
    # Mat=np.dot(Bm1,MatrixA)
    EigVS=ssl.eigs(MatrixA,k=NE,M=MatrixB,which='SM',return_eigenvectors=False)
    #EigVS2=spl.eig(a=MatrixA,b=MatrixB,left=False,right=False)
    # EigV=spl.eigvals(Mat)
    VP[:,kk]=np.sqrt(EigVS)
    # VP2[:,kk]=np.sqrt(EigV)
    
recoE=np.zeros(Ni)
recoR=np.zeros(Ni)
x=np.linspace(-h/2,h/2,Ni)
dx=x[1]-x[0]
for n in range(-N,N+1):
    nn=n
    #pp=p-N
    SFpn=SF(rhoh,nn,h,Ni)
    SFEn=SF(Eh,nn,h,Ni)
    for xi in range(Ni):
        recoE[xi]=recoE[xi]+SFEn*np.exp(2j*nn*np.pi/h*(x[xi]))
        recoR[xi]=recoR[xi]+SFpn*np.exp(2j*nn*np.pi/h*(x[xi]))
# # plt.figure()
# # plt.plot(recoE)


VP.sort(axis=0)
#VP2.sort(axis=0)
# plt.figure()
# for i in range(NE):
#     plt.plot(k,(VP2[i,:]))


##################################
#Comparaison


rho1=rho*(1-Delta)
E1=1/(1/E*(1-delta))

rho2=rho*(1+Delta)
E2=1/(1/E*(1+delta))

c1=np.sqrt(E1/rho1)
c2=np.sqrt(E2/rho2)

cr=(phi)*c1+(1-phi)*c2

plt.figure()
plt.plot(x,rhoh(x),label=r'$\rho(x)$')
plt.plot(x,recoR,label='PWE - 32 modes')
plt.legend(loc="lower right",fontsize=16)
plt.xlabel(r'$x$ (s)',fontsize=16)
plt.ylabel(r'$\rho$ (kg/m$^3$)',fontsize=16)
plt.tight_layout()

plt.figure()
plt.plot(x,np.sqrt(recoE/recoR),label='PWE - 32 modes')
plt.legend(loc="lower right",fontsize=16)
plt.xlabel(r'$x$ (s)',fontsize=16)
plt.ylabel(r'$\rho$ (kg/m$^3$)',fontsize=16)
plt.tight_layout()
# plt.plot(np.linspace(-h/2,-h/2+phi*h,Ni),np.linspace(rho1,rho1,Ni))
# plt.plot(np.linspace(-h/2+phi*h,h/2,Ni),np.linspace(rho2,rho2,Ni))#,color='gray',linestyle='--')
plt.figure()
plt.plot(x,recoE)
plt.plot(x,Eh(x))

# plt.plot(np.linspace(-h/2,-h/2+phi*h,Ni),np.linspace(E1,E1,Ni))
# plt.plot(np.linspace(-h/2+phi*h,h/2,Ni),np.linspace(E2,E2,Ni))#,color='gray',linestyle='--')


Z1=rho1*c1
Z2=rho2*c2

L1 = phi*h
L2 = (1-phi)*h 

def D_micro(k,omega):
    Dmicro=np.cos(k*h)-np.cos(omega/c1*L1)*np.cos(omega/c2*L2)+1/2*(Z1/Z2+Z2/Z1)*np.sin(omega/c1*L1)*np.sin(omega/c2*L2)
    return Dmicro

Nome=1000
Nkk=500

Om_plot = np.linspace(0,6*np.pi/h*cr,Nome) 
k_map = np.linspace(-np.pi/h,np.pi/h,Nkk)


DD_micro = np.zeros([Nome,Nkk])
DD_det = np.zeros([Nome,Nkk])


for ind in range(Nkk):
    for ind2 in range(Nome):
        DD_micro[ind2,ind] = D_micro(k_map[ind],Om_plot[ind2])
        #DD_det(ind2,ind) = Det_micro(k_map(ind),Omega_plot(ind2)


plt.figure()
plt.pcolor(k_map,Om_plot,np.log10(np.abs(DD_micro)),vmin=-2.5,vmax=0.5)
plt.colorbar()
plt.xlabel(r"$k$ (1/m)",fontsize=16)
plt.ylabel(r'$\omega$ (rad/s)',fontsize=16)
plt.tight_layout()
for i in range(NE):
    plt.plot(k,(VP[i,:]),'red')
# for i in range(NE):    
#     plt.plot(k*h/2/np.pi,(VP2[i,:])*h/cr/2/np.pi,'--')