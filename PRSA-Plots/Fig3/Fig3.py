#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 16:53:54 2026

@author: michael
"""


import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

import sys
import os

sys.path.append(os.path.abspath(".."))

import functions as func

##################
## Load parameters
##################
configuration=func.readfile('Demarrer.txt')
mat=func.readfile('Milieu.txt')
frontiere=func.readfile('Frontiere.txt')     
sourcef=func.readfile('Source.txt')  

X,xmin,xmax,Nx,dx,ESIM,Npt=func.geometric(configuration)
fs,CFL,x_0=func.source(sourcef)
Nmat,rho,cm=func.material(mat)
alpha=func.frontiere(frontiere)
rhox,Celx=func.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
fm=float(frontiere["Freq_0"]) 



##################
## Load datas
##################

nSveRes='Vfilm'+str(fs)+'fm'+str(fm)+'.txt'
nSveHomo='Uhomo'+str(fs)+'fm'+str(fm)+'.txt'
nSvetp='tp'+str(fs)+'fm'+str(fm)+'.txt'
 
Vfilm=np.loadtxt(nSveRes, delimiter='\t')
U=np.loadtxt(nSveHomo, delimiter='\t')
tp=np.loadtxt(nSvetp,delimiter='\t')
dt=tp[1]-tp[0]
Ufilm=func.VtoU(Vfilm,dt)


##################
## Make plots
##################
Ntfig=int(np.floor(0.2/dt))

psmax=np.where(Ufilm[Ntfig-1,:Nx//2]==max(Ufilm[Ntfig-1,:Nx//2]))
xzoommin=X[psmax[0][0]-200]
xzoommax=X[psmax[0][0]+600]
zoomlim=[xzoommin,xzoommax]


plt.figure()
plt.axvspan(xzoommin, xzoommax, facecolor='g', alpha=0.3)
plt.scatter(X[::1],Ufilm[Ntfig-1,::1],marker=r'$\bigcirc$',label=r'$U_h$',color='red',s=30)
plt.plot(X,U[Ntfig-1,:],linewidth=3,label=r'$U_0$')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.title(r'$t=$ '+str(round(tp[Ntfig-1], 3))+" s",fontsize=18)
plt.legend(fontsize=16)
plt.tight_layout()
plt.show()



plt.figure()
for i in range(alpha.size):
    if frontiere["Contact_"+str(i)]!="parfait":
        plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(X[::1],Ufilm[Ntfig-1,::1],marker=r'$\bigcirc$',label=r'$U_h$',color='red',s=30)
plt.plot(X,U[Ntfig-1,:],linewidth=3,label=r'$U_0$')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.title(r'$t=$ '+str(round(tp[Ntfig-1], 3))+" s",fontsize=18)
plt.xlim(zoomlim)
plt.legend(fontsize=16)#loc="lower right",
plt.tight_layout()
plt.show()

