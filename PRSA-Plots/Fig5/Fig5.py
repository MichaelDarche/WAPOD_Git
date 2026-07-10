#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 18:01:38 2025

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

for isource in range(3):
    configuration=func.readfile('Demarrer.txt')
    mat=func.readfile('Milieu.txt')
    frontiere=func.readfile('Frontiere.txt')     
    sourcef=func.readfile('Source'+str(isource)+'.txt')   
    modulation=func.readfile('Modulation.txt')
    
    X,xmin,xmax,Nx,dx,ESIM,Npt=func.geometric(configuration)
    fs,CFL,x_0=func.source(sourcef)
    Nmat,rho,cm=func.material(mat)
    alpha=func.frontiere(frontiere)
    rhox,Celx=func.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
    fm=float(modulation["Freq"]) 
    
    skipSpacePoints=2
    
    nSveUz='Uz_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.txt'
    nSveU='U_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.txt'
    nSveP1='P1_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.txt'
    nSveP2='P2_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.txt'
    nSvetp='tp_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.txt'
    nSveVfilm='Vfilm_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.txt'
     
    U=np.loadtxt(nSveU, delimiter='\t')
    Uz=np.loadtxt(nSveUz, delimiter='\t')
    Vfilm=np.loadtxt(nSveVfilm,delimiter='\t')
    tp=np.loadtxt(nSvetp,delimiter='\t')
    P1=np.loadtxt(nSveP1,delimiter='\t')
    P2=np.loadtxt(nSveP2,delimiter='\t')
    

    dt=tp[1]-tp[0]
    testcfl=Celx[0]*dt/dx
    
    Ufilm=func.VtoU(Vfilm, dt)
    
    Tfig=0.2 # Time of the figure
    Ntfig=int(np.ceil(Tfig/dt)) # Corresponding index
    U2C=func.postCorrector(Ntfig-1,dx,Nx,U,P1,P2) # Modification with the order2 correctors
    
    # Position of the zoom
    psmax=np.where(Ufilm[Ntfig-1,:Nx//2]==max(Ufilm[Ntfig-1,:Nx//2]))
    xzoommin=X[psmax[0][0]-200]
    xzoommax=X[psmax[0][0]+600]
    zoomlim=[xzoommin,xzoommax]
    
    # Displacement fields
    plt.figure()
    plt.axvspan(xzoommin, xzoommax, facecolor='g', alpha=0.3)
    plt.scatter(X[::skipSpacePoints],Ufilm[Ntfig,::skipSpacePoints],marker=r'$\bigcirc$',label=r'$U_h$',color='red',s=40)
    plt.plot(X[::skipSpacePoints],Uz[Ntfig-1,::skipSpacePoints],'grey',linewidth=3,label=r'$U_0$')
    plt.plot(X[::skipSpacePoints],U2C[::skipSpacePoints],color='blue',linestyle='None',marker='o',markersize=3,label=r'$U^{(2)}$')
    plt.xlabel(r'$X$ (m)',fontsize=16)
    plt.ylabel(r'$U$ (m)',fontsize=16)
    plt.title("$T=$"+str(round(tp[Ntfig], 3))+" s",fontsize=16)
    plt.legend(fontsize=14)#loc="lower right",)
    plt.tight_layout()
    plt.show()

    # Zoom
    plt.figure()
    for i in range(alpha.size):
        plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
    plt.scatter(X[::skipSpacePoints],Ufilm[Ntfig,::skipSpacePoints],marker=r'$\bigcirc$',label=r'$U_h$',color='red',s=40)
    plt.plot(X[::skipSpacePoints],Uz[Ntfig-1,::skipSpacePoints],'grey',linewidth=3,label=r'$U_0$')
    plt.plot(X[::skipSpacePoints],U2C[::skipSpacePoints],color='blue',linestyle='None',marker='o',markersize=3,label=r'$U^{(2)}$')
    plt.xlabel(r'$X$ (m)',fontsize=16)
    plt.ylabel(r'$U$ (m)',fontsize=16)
    plt.xlim(zoomlim)
    plt.title("$T=$"+str(round(tp[Ntfig], 3))+" s",fontsize=16)
    plt.legend(fontsize=14)#loc="lower right",)
    plt.tight_layout()
    plt.show()



    
