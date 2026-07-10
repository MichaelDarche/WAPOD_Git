#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 13:55:10 2025

@author: michael
"""


import os
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/1-Interface_Imerged_Cauchy_modulated_enerfm075")
#os.chdir("/home/darche/Documents/ModulatedMetamaterials/WAPRO_TD/1-Interface_Imerged_Cauchy_modulated")

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../WAPOD/")
sys.path.append("/home/darche/Documents/ModulatedMetamaterials/WAPRO_TD/WAPOD/")
import scheme1D_TDold as scheme1D
import source
import inputReading

import modulation as modul
import modulation as modulateP

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


mat=inputReading.readfile('Milieu.txt')
frontiere=inputReading.readfile('Frontiere.txt')     
sourcef=inputReading.readfile('Source.txt')   
# sourceRev=inputReading.readfile('SourceRev.txt')
# nomodulation=inputReading.readfile('NoModulation.txt')
################ Parameters of the simulation   

fs,CFL,x_0=inputReading.source(sourcef)
# fsR,CFLR,x_0R=inputReading.source(sourceRev)
Nmat,rho,cm=inputReading.material(mat)
alpha=inputReading.frontiere(frontiere)


###### Mise en oeuvre
# Nt=1000
# Tmax=Nt*dt
# tmax=str(Tmax)+"s"
K=4
# t=np.linspace(0,Tmax,Nt+1)




pFilm=1
Nmod=10
mod075=[0.]
diffener075=[0.0]

#for i in range(1,25):
for i in [24]:
    namecof="Demarrer"+str(i-1)+".txt"
    configuration=inputReading.readfile(namecof)
    namemod="Modulation"+str(i-1)+".txt"
    modulation=inputReading.readfile(namemod)
    X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(configuration)
    rhox,Celx=inputReading.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
    dt=min(scheme1D.timeStep(CFL,dx,cm))
    
    Tmax=0.04
    Ntfig=int(np.floor(Tmax/dt))
    t=np.linspace(0,Tmax,Ntfig+1)
    Nt=Ntfig
    sce=np.zeros(Nt+1)
    for i in range(Nt+1):
        sce[i]=source.choice_timefct(sourcef,t[i])
    # lambdaOnde=cm/fs
    # xsmin=cm*tshift-lambdaOnde
    # xsmax=cm*tshift
    #%

    #% Initial conditions
    U0=scheme1D.initCauchyProblem(configuration,frontiere,mat,sourcef)
    import time 
    time0=time.time()
    Res,Vfilm,Sfilm,CSVM=scheme1D.FD_cauchyII(U0,Ntfig+1,Nx,K,rhox,Celx,dt,dx,configuration,mat,frontiere,modulation,pFilm,rCSVM="yes")
    #U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
    timeM=time.time()-time0
    maxV=np.max([np.max(Vfilm),-np.min(Vfilm)])
    print(timeM)  
    if float(frontiere["Kn_0"])==0:
        K0=np.infty
    else:
        K0=float(frontiere["Kn_0"])
    frK=float(modulation["Freq"])
    delta=float(modulation["DeltaF"])

    M0=float(frontiere["Mn_0"])
    frM=float(modulation["Freq"])
    deltaM=float(modulation["DeltaM"])

    modulb=modulateP.fonctionModM(t,alpha[0],modulation,frontiere,0)
    modulM=M0*(1+modulb)
    dmodulb=np.zeros(int(Nt/pFilm))
    # for i in range(int(Nt/pFilm)):
    #     modulb[i]=M0*(1+deltaM*np.sin(2*np.pi*frM*t[pFilm*i]))
    #     dmodulb[i]=2*np.pi*frM*M0*(deltaM*np.cos(2*np.pi*frM*t[pFilm*i]))
        
    modula=modulateP.fonctionModK(t,alpha[0],modulation,frontiere,0)
    modulK=1/K0*(1+modula)
    dmodula=np.zeros(int(Ntfig/pFilm))

    energyV=np.zeros(int(Ntfig/pFilm))
    energyI=np.zeros(int(Ntfig/pFilm))
    for ti in range(int(Ntfig/pFilm)):
        energyV[ti]=dx/2*np.sum(rhox*Vfilm[ti]*Vfilm[ti])+dx/2*np.sum(1/rhox/Celx/Celx*Sfilm[ti]*Sfilm[ti])
        energyI[ti]=1/2*modulM[ti]*CSVM[ti,2]*CSVM[ti,2]+1/2*modulK[ti]*CSVM[ti,3]*CSVM[ti,3]     
        
    energyM=energyV+energyI
    mod075.append(frM)
    diffe075=energyM[-1]-energyM[0]
    diffener075.append(diffe075)