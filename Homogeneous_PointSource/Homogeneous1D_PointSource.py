#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:50:43 2023

@author: michael
"""

import os
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/Homogeneous_PointSource/")

import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../WAPOD/")

import scheme1D_TDold as scheme1D
import source
import inputReading


################# Input files 
#% Demarrer file
configuration=inputReading.readfile('Demarrer.txt')
mat=inputReading.readfile('Milieu.txt')
frontiere=inputReading.readfile('Frontiere.txt')     
sourcef=inputReading.readfile('Source.txt')   
################# Parameters of the simulation   
X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(configuration)
fs,CFL,x_0=inputReading.source(sourcef)
Nmat,rho,cm=inputReading.material(mat)
alpha=inputReading.frontiere(frontiere)
rhox,Celx=inputReading.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
dt=scheme1D.timeStep(CFL,dx,cm)
# lambdaOnde=cm/fs
# xsmin=cm*tshift-lambdaOnde
# xsmax=cm*tshift

#%

#% Initial conditions
U0=scheme1D.initPointSourceProblem(configuration,frontiere,mat,sourcef)



###### Mise en oeuvre
Nt=10


Tmax=Nt*dt
tmax=str(Tmax)+"s"

t=np.linspace(0,Tmax,Nt+1)
# x_0=120
sce=np.zeros(Nt+1)
for i in range(Nt+1):
    sce[i]=source.choice_timefct(sourcef,t[i])
    
K=2
Res,V,Sig=scheme1D.FD_sourceV(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0,sce,5)

# plt.figure()
# plt.plot(X,U0[0,:])
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v$ (m/s)")
# plt.title("LW : $t=$0 s")
# # plt.savefig('Figures/LW-V_0.eps', format='eps')

# plt.figure()
# plt.plot(X,U0[1,:])
# plt.xlabel("$x$ (m)")
# plt.ylabel("$\sigma$ (Pa)")
# plt.title("LW : $t=$0 s")
# # plt.savefig('Figures/LW-S_0.eps', format='eps')

plt.figure()
plt.plot(X,Res[0,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$v$ (m/s)")
plt.title("LW : $t=$"+tmax)
# plt.savefig('Figures/LW-V_100.eps', format='eps')

plt.figure()
plt.plot(X,Res[1,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$\sigma$ (Pa)")
plt.title("LW : $t=$"+tmax)
# plt.savefig('Figures/LW-S_100.eps', format='eps')


# ###############################

THV0y=np.zeros(Nx)
THV100y=np.zeros(Nx)
THS0y=np.zeros(Nx)
THS100y=np.zeros(Nx)
for x in range(Nx):
    THV0y[x]=1/2/cm*source.choice_timefct(sourcef,-np.abs(x*dx-x_0)/cm)
    THS0y[x]=-rho/2*source.choice_timefct(sourcef,-np.abs(x*dx-x_0)/cm)*np.sign(x*dx-x_0)
    THV100y[x]=1/2/cm*source.choice_timefct(sourcef,Tmax-np.abs(x*dx-x_0)/cm)
    THS100y[x]=-rho/2*source.choice_timefct(sourcef,Tmax-np.abs(x*dx-x_0)/cm)*np.sign(x*dx-x_0)
    

# plt.figure()
# plt.plot(X,U0[0,:])
# plt.plot(X,THV0y.T,'--')
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v$ (m/s)")
# plt.title("LW : $t=0 s$")
# # plt.savefig('Figures/Comp_LW-V_0.eps', format='eps')

# plt.figure()
# plt.plot(X,U0[0,:]-THV0y.T)
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v-v_{th}$ (m/s)")
# plt.title("LW : $t=0 s$")
# # plt.savefig('Figures/Diff_LW-V_0.eps', format='eps')

# plt.figure()
# plt.plot(X,U0[1,:])
# plt.plot(X,THS0y.T,'--')
# plt.xlabel("$x$ (m)")
# plt.ylabel("$\sigma$ (Pa)")
# plt.title("LW : $t=0 s$")
# # plt.savefig('Figures/Comp_LW-S_0.eps', format='eps')

# plt.figure()
# plt.plot(X,U0[1,:]-THS0y.T)
# plt.xlabel("$x$ (m)")
# plt.ylabel("$\sigma-\sigma_{th}$ (Pa)")
# plt.title("LW : $t=0 s$")
# # plt.savefig('Figures/Diff_LW-S_0.eps', format='eps')

plt.figure()
plt.plot(X,Res[0,:])
plt.plot(X,THV100y.T,'--')
plt.xlabel("$x$ (m)")
plt.ylabel("$v$ (m/s)")
plt.title("LW : $t=$"+tmax)
# plt.savefig('Figures/Comp_LW-V_100.eps', format='eps')

plt.figure()
plt.plot(X,Res[0,:]-THV100y.T)
plt.xlabel("$x$ (m)")
plt.ylabel("$v-v_{th}$ (m/s)")
plt.title("LW : $t=$"+tmax)
# plt.savefig('Figures/Diff_LW-V_100.eps', format='eps')

plt.figure()
plt.plot(X,Res[1,:])
plt.plot(X,THS100y.T,'--')
plt.xlabel("$x$ (m)")
plt.ylabel("$\sigma$ (m/s)")
plt.title("LW : $t=$"+tmax)
# plt.savefig('Figures/Comp_LW-S_100.eps', format='eps')

plt.figure()
plt.plot(X,Res[1,:]-THS100y.T)
plt.xlabel("$x$ (m)")
plt.ylabel("$\sigma-\sigma_{th}$ (m/s)")
plt.title("LW : $t=$"+tmax)
# plt.savefig('Figures/Diff_LW-S_100.eps', format='eps')


K=4
Res,V0,S0=scheme1D.FD_sourceV(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0,sce,1)

# plt.figure()
# plt.plot(X,U0[0,:])
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v$ (m/s)")
# plt.title("ADER4 : $t=$0 s")
# # plt.savefig('Figures/ADER4-V_0.eps', format='eps')

# plt.figure()
# plt.plot(X,U0[1,:])
# plt.xlabel("$x$ (m)")
# plt.ylabel("$\sigma$ (Pa)")
# plt.title("ADER4 : $t=$0 s")
# # plt.savefig('Figures/ADER4-S_0.eps', format='eps')


plt.figure()
plt.plot(X,Res[0,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$v$ (m/s)")
plt.title("ADER4 : $t=$"+tmax)
# plt.savefig('Figures/ADER4-V_100.eps', format='eps')

plt.figure()
plt.plot(X,Res[1,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$\sigma$ (Pa)")
plt.title("ADER4 : $t=$"+tmax)
# plt.savefig('Figures/ADER4-S_100.eps', format='eps')



# plt.figure()
# plt.plot(X,U0[0,:])
# plt.plot(X,THV0y.T,'--')
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v$ (m/s)")
# plt.title("ADER4 : $t=0 s$")
# # plt.savefig('Figures/Comp_ADER4-V_0.eps', format='eps')

# plt.figure()
# plt.plot(X,U0[0,:]-THV0y.T)
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v-v_{th}$ (m/s)")
# plt.title("ADER4 : $t=0 s$")
# # plt.savefig('Figures/Diff_ADER4-V_0.eps', format='eps')

# plt.figure()
# plt.plot(X,U0[1,:])
# plt.plot(X,THS0y.T,'--')
# plt.xlabel("$x$ (m)")
# plt.ylabel("$\sigma$ (Pa)")
# plt.title("ADER4 : $t=0 s$")
# # plt.savefig('Figures/Comp_ADER4-S_0.eps', format='eps')

# plt.figure()
# plt.plot(X,U0[1,:]-THS0y.T)
# plt.xlabel("$x$ (m)")
# plt.ylabel("$\sigma-\sigma_{th}$ (Pa)")
# plt.title("ADER4 : $t=0 s$")
# # plt.savefig('Figures/Diff_ADER4-S_0.eps', format='eps')

plt.figure()
plt.plot(X,Res[0,:])
plt.plot(X,THV100y.T,'--')
plt.xlabel("$x$ (m)")
plt.ylabel("$v$ (m/s)")
plt.title("ADER4 : $t=$"+tmax)
# plt.savefig('Figures/Comp_ADER4-V_100.eps', format='eps')

plt.figure()
plt.plot(X,Res[0,:]-THV100y.T)
plt.xlabel("$x$ (m)")
plt.ylabel("$v-v_{th}$ (m/s)")
plt.title("ADER4 : $t=$"+tmax)
# plt.savefig('Figures/Diff_ADER4-V_100.eps', format='eps')

plt.figure()
plt.plot(X,Res[1,:])
plt.plot(X,THS100y.T,'--')
plt.xlabel("$x$ (m)")
plt.ylabel("$\sigma$ (m/s)")
plt.title("ADER4 : $t=$"+tmax)
# plt.savefig('Figures/Comp_ADER4-S_100.eps', format='eps')

plt.figure()
plt.plot(X,Res[1,:]-THS100y.T)
plt.xlabel("$x$ (m)")
plt.ylabel("$\sigma-\sigma_{th}$ (m/s)")
plt.title("ADER4 : $t=$"+tmax)
# plt.savefig('Figures/Diff_ADER4-S_100.eps', format='eps')
