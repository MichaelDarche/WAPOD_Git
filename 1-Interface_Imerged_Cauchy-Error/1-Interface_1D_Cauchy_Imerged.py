#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:50:43 2023

@author: michael
"""

import os
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/1-Interface_Imerged_Cauchy-Error/")

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../WAPOD/")

import scheme1D
import source
import inputReading



################# Input files 
#%
configuration=inputReading.readfile('Demarrer.txt')
mat=inputReading.readfile('Milieu.txt')
frontiere=inputReading.readfile('Frontiere.txt')     
sourcef=inputReading.readfile('Source.txt')   
################ Parameters of the simulation   
X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(configuration)
fs,CFL,tshift=inputReading.source(sourcef)
Nmat,rho,cm=inputReading.material(mat)
alpha=inputReading.frontiere(frontiere)
rhox,Celx=inputReading.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
dt=min(scheme1D.timeStep(CFL,dx,cm))
#%

#% Initial conditions
U0=scheme1D.initCauchyProblem(configuration,frontiere,mat,sourcef)

###### Mise en oeuvre
Nt=250
Tmax=Nt*dt
tmax=str(Tmax)+"s"
K=2
Res,Vfilm,Sfilm=scheme1D.FD_cauchyII(U0,Nt,Nx,K,rhox,Celx,dt,dx,configuration,mat,frontiere,5)

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

# plt.figure()
# plt.plot(X,Res[0,:],"*-")
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v$ (m/s)")
# plt.title("LW : $t=$"+tmax)
# plt.axvline(x=alpha[0],color='gray',linestyle='--')
# # plt.savefig('Figures/LW-V_300.eps', format='eps')

# plt.figure()
# plt.plot(X,Res[1,:],"*-")
# plt.xlabel("$x$ (m)")
# plt.ylabel("$\sigma$ (Pa)")
# plt.title("LW : $t=$"+tmax)
# plt.axvline(x=alpha[0],color='gray',linestyle='--')
# # plt.savefig('Figures/LW-S_300.eps', format='eps')

# from matplotlib.animation import FuncAnimation

# # Fonction pour mettre à jour le contenu du graphique à chaque image
# def update(frame):
#     plt.clf()
#     plt.plot(X, Vfilm[frame,:])#, label=f'Frame {frame}')
#     plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0008,0.0008], aspect='auto', alpha=0.35)
#     # plt.axvline(x=alpha[0],color='gray',linestyle='--')
#     plt.colorbar()
#     # plt.plot(X, VfilmH[frame,:])
#     plt.xlabel("x (m)")
#     plt.ylabel("v (m/s)")
#     plt.ylim([-0.0008,0.0008])
#     # plt.legend()


# total_frames = Vfilm.shape[0]
# fig, ax = plt.subplots()

# # Cree l'animation
# animation = FuncAnimation(fig, update, frames=total_frames, interval=100)
# # Affiche l'animation
# plt.show()
# nomAnim="Film_1I_"+sourcef["Frequence"]+"Hz.mp4"
# animation.save(nomAnim, writer='ffmpeg', fps=20)

# ###############################

THV0y=np.zeros(Nx)
THVTMy=np.zeros(Nx)
THS0y=np.zeros(Nx)
THSTMy=np.zeros(Nx)
R01=(rho[1]*cm[1]-rho[0]*cm[0])/(rho[1]*cm[1]+rho[0]*cm[0])
T01=2*rho[0]*cm[1]/(rho[1]*cm[1]+rho[0]*cm[0])
Nxalpha=int(alpha[0]/dx)
for x in range(Nxalpha):
    f0=source.choice_timefct(sourcef,tshift-x/cm[0])
    fr=source.choice_timefct(sourcef,tshift+x/cm[0]-2*alpha[0]/cm[0])
    THV0y[x]=1/cm[0]*f0-R01/cm[0]*fr
    THS0y[x]=-rho[0]*f0-R01*rho[0]*fr
    fT=source.choice_timefct(sourcef,Tmax+tshift-x/cm[0])
    frT=source.choice_timefct(sourcef,Tmax+tshift+x/cm[0]-2*alpha[0]/cm[0])
    THVTMy[x]=1/cm[0]*fT-R01/cm[0]*frT
    THSTMy[x]=-rho[0]*fT-R01*rho[0]*frT
for x in range(Nxalpha,Nx):
    ft=source.choice_timefct(sourcef,tshift-x/cm[1]+(1/cm[1]-1/cm[0])*alpha[0])
    ftT=source.choice_timefct(sourcef,Tmax+tshift-x/cm[1]+(1/cm[1]-1/cm[0])*alpha[0])
    THV0y[x]=1/cm[1]*T01*ft
    THS0y[x]=-rho[1]*T01*ft
    THVTMy[x]=1/cm[1]*T01*ftT
    THSTMy[x]=-rho[1]*T01*ftT



plt.figure()
plt.plot(X,U0[0,:])
plt.plot(X,THV0y,'--')
plt.xlabel("$x$ (m)")
plt.ylabel("$v_{th}$ (m/s)")
plt.title("LW : $t=T_{max}$")

plt.figure()
plt.plot(X,U0[1,:])
plt.plot(X,THS0y,'--')
plt.xlabel("$x$ (m)")
plt.ylabel("$v_{th}$ (m/s)")
plt.title("LW : $t=T_{max}$")

plt.figure()
plt.plot(X,Res[0,:])
plt.plot(X,THVTMy,'--')
plt.xlabel("$x$ (m)")
plt.ylabel("$v_{th}$ (m/s)")
plt.title("LW : $t=T_{max}$")

plt.figure()
plt.plot(X,Res[1,:])
plt.plot(X,THSTMy,'--')
plt.xlabel("$x$ (m)")
plt.ylabel("$v_{th}$ (m/s)")
plt.title("LW : $t=T_{max}$")

def errL2(A,B):
    err=0
    for l in range(max(B.shape)):
        eps=A[l]-B[l]
        err=err+eps**2
    return err

ndx=7
err=np.zeros([4,ndx])
k0=2**((ndx-1)/2)
K=2
for k in range(ndx):
    k0=k0/2
    Res=scheme1D.FD_cauchyII(U0,Nt,Nx,K,rhox,Celx,dt,k0*dx,configuration,mat,frontiere,"no")
    err[0,k]=errL2(Res[0,:], THVTMy)
    err[1,k]=errL2(Res[1,:], THSTMy)
K=4
k0=2**((ndx-1)/2)
for k in range(ndx):
    k0=k0/2
    Res=scheme1D.FD_cauchyII(U0,Nt,Nx,K,rhox,Celx,dt,k0*dx,configuration,mat,frontiere,"no")
    err[2,k]=errL2(Res[0,:], THVTMy)
    err[3,k]=errL2(Res[1,:], THSTMy)
    

    

# LWV300y = np.array(colonne2)

# plt.figure()
# plt.plot(X,U0[0,:])
# plt.plot(X,LWV0y.T,'--')
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v$ (m/s)")
# plt.title("LW : $t=0 s$")
# plt.savefig('Figures/Comp_LW-V_0.eps', format='eps')

# plt.figure()
# plt.plot(X,U0[0,:]-LWV0y.T)
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v-v_{Bruno}$ (m/s)")
# plt.title("LW : $t=0 s$")
# plt.savefig('Figures/Diff_LW-V_0.eps', format='eps')

# plt.figure()
# plt.plot(X,U0[1,:])
# plt.plot(X,LWS0y.T,'--')
# plt.xlabel("$x$ (m)")
# plt.ylabel("$\sigma$ (Pa)")
# plt.title("LW : $t=0 s$")
# plt.savefig('Figures/Comp_LW-S_0.eps', format='eps')

# plt.figure()
# plt.plot(X,U0[1,:]-LWS0y.T)
# plt.xlabel("$x$ (m)")
# plt.ylabel("$\sigma-\sigma_{Bruno}$ (Pa)")
# plt.title("LW : $t=0 s$")
# plt.savefig('Figures/Diff_LW-S_0.eps', format='eps')

# plt.figure()
# plt.plot(X,Res[0,:])
# plt.plot(X,LWV300y.T,'--')
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v$ (m/s)")
# plt.title("LW : $t=T_{max}$")
# plt.savefig('Figures/Comp_LW-V_300.eps', format='eps')

# plt.figure()
# plt.plot(X,Res[0,:]-LWV300y.T)
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v-v_{Bruno}$ (m/s)")
# plt.title("LW : $t=T_{max}$")
# plt.savefig('Figures/Diff_LW-V_300.eps', format='eps')


K=4
Res,Vfilm,Sfilm=scheme1D.FD_cauchyII(U0,Nt,Nx,K,rhox,Celx,dt,dx,configuration,mat,frontiere,5)

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


# plt.figure()
# plt.plot(X,Res[0,:])
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v$ (m/s)")
# plt.title("ADER4 : $t=T_{max}$")
# plt.axvline(x=alpha[0],color='gray',linestyle='--')
# # plt.savefig('Figures/ADER4-V_300.eps', format='eps')

# plt.figure()
# plt.plot(X,Res[1,:])
# plt.xlabel("$x$ (m)")
# plt.ylabel("$\sigma$ (Pa)")
# plt.title("ADER4 : $t=T_{max}$")
# plt.axvline(x=alpha[0],color='gray',linestyle='--')
# plt.savefig('Figures/ADER4-S_300.eps', format='eps')


# plt.figure()
# plt.plot(X,Res[0,:])
# plt.plot(X,A4V300y.T,'--')
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v$ (m/s)")
# plt.title("ADER4 : $t=T_{max}$")
# plt.savefig('Figures/Comp_ADER4-V_300.eps', format='eps')

# plt.figure()
# plt.plot(X,Res[0,:]-A4V300y.T)
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v-v_{Bruno}$ (m/s)")
# plt.title("ADER4 : $t=T_{max}$")
# plt.savefig('Figures/Diff_ADER4-V_300.eps', format='eps')
