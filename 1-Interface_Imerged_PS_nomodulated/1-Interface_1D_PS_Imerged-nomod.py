#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:50:43 2023

@author: michael
"""

import os
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/1-Interface_Imerged_PS_nomodulated/")

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../WAPOD/")

import scheme1D as scheme1D
import source
import inputReading



################# Input files 
#%
configuration=inputReading.readfile('Demarrer.txt')
mat=inputReading.readfile('Milieu.txt')
frontiere=inputReading.readfile('Frontiere.txt')     
sourcef=inputReading.readfile('Source.txt')   
sourceRev=inputReading.readfile('SourceRev.txt')
modulation=inputReading.readfile('Modulation.txt')
nomodulation=inputReading.readfile('NoModulation.txt')
################ Parameters of the simulation   
X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(configuration)
fs,CFL,x_0=inputReading.source(sourcef)
fs,CFL,x_0R=inputReading.source(sourceRev)
Nmat,rho,cm=inputReading.material(mat)
alpha=inputReading.frontiere(frontiere)
rhox,Celx=inputReading.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
dt=min(scheme1D.timeStep(CFL,dx,cm))
# lambdaOnde=cm/fs
# xsmin=cm*tshift-lambdaOnde
# xsmax=cm*tshift
#%

#% Initial conditions
U0=scheme1D.initPointSourceProblem(configuration,frontiere,mat,sourcef)

###### Mise en oeuvre
Nt=2000
Tmax=Nt*dt
tmax=str(Tmax)+"s"
K=4

t=np.linspace(0,Tmax,Nt+1)

sce=np.zeros(Nt+1)
for i in range(Nt+1):
    sce[i]=source.choice_timefct(sourcef,t[i])
    
sceR=np.zeros(Nt+1)
for i in range(Nt+1):
    sceR[i]=source.choice_timefct(sourceRev,t[i])
    
    
from matplotlib.animation import FuncAnimation

pFilm=1
############### MODULATION DIRECT PROPAGATION ##########################################
import time
time0=time.time()
Res0,Vfilm0,Sfilm0=scheme1D.FD_sourceVII(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0,sce,configuration,mat,frontiere,pFilm)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
timeM=time.time()-time0
print(timeM)  
# Fonction pour mettre à jour le contenu du graphique à chaque image
def update0(frame):
    plt.clf()
    plt.plot(X, Vfilm0[frame,:])#, label=f'Frame {frame}')
    plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0008,0.0008], aspect='auto', alpha=0.35)
    plt.axvline(x=alpha[0],color='gray',linestyle='--')
    plt.colorbar()
    # plt.plot(X, VfilmH[frame,:])
    plt.xlabel("x (m)")
    plt.ylabel("v (m/s)")
    plt.ylim([-0.0008,0.0008])
    plt.title(t[frame]*2)
    # plt.legend()


total_frames = Vfilm0.shape[0]
fig0, ax0 = plt.subplots()

# Cree l'animation
animation0 = FuncAnimation(fig0, update0, frames=total_frames, interval=100)
# Affiche l'animation
plt.show()
# nomAnim="Film_SM_"+sourcef["Frequence"]+"Hz.mp4"
# animation.save(nomAnim, writer='ffmpeg', fps=20)

################ REVERSE PROPAGATION ####################################
time0=time.time()
ResR0,VfilmR0,SfilmR0=scheme1D.FD_sourceVII(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0R,sceR,configuration,mat,frontiere,pFilm)
timeR=time.time()-time0
print(timeR)  
def updateR0(frame):
    plt.clf()
    plt.plot(X, VfilmR0[frame,:])#, label=f'Frame {frame}')
    plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0008,0.0008], aspect='auto', alpha=0.35)
    plt.axvline(x=alpha[0],color='gray',linestyle='--')
    plt.colorbar()
    # plt.plot(X, VfilmH[frame,:])
    plt.xlabel("x (m)")
    plt.ylabel("v (m/s)")
    plt.ylim([-0.0008,0.0008])
    plt.title(t[frame]*2)
    # plt.legend()

total_frames = VfilmR0.shape[0]
figR0, axR0 = plt.subplots()

# Cree l'animation
animationR0 = FuncAnimation(figR0, updateR0, frames=total_frames, interval=100)
# Affiche l'animation
plt.show()
# nomAnimR="Film_SMR_"+sourcef["Frequence"]+"Hz.mp4"
# animationR.save(nomAnimR, writer='ffmpeg', fps=20)


###################### NO MODULATION ##########################################################################
# time0=time.time()
# ResNM,VfilmNM,SfilmNM=scheme1D.FD_sourceVII(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0,sce,configuration,mat,frontiere,2)
# timeNM=time.time()-time0
# print(timeNM)  

# # Fonction pour mettre à jour le contenu du graphique à chaque image
# def updateNM(frame):
#     plt.clf()
#     plt.plot(X, VfilmNM[frame,:])#, label=f'Frame {frame}')
#     plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0008,0.0008], aspect='auto', alpha=0.35)
#     plt.axvline(x=alpha[0],color='gray',linestyle='--')
#     plt.colorbar()
#     # plt.plot(X, VfilmH[frame,:])
#     plt.xlabel("x (m)")
#     plt.ylabel("v (m/s)")
#     plt.ylim([-0.0008,0.0008])
#     # plt.legend()


# total_frames = VfilmNM.shape[0]
# figNM, axNM = plt.subplots()

# # Cree l'animation
# animationNM = FuncAnimation(figNM, updateNM, frames=total_frames, interval=100)
# # Affiche l'animation
# plt.show()
# nomAnimNM="Film_SMNM_"+sourcef["Frequence"]+"Hz.mp4"
# animationNM.save(nomAnimNM, writer='ffmpeg', fps=20)


#######################################
# plt.figure()
# plt.clf()
# plt.plot(Vfilm[:,int(x_0R)])
# plt.plot(VfilmNM[:,int(x_0R)])
# plt.xlabel("t (s)")

plt.figure()
plt.clf()
plt.plot(VfilmR0[:,int(x_0)],label="Xsource<-Xrecepteur")
plt.plot(Vfilm0[:,int(x_0R)],"--",label="Xsource->Xrecepteur")
plt.xlabel("t (s)")
plt.ylabel("v (m/s)")
plt.legend(loc= 'upper right')

plt.figure()
plt.clf()
plt.plot(VfilmNM[:,int(x_0R)], label="Ancienne version non-module")
plt.plot(Vfilm0[:,int(x_0R)],"--", label= "Modulation nulle")
plt.xlabel("t (s)")
plt.ylabel("v (m/s)")
plt.legend(loc= 'upper right')


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

# plt.figure()
# plt.plot(X,U0[0,:])
# plt.plot(X,THV0y,'--')
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v_{th}$ (m/s)")
# plt.title("LW : $t=T_{max}$")

# plt.figure()
# plt.plot(X,U0[1,:])
# plt.plot(X,THS0y,'--')
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v_{th}$ (m/s)")
# plt.title("LW : $t=T_{max}$")

# plt.figure()
# plt.plot(X,Res[0,:])
# plt.plot(X,THVTMy,'--')
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v_{th}$ (m/s)")
# plt.title("LW : $t=T_{max}$")

# plt.figure()
# plt.plot(X,Res[1,:])
# plt.plot(X,THSTMy,'--')
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v_{th}$ (m/s)")
# plt.title("LW : $t=T_{max}$")


# K=2
# for k in range(ndx):
#     k0=k0/2
#     Res=scheme1D.FD_cauchyII(U0,Nt,Nx,K,rhox,Celx,dt,k0*dx,configuration,mat,frontiere,"no")
    
# K=4
# k0=2**((ndx-1)/2)
# for k in range(ndx):
#     k0=k0/2
#     Res=scheme1D.FD_cauchyII(U0,Nt,Nx,K,rhox,Celx,dt,k0*dx,configuration,mat,frontiere,"no")
#     err[2,k]=errL2(Res[0,:], THVTMy)
#     err[3,k]=errL2(Res[1,:], THSTMy)
    

    

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


# K=4
# Res,Vfilm,Sfilm=scheme1D.FD_cauchyII(U0,Nt,Nx,K,rhox,Celx,dt,dx,configuration,mat,frontiere,5)

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
