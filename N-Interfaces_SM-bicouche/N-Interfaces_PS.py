#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 13:44:33 2023

@author: michael
"""

import os
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/N-Interfaces_SM-bicouche/")

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../WAPOD")

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
fs,CFL,x_0=inputReading.source(sourcef)
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
Nt=2**12
Tmax=Nt*dt
tmax=str(Tmax)+"s"
K=4

t=np.linspace(0,Tmax,Nt+1)

sce=np.zeros(Nt+1)
for i in range(Nt+1):
    sce[i]=source.choice_timefct(sourcef,t[i])
    
pFilm=1
Res,Vfilm,Sfilm=scheme1D.FD_sourceVII(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0,sce,configuration,mat,frontiere,pFilm)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
ResD,VfilmD,SfilmD=scheme1D.FD_sourceV(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0,sce,pFilm)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")


energyM=np.zeros(int(Nt/pFilm))
for ti in range(int(Nt/pFilm)):

     energyM[ti]=1/2*np.sum(rhox*Vfilm[ti]*Vfilm[ti])+1/2*np.sum(1/rhox/Celx/Celx*Sfilm[ti]*Sfilm[ti])
energyN=np.zeros(int(Nt/pFilm))
for ti in range(int(Nt/pFilm)):

     energyN[ti]=1/2*np.sum(rhox*VfilmD[ti]*VfilmD[ti])+1/2*np.sum(1/rhox/Celx/Celx*SfilmD[ti]*SfilmD[ti])
plt.figure()
plt.plot
plt.plot(energyM)
plt.plot(energyN)
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
# plt.plot(X,Res[0,:])
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v$ (m/s)")
# plt.title("LW : $t=$"+tmax)
# # plt.savefig('Figures/LW-V_300.eps', format='eps')

# plt.figure()
# plt.plot(X,Res[1,:])
# plt.xlabel("$x$ (m)")
# plt.ylabel("$\sigma$ (Pa)")
# plt.title("LW : $t=$"+tmax)
# plt.savefig('Figures/LW-S_300.eps', format='eps')

nom_fichier = "Vtot_"+sourcef["Frequence"]+"Hz.txt"

# Exportez le tableau en tant que fichier texte
np.savetxt(nom_fichier, Vfilm, fmt='%d', delimiter='\t')

nom_fichier = "Stot_"+sourcef["Frequence"]+"Hz.txt"


# Exportez le tableau en tant que fichier texte
np.savetxt(nom_fichier, Sfilm, fmt='%d', delimiter='\t')


# nameVFH="Veq_"+sourcef["Frequence"]+"Hz.txt"
# VfilmH = np.loadtxt(nameVFH, dtype=int, delimiter='\t')

from matplotlib.animation import FuncAnimation

# Fonction pour mettre à jour le contenu du graphique à chaque image
def update(frame):
    plt.clf()
    plt.plot(X, Vfilm[frame,:])#, label=f'Frame {frame}')
    plt.plot(X, VfilmD[frame,:])
    plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0008,0.0008], aspect='auto', alpha=0.35)
    plt.colorbar()
    # plt.plot(X, VfilmH[frame,:])
    plt.xlabel("x (m)")
    plt.ylabel("v (m/s)")
    plt.ylim([-0.0008,0.0008])
    # plt.legend()


total_frames = Vfilm.shape[0]
fig, ax = plt.subplots()

# Cree l'animation
animation = FuncAnimation(fig, update, frames=total_frames, interval=100)
# Affiche l'animation
plt.show()
nomAnim="Film_NInterfaces_"+sourcef["Frequence"]+"Hz.mp4"
animation.save(nomAnim, writer='ffmpeg', fps=20)


def specfkR(U,dx,dt,Nt,Nx):
    """
    FK spectrum of traces using the numpy.fft functions.
    """
    # Time padding
    zt=2**12
    zx=2**12
    Uzp=np.zeros([zt,zx])
    Uzp[:Nt-1,int(zx/2-Nx/2):int(zx/2+Nx/2)]=U
    #Fourier in time domain 
    FFTx = np.fft.fft(Uzp, axis=1)
    #Fourier in sapce
    FFTt = np.fft.rfft(FFTx, axis=0)
    FFT = np.flip(np.fft.fftshift(FFTt, axes=1), axis=1)
    FFK = np.absolute(FFT)

    # Get the frequency and K vectors
    frqv = np.fft.rfftfreq(zt, dt)
    wavv = np.fft.fftfreq(zx, dx)
    # Centering
    wavv = np.fft.fftshift(wavv)
    return frqv,wavv,FFK

# pFilm=1
# # frqv,wavv,FFK=specfk(Vfilm, dx, dt2p*pFilm, int(Ntp/pFilm)+1, Nx)
# # frqvS,wavvS,FFKS=specfkR(VfilmSS, dx, dt2p*pFilm, int(Ntp/pFilm)+1, Nx)
# frqvG,wavvG,FFKG=specfkR(Vfilm[::,::], dx, dt*pFilm, int(Nt/pFilm)+1, int((Nx)))
# # plt.figure()
# # # plt.pcolor(wavv[::]*cmet*tau,frqv[::2]*tau,np.log(FFK[::2,::]),vmin=-12, vmax=2)
# # # plt.ylim([0,2])
# # # plt.xlim([0,12])
# # plt.colorbar()
# # plt.xlabel("k (1/m)")
# # plt.ylabel("f (Hz)")
# # plt.tight_layout()

# plt.figure()
# plt.pcolor(wavvG[::],frqvG[::],np.log(FFKG[::,::]),vmin=-5, vmax=4)
# # plt.ylim([0,8*np.pi])
# # plt.xlim([-0,8*np.pi])
# plt.colorbar()
# plt.xlabel("kc^*\tau")
# plt.ylabel("f*\tau")
# plt.tight_layout()


# frqvD,wavvD,FFKD=specfkR(VfilmD[::,::], dx, dt*pFilm, int(Nt/pFilm)+1, int((Nx)))
# # plt.figure()
# # # plt.pcolor(wavv[::]*cmet*tau,frqv[::2]*tau,np.log(FFK[::2,::]),vmin=-12, vmax=2)
# # # plt.ylim([0,2])
# # # plt.xlim([0,12])
# # plt.colorbar()
# # plt.xlabel("k (1/m)")
# # plt.ylabel("f (Hz)")
# # plt.tight_layout()

# plt.figure()
# plt.pcolor(wavvD[::],frqvD[::],np.log(FFKD[::,::]),vmin=-5, vmax=4)
# # plt.ylim([0,8*np.pi])
# # plt.xlim([-0,8*np.pi])
# plt.colorbar()
# plt.xlabel("kc^*\tau")
# plt.ylabel("f*\tau")
# plt.tight_layout()
# plt.figure()
# ###############################

# THV0y=np.zeros(Nx)
# THV100y=np.zeros(Nx)
# THS0y=np.zeros(Nx)
# THS100y=np.zeros(Nx)
# R01=(rho[1]*cm[1]-rho[0]*cm[0])/(rho[1]*cm[1]+rho[0]*cm[0])
# T01=2*rho[0]*cm[1]/(rho[1]*cm[1]+rho[0]*cm[0])
# Nxalpha=int(alpha[0]/dx)
# for x in range(Nxalpha):
#     xabs=(x*dx)
#     THV0y[x]=1*cm[0]*source.choice_timefct(sourcef,-xabs/cm[0])+1/cm[0]*R01*source.choice_timefct(sourcef,xabs/cm[0]-2*alpha[0]/cm[0])
#     THS0y[x]=-rho[0]*source.choice_timefct(sourcef,-xabs/cm[0])+rho[0]*R01*source.choice_timefct(sourcef,xabs/cm[0]-2*alpha[0]/cm[0])
#     THV100y[x]=1/cm[0]*source.choice_timefct(sourcef,Tmax-xabs/cm[0])+1/cm[0]*R01*source.choice_timefct(sourcef,Tmax+xabs/cm[0]-2*alpha[0]/cm[0])
#     THS100y[x]=-rho[0]*source.choice_timefct(sourcef,Tmax-xabs/cm[0])+rho[0]*R01*source.choice_timefct(sourcef,Tmax+xabs/cm[0]-2*alpha[0]/cm[0])
# for x in range(Nxalpha,Nx):
#     xabs=(x*dx)
#     THV0y[x]=1/cm[1]*T01*source.choice_timefct(sourcef,-xabs/cm[1]+(1/cm[1]-1/cm[0])*alpha[0])
#     THS0y[x]=-rho[1]*T01*source.choice_timefct(sourcef,-xabs/cm[1]+(1/cm[1]-1/cm[0])*alpha[0])
#     THV100y[x]=1/cm[1]*T01*source.choice_timefct(sourcef,Tmax-xabs/cm[1]+(1/cm[1]-1/cm[0])*alpha[0])
#     THS100y[x]=-rho[1]*T01*source.choice_timefct(sourcef,Tmax-xabs/cm[1]+(1/cm[1]-1/cm[0])*alpha[0])



# plt.figure()
# plt.plot(X,THV0y)
# plt.plot(X,U0[0,:])
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v_{th}$ (m/s)")
# plt.title("LW : $t=T_{max}$")

# plt.figure()
# plt.plot(X,THS0y)
# plt.plot(X,U0[1,:])
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v_{th}$ (m/s)")
# plt.title("LW : $t=T_{max}$")

# plt.figure()
# plt.plot(X,THV100y)
# plt.plot(X,Res[0,:])
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v_{th}$ (m/s)")
# plt.title("LW : $t=T_{max}$")

# plt.figure()
# plt.plot(X,THS100y)
# plt.plot(X,Res[1,:])
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v_{th}$ (m/s)")
# plt.title("LW : $t=T_{max}$")
    
# with open("DataBruno/LW-S-Iter0.txt", "r") as fichier:
#     lignes = fichier.readlines()
# colonne1 = []
# colonne2 = []
# for ligne in lignes:
#     mots = ligne.split()
#     if len(mots) >= 2:
#         try:
#             valeur_colonne1 = float(mots[0])
#             valeur_colonne2 = float(mots[1])
#             colonne1.append(valeur_colonne1)
#             colonne2.append(valeur_colonne2)
#         except ValueError:
#             pass
# LWS0x = np.array(colonne1)
# LWS0y = np.array(colonne2)        


# with open("DataBruno/LW-V-Iter0.txt", "r") as fichier:
#     lignes = fichier.readlines()
# colonne1 = []
# colonne2 = []
# for ligne in lignes:
#     mots = ligne.split()
#     if len(mots) >= 2:
#         try:
#             valeur_colonne1 = float(mots[0])
#             valeur_colonne2 = float(mots[1])
#             colonne1.append(valeur_colonne1)
#             colonne2.append(valeur_colonne2)
#         except ValueError:
#             pass
# LWV0x = np.array(colonne1)
# LWV0y = np.array(colonne2)

# with open("DataBruno/LW-V-Iter300.txt", "r") as fichier:
#     lignes = fichier.readlines()
# colonne1 = []
# colonne2 = []
# for ligne in lignes:
#     mots = ligne.split()
#     if len(mots) >= 2:
#         try:
#             valeur_colonne1 = float(mots[0])
#             valeur_colonne2 = float(mots[1])
#             colonne1.append(valeur_colonne1)
#             colonne2.append(valeur_colonne2)
#         except ValueError:
#             pass
# LWV300x = np.array(colonne1)
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
# Res=scheme1D.FD_cauchy(U0,Nt,Nx,K,A,dt,dx,"no")

# plt.figure()
# plt.plot(X,U0[0,:])
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v$ (m/s)")
# plt.title("ADER4 : $t=$0 s")
# plt.savefig('Figures/ADER4-V_0.eps', format='eps')

# plt.figure()
# plt.plot(X,U0[1,:])
# plt.xlabel("$x$ (m)")
# plt.ylabel("$\sigma$ (Pa)")
# plt.title("ADER4 : $t=$0 s")
# plt.savefig('Figures/ADER4-S_0.eps', format='eps')


# plt.figure()
# plt.plot(X,Res[0,:])
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v$ (m/s)")
# plt.title("ADER4 : $t=T_{max}$")
# plt.savefig('Figures/ADER4-V_300.eps', format='eps')

# plt.figure()
# plt.plot(X,Res[1,:])
# plt.xlabel("$x$ (m)")
# plt.ylabel("$\sigma$ (Pa)")
# plt.title("ADER4 : $t=T_{max}$")
# plt.savefig('Figures/ADER4-S_300.eps', format='eps')

# with open("DataBruno/ADER4-V-Iter300.txt", "r") as fichier:
#     lignes = fichier.readlines()
# colonne1 = []
# colonne2 = []
# for ligne in lignes:
#     mots = ligne.split()
#     if len(mots) >= 2:
#         try:
#             valeur_colonne1 = float(mots[0])
#             valeur_colonne2 = float(mots[1])
#             colonne1.append(valeur_colonne1)
#             colonne2.append(valeur_colonne2)
#         except ValueError:
#             pass
# A4V300x = np.array(colonne1)
# A4V300y = np.array(colonne2)

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
