#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 13:44:33 2023

@author: michael
"""

import os
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/N-Interfaces_SM_Cracks/")

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
Nt=4000
Tmax=Nt*dt
tmax=str(Tmax)+"s"
K=4

t=np.linspace(0,Tmax,Nt+1)

sce=np.zeros(Nt+1)
for i in range(Nt+1):
    sce[i]=source.choice_timefct(sourcef,t[i])
import time
time0=time.time()
Res,Vfilm,Sfilm=scheme1D.FD_sourceVII(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0,sce,configuration,mat,frontiere,5)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
timef=time.time()-time0
print(timef)  
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
    plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0008,0.0008], aspect='auto', alpha=0.35)
    plt.colorbar()
    # plt.plot(X, VfilmH[frame,:])
    plt.xlabel("x (m)")
    plt.ylabel("v (m/s)")
    plt.ylim([-0.0008,0.0008])
    for i in range(alpha.size):
        plt.axvline(x=alpha[i],color='gray',linestyle='--')
    # plt.legend()


total_frames = Vfilm.shape[0]
fig, ax = plt.subplots()

# Cree l'animation
animation = FuncAnimation(fig, update, frames=total_frames, interval=100)
# Affiche l'animation
plt.show()
nomAnim="Film_NInterfaces_"+sourcef["Frequence"]+"Hz.mp4"
animation.save(nomAnim, writer='ffmpeg', fps=20)


import homogenize
l=alpha[0]-xmin
matH=inputReading.readfile('EqH/Milieu.txt')
sourceH=inputReading.readfile('Source.txt')   
fs,CFL,x_0=inputReading.source(sourceH)
Kn=float(frontiere["Kn_0"])
Mn=float(frontiere["Mn_0"])
cmH,rhoH,fCFLh=homogenize.hPropCracks(cm[0],rho[0],Kn,Mn,l)
sourceH["Force"]=str(rhox[int(x_0/dx)]/rhoH)
matH["Rho_0 "]=str(rhoH)
matH["Cel_0 "]=str(cmH)
sourceH["CFL"]=str(CFL)
sce=np.zeros(Nt+1)
for i in range(Nt+1):
    sce[i]=source.choice_timefct(sourceH,t[i])

################ Parameters of the simulation

CFLH=fCFLh*CFL
rhohx,Celhx=inputReading.affectMaterialsH(rhoH,cmH,xmin,xmax,Nx,dx)
dt=np.min(scheme1D.timeStep(CFLH,dx,cmH))
# lambdaOnde=cm/fs
# xsmin=cm*tshift-lambdaOnde
# xsmax=cm*tshift
#%

#% Initial conditions
U0h=scheme1D.initPointSourceProblem(configuration,frontiere,matH,sourceH)

ResH,Vfilmh,Sfilmh=scheme1D.FD_sourceV(U0h,Nt,Nx,K,rhohx,Celhx,dt,dx,x_0,sce,5)

nom_fichier = "Vtot_"+sourcef["Frequence"]+"Hz.txt"

# Exportez le tableau en tant que fichier texte
np.savetxt(nom_fichier, Vfilmh, fmt='%d', delimiter='\t')

nom_fichier = "Stot_"+sourcef["Frequence"]+"Hz.txt"


# Exportez le tableau en tant que fichier texte
np.savetxt(nom_fichier, Sfilmh, fmt='%d', delimiter='\t')


# nameVFH="Veq_"+sourcef["Frequence"]+"Hz.txt"
# VfilmH = np.loadtxt(nameVFH, dtype=int, delimiter='\t')


# Fonction pour mettre à jour le contenu du graphique à chaque image
def updateH(frame):
    plt.clf()
    plt.plot(X, Vfilmh[frame,:])#, label=f'Frame {frame}')
    plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0008,0.0008], aspect='auto', alpha=0.35)
    plt.colorbar()
    # plt.plot(X, VfilmH[frame,:])
    plt.xlabel("x (m)")
    plt.ylabel("v (m/s)")
    plt.ylim([-0.0008,0.0008])
    # plt.legend()


total_frames = Vfilmh.shape[0]
fig, ax = plt.subplots()

# Cree l'animation
animationH = FuncAnimation(fig, updateH, frames=total_frames, interval=100)
# Affiche l'animation
plt.show()
nomAnim="Film_NInterfaces_"+sourcef["Frequence"]+"Hz_eqH.mp4"
animationH.save(nomAnim, writer='ffmpeg', fps=20)

def updateComp(frame):
    plt.clf()
    plt.plot(X, Vfilm[frame,:])#, label=f'Frame {frame}')
    plt.plot(X, Vfilmh[frame,:],"--")#, label=f'Frame {frame}')
    plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0008,0.0008], aspect='auto', alpha=0.35)
    plt.colorbar()
    # plt.plot(X, VfilmH[frame,:])
    plt.xlabel("x (m)")
    plt.ylabel("v (m/s)")
    plt.ylim([-0.0008,0.0008])
    # plt.legend()


total_frames = Vfilmh.shape[0]
fig, ax = plt.subplots()

# Cree l'animation
animationComp = FuncAnimation(fig, updateComp, frames=total_frames, interval=100)
plt.show()
nomAnim="Film_NInterfaces_"+sourcef["Frequence"]+"Hz_Comp.mp4"
animationComp.save(nomAnim, writer='ffmpeg', fps=20)


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
