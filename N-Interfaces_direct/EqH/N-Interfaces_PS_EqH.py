#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 13:44:33 2023

@author: michael
"""

import os
os.chdir("/home/michael/Documents/Post-doc IMI/Simulations/1D/N-Interfaces_direct/EqH/")

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../../../1D/WAPOD/")

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
dt=0.000475
# lambdaOnde=cm/fs
# xsmin=cm*tshift-lambdaOnde
# xsmax=cm*tshift
#%

#% Initial conditions
U0=scheme1D.initPointSourceProblem(configuration,frontiere,mat,sourcef)

###### Mise en oeuvre
Nt=500
Tmax=Nt*dt
tmax=str(Tmax)+"s"
K=4

t=np.linspace(0,Tmax,Nt+1)

sce=np.zeros(Nt+1)
for i in range(Nt+1):
    sce[i]=source.choice_timefct(sourcef,t[i])

ResH,VfilmH,SfilmH=scheme1D.FD_sourceVII(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0,sce,configuration,mat,frontiere,5)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")

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

nom_fichier = "../Veq_"+sourcef["Frequence"]+"Hz.txt"

# Exportez le tableau en tant que fichier texte
np.savetxt(nom_fichier, VfilmH, fmt='%d', delimiter='\t')

nom_fichier = "../Seq_"+sourcef["Frequence"]+"Hz.txt"

# Exportez le tableau en tant que fichier texte
np.savetxt(nom_fichier, SfilmH, fmt='%d', delimiter='\t')

from matplotlib.animation import FuncAnimation

# Fonction pour mettre à jour le contenu du graphique à chaque image
def update(frame):
    plt.clf()
    plt.plot(X, VfilmH[frame,:])#, label=f'Frame {frame}')
    plt.xlabel("x (m)")
    plt.xlabel("v (m/s)")


total_frames = VfilmH.shape[0]
fig, ax = plt.subplots()

# Cree l'animation
animation = FuncAnimation(fig, update, frames=total_frames, interval=100)

# Affiche l'animation
plt.show()
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
