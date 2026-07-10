#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 14:24:41 2023

@author: michael
"""


import os
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/1-Interface_PointSource/")

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
fs,CFL,x0=inputReading.source(sourcef)
Nmat,rho,cm=inputReading.material(mat)
alpha=inputReading.frontiere(frontiere)
rhox,Celx=inputReading.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
dt=min(scheme1D.timeStep(CFL,dx,cm))


#% Initial conditions
U0=scheme1D.initPointSourceProblem(configuration,frontiere,mat,sourcef)

###### Mise en oeuvre
Nt=250
Tmax=(Nt)*dt
tmax=str(Tmax)+"s"
K=4

t=np.linspace(0,Tmax,Nt+1)

sce=np.zeros(Nt+1)
for i in range(Nt+1):
    sce[i]=source.choice_timefct(sourcef,t[i])

Res,Vfilm,Sfilm=scheme1D.FD_sourceV(U0,Nt,Nx,K,rhox,Celx,dt,dx,x0,sce,5)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")

plt.figure()
plt.plot(X,U0[0,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$v$ (m/s)")
plt.title("LW : $t=$0 s")
# plt.savefig('Figures/LW-V_0.eps', format='eps')

plt.figure()
plt.plot(X,U0[1,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$\sigma$ (Pa)")
plt.title("LW : $t=$0 s")
# plt.savefig('Figures/LW-S_0.eps', format='eps')

plt.figure()
plt.plot(X,Res[0,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$v$ (m/s)")
plt.title("LW : $t=$"+tmax)
# plt.savefig('Figures/LW-V_300.eps', format='eps')

plt.figure()
plt.plot(X,Res[1,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$\sigma$ (Pa)")
plt.title("LW : $t=$"+tmax)
# plt.savefig('Figures/LW-S_300.eps', format='eps')


# ###############################
    
THV0y=np.zeros(Nx)
THVTMy=np.zeros(Nx)
THS0y=np.zeros(Nx)
THSTMy=np.zeros(Nx)
R01=(rho[1]*cm[1]-rho[0]*cm[0])/(rho[1]*cm[1]+rho[0]*cm[0])
T01=2*rho[0]*cm[1]/(rho[1]*cm[1]+rho[0]*cm[0])
Nxalpha=int(alpha[0]/dx)
for x in range(Nxalpha):
    f0p=source.choice_timefct(sourcef,-(x*dx-x0)/cm[0])
    f0m=+source.choice_timefct(sourcef,(x*dx-x0)/cm[0])
    fr=source.choice_timefct(sourcef,(x*dx+x0)/cm[0]-2*alpha[0]/cm[0])
    THV0y[x]=1/2/cm[0]*(f0p+f0m)-R01/cm[0]/2*fr
    THS0y[x]=(-rho[0]/2*(f0p+f0m)*np.sign(x*dx-x0)-R01*rho[0]/2*fr)#
    fTp=source.choice_timefct(sourcef,Tmax-(x*dx-x0)/cm[0])
    fTm=source.choice_timefct(sourcef,Tmax+(x*dx-x0)/cm[0])
    frT=source.choice_timefct(sourcef,Tmax+(x*dx+x0)/cm[0]-2*alpha[0]/cm[0])
    THVTMy[x]=1/2/cm[0]*(fTp+fTm)-R01/cm[0]/2*frT
    THSTMy[x]=(-rho[0]/2*(fTp+fTm)*np.sign(x*dx-x0)-R01*rho[0]/2*frT)#
for x in range(Nxalpha,Nx):
    ft=source.choice_timefct(sourcef,-(-x0)/cm[0]-x*dx/cm[1]+(1/cm[1]-1/cm[0])*alpha[0])
    ftT=source.choice_timefct(sourcef,Tmax-(-x0)/cm[0]-x*dx/cm[1]+(1/cm[1]-1/cm[0])*alpha[0])
    THV0y[x]=1/2/cm[1]*T01*ft
    THS0y[x]=(-rho[1]/2*T01*ft)*np.sign(x*dx-x0)
    THVTMy[x]=1/cm[1]/2*T01*ftT
    THSTMy[x]=(-rho[1]/2*T01*ftT)*np.sign(x*dx-x0)


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

plt.figure()
plt.plot(X,THVTMy)
plt.plot(X,Res[0,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$v_{th}$ (m/s)")
plt.title("LW : $t=T_{max}$")

plt.figure()
plt.plot(X,THSTMy)
plt.plot(X,Res[1,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$v_{th}$ (m/s)")
plt.title("LW : $t=T_{max}$")
    
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
