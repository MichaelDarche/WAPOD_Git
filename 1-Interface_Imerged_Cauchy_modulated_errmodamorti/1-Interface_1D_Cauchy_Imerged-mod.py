#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:50:43 2023

@author: michael
"""

import os
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/1-Interface_Imerged_Cauchy_modulated_errmodamorti")
#os.chdir("/home/darche/Documents/ModulatedMetamaterials/WAPRO_TD/1-Interface_Imerged_Cauchy_modulated")

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../WAPOD/")
sys.path.append("/home/darche/Documents/ModulatedMetamaterials/WAPRO_TD/WAPOD/")

import scheme1D_TDold as scheme1D
import source
import inputReading

import modulation as modulateP
import analytics
NOde=20001

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
################# Input files 
#%
configuration=inputReading.readfile('Demarrer.txt')
mat=inputReading.readfile('Milieu.txt')
frontiere=inputReading.readfile('Frontiere.txt')     
sourcef=inputReading.readfile('Source.txt')   
# sourceRev=inputReading.readfile('SourceRev.txt')
modulation=inputReading.readfile('Modulation.txt')
nomodulation=inputReading.readfile('NoModulation.txt')
################ Parameters of the simulation   
X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(configuration)
fs,CFL,x_0=inputReading.source(sourcef)
# fsR,CFLR,x_0R=inputReading.source(sourceRev)
Nmat,rho,cm=inputReading.material(mat)
alpha=inputReading.frontiere(frontiere)
rhox,Celx=inputReading.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
dt=min(scheme1D.timeStep(CFL,dx,cm))
# lambdaOnde=cm/fs
# xsmin=cm*tshift-lambdaOnde
# xsmax=cm*tshift
#%

pFilm=1


#% Initial conditions
U0=scheme1D.initCauchyProblem(configuration,frontiere,mat,sourcef)
#U0=scheme1D.initPointSourceProblem(configuration,frontiere,mat,sourcef)
###### Mise en oeuvre
# Nt=600
# Tmax=Nt*dt
# tmax=str(Tmax)+"s"
K=4

# t=np.linspace(0,Tmax,Nt+1)

# sce=np.zeros(Nt+1)
# for i in range(Nt+1):
#     sce[i]=source.choice_timefct(sourcef,t[i])
    
    
    
Ntfig=int(np.floor(0.045/dt))
# Ntfig=600
Tmax=Ntfig*dt
tmax=str(Tmax)+"s"
t=np.linspace(0,Tmax,Ntfig+1)
sce=np.zeros(Ntfig+1)
for i in range(Ntfig+1):
    sce[i]=source.choice_timefct(sourcef,t[i])
# THV0y=np.zeros(Nx)
# THVTMy=np.zeros(Nx)
# THS0y=np.zeros(Nx)
# THSTMy=np.zeros(Nx)
# Nxalpha=int(alpha[0]/dx)

SolV,SolS=analytics.analytics_1IMod(Tmax,X,sourcef,modulation,frontiere,mat,Ntfig+1,dt,Ntfig+1)
# analytics_Cauchy_1IModK

frontierep=inputReading.readfile('FrontierePlus.txt')  
SolVp,SolSp=analytics.analytics_1IMod(Tmax,X,sourcef,nomodulation,frontierep,mat,Ntfig+1,dt,Ntfig+1)

frontierem=inputReading.readfile('FrontiereMoins.txt')  
SolVm,SolSm=analytics.analytics_1IMod(Tmax,X,sourcef,nomodulation,frontierem,mat,Ntfig+1,dt,Ntfig+1)

# Res=np.zeros([2,Nx])
# Res[0,:]=np.loadtxt(nSve, delimiter='\t')
# Res[1,:]=np.loadtxt(nSveS, delimiter='\t')
Res,Vfilm,Sfilm,CSVM=scheme1D.FD_cauchyII(U0,Ntfig,Nx,K,rhox,Celx,dt,dx,configuration,mat,frontiere,modulation,pFilm,rCSVM="yes")
#Res,Vfilm,Sfilm,CSVM,traces=scheme1D.FD_sourceVII(U0,Ntfig,Nx,K,rhox,Celx,dt,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm,rCSVM="yes",rtraces="yes")

plt.figure()
plt.axvline(x=alpha[0],color='black')
plt.plot(X, Res[0,:],linestyle=None,marker='o', label="numerics")
plt.plot(X,SolV,linewidth=2,label='analytics')
plt.xlabel(r'$x$ (m)',fontsize=16)
plt.ylabel(r'$v$ (m/s)',fontsize=16)
plt.legend(loc="lower left",fontsize=16)
plt.tight_layout()
# plt.figure()
# plt.plot(X,SolS)
# plt.plot(X,Res[1,:],'*')

plt.figure()
plt.axvline(x=alpha[0],color='black')
plt.plot(X, Res[0,:],linestyle=None,marker=None, label="numerics")
plt.plot(X,SolV,linewidth=2,label='analytics')
plt.plot(X,SolVp,"--",linewidth=3,label='static $$')
plt.plot(X,SolVm,"--",linewidth=3,label='static $$')
plt.xlabel(r'$x$ (m)',fontsize=16)
plt.ylabel(r'$v$ (m/s)',fontsize=16)
plt.tight_layout()

def cosine_taper(signal, alpha=0.00,zp=18):
    zeropadding=2**zp
    Size=signal.shape[0]
    M=np.floor((Size*alpha)/2+0.5);
    tapered=np.zeros(zeropadding)
    for j in range(int(Size)):
        if j<=M+1:
            tapered[j]=signal[j] * (0.5 * ( 1-np.cos(j*np.pi/(M+1))))
        elif j<Size - M-1:
            tapered[j]=signal[j]
        elif j<=Size:
            tapered[j]=signal[j]* (0.5 * (1-np.cos((Size-j)*np.pi/(M+1))))
    return tapered

N=2**19
sigcor = cosine_taper(sce)
fftSig=1/Ntfig*np.fft.fftshift(np.fft.fft(sigcor,n=N))
freq = np.fft.fftshift(np.fft.fftfreq(N,d=dt))
plt.figure()
plt.plot(np.abs(freq),np.abs(fftSig),linewidth=3)
plt.xlim([0,300])
plt.ylim([0,0.125])
plt.xlabel(r'$f$ (Hz)',fontsize=16)
# plt.ylabel(r'$\overline{S}$',fontsize=16)

plt.ylabel(r'$\hat{S}$',fontsize=16)
plt.axvline(x=fs,color='black',linestyle='--')
plt.tight_layout()



if float(frontiere["Kn_0"])==0:
    K0=np.inf
else:
    K0=float(frontiere["Kn_0"])
frK=float(modulation["Freq"])
delta=float(modulation["DeltaF"])

M0=float(frontiere["Mn_0"])
frM=float(modulation["Freq"])
deltaM=float(modulation["DeltaM"])

modulb=modulateP.fonctionModM(t,alpha[0],modulation,frontiere,0)
modulM=M0*(1+modulb)
dmodulb=np.zeros(int(Ntfig/pFilm))
# for i in range(int(Nt/pFilm)):
#     modulb[i]=M0*(1+deltaM*np.sin(2*np.pi*frM*t[pFilm*i]))
#     dmodulb[i]=2*np.pi*frM*M0*(deltaM*np.cos(2*np.pi*frM*t[pFilm*i]))
    
modula=modulateP.fonctionModK(t,alpha[0],modulation,frontiere,0)
if float(frontiere["Kn_0"])==0:
    modulK=0*modula
else:
    modulK=1/K0*(1+modula)

dmodula=np.zeros(int(Ntfig/pFilm))

energyV=np.zeros(int(Ntfig/pFilm))
energyI=np.zeros(int(Ntfig/pFilm))
for ti in range(int(Ntfig/pFilm)):
    energyV[ti]=dx/2*np.sum(rhox*Vfilm[ti]*Vfilm[ti])+dx/2*np.sum(1/rhox/Celx/Celx*Sfilm[ti]*Sfilm[ti])
    energyI[ti]=1/2*modulM[ti]*CSVM[ti,2]*CSVM[ti,2]+1/2*modulK[ti]*CSVM[ti,3]*CSVM[ti,3]     
    
energyM=energyV+energyI
# plt.figure()
# plt.plot(t[:-1],energyV,label=r'$\mathcal{E}_b$',linewidth=3)
# plt.plot(t[:-1],energyI,label=r'$\mathcal{E}_i$',linewidth=3)
# plt.plot(t[:-1],energyM,label=r'$\mathcal{E}_m$',linewidth=3)
# plt.legend(loc="lower right",fontsize=16)
# plt.xlabel(r'$t$ (s)',fontsize=16)
# plt.ylabel(r'$\mathcal{E}$ (J)',fontsize=16)
# plt.imshow([1+modula], cmap='gray', extent=[0, t[-1], 0,max(energyM)], aspect='auto', alpha=0.35)
# plt.colorbar().set_label(label=r'$1+\varepsilon~\sin (\Omega t)$',size=16)
# plt.xlim([0,0.06])
# plt.tight_layout()
# plt.figure()


# with open("Compare.txt", "r") as fichier:
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
# plt.plot(X,SolV)
# plt.plot(X,Res[0,:],'*')
# plt.plot(LWV300x[1601:3201],LWV300y[1601:3201],"o")
# plt.figure()
# plt.plot(X,SolV)
# plt.plot(X,Res[0,:],'*')
# plt.plot(LWV300x[0:1600],LWV300y[:1600],"o")
# plt.figure()
# plt.plot(X,SolS,'*')
# plt.plot(X,Res[1,:])

# sceR=np.zeros(Nt+1)
# for i in range(Nt+1):
#     sceR[i]=source.choice_timefct(sourceRev,t[i])
    
    
# from matplotlib.animation import FuncAnimation

# pFilm=1
############### MODULATION DIRECT PROPAGATION ##########################################
# import time
# time0=time.time()
# # Res,Vfilm,Sfilm,CSVM=scheme1D.FD_cauchyII(U0,Nt,Nx,K,rhox,Celx,dt,dx,configuration,mat,frontiere,modulation,pFilm,rCSVM="yes")
# #U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
# timeM=time.time()-time0
# # maxV=np.max([np.max(Vfilm),-np.min(Vfilm)])
# print(timeM)  
# Fonction pour mettre à jour le contenu du graphique à chaque image
# def update(frame):
#     plt.clf()xlabel, kwargs)
#     plt.plot(X, Vfilm[frame,:])#, label=f'Frame {frame}')
#     plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0008,0.0008], aspect='auto', alpha=0.35)
#     for i in range(alpha.size-1):
#         if modul.fonctionModK(pFilm*dt*frame,alpha[i],modulation,frontiere,i)<0.:
#             plt.axvline(x=alpha[i],color='red')#,linestyle='--')
#         else:
#             plt.axvline(x=alpha[i],color='lime')#,linestyle='--')
#     plt.colorbar()
#     # plt.plot(X, VfilmH[frame,:])
#     plt.xlabel("x (m)")
#     plt.ylabel("v (m/s)")
#     #plt.title(ti[frame]*2)
#     plt.ylim([-0.0008,0.0008])
    # plt.legend()

# def update(frame):
#     plt.clf()
#     plt.plot(X, Vfilm[frame,:])#, label=f'Frame {frame}')
#     plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -maxV,maxV], aspect='auto', alpha=0.25)
#     for i in range(alpha.size-1):
#         if modulateP.fonctionModK(pFilm*dt*frame,alpha[i],modulation,frontiere,i)<0.:
#             plt.axvline(x=alpha[i],color='red')#,linestyle='--')
#         else:
#             plt.axvline(x=alpha[i],color='lime')#,linestyle='--')
#     #plt.axvline(x=alpha[0],color='gray',linestyle='--')
#     # plt.colorbar()
#     # plt.plot(X, VfilmH[frame,:])
#     plt.xlabel("x (m)",fontsize=16)
#     plt.ylabel("v (m/s)",fontsize=16)
#     #plt.title(t[frame]*2)
#     plt.ylim([-maxV,maxV])
#     plt.tight_layout()
#     # plt.legend()


# total_frames = Vfilm.shape[0]
# fig, ax = plt.subplots()

# # Cree l'animation
# animation = FuncAnimation(fig, update, frames=total_frames, interval=2)
# # Affiche l'animation
# plt.show()
# # nomAnim="Film_SM_"+sourcef["Frequence"]+"Hz_"+modulation["Freq"]+"Hz.mp4"
# # animation.save(nomAnim, writer='ffmpeg', fps=100)

# plt.figure()
# plt.plot(X, Vfilm[0,:],linewidth=3)
# plt.axvline(x=alpha[0],color='black')
# plt.xlabel(r'$x$ (m)',fontsize=16)
# plt.ylabel(r'$v$ (m/s)',fontsize=16)
# plt.ylim([-maxV,maxV])
# plt.title("$t=$0 s",fontsize=16)
# plt.tight_layout()

# Ntfig=int(np.floor(0.045/dt))
# plt.figure()
# plt.plot(X, Vfilm[Ntfig,:],linewidth=3)
# plt.axvline(x=alpha[0],color='black')
# plt.xlabel(r'$x$ (m)',fontsize=16)
# plt.ylabel(r'$v$ (m/s)',fontsize=16)
# plt.ylim([-maxV,maxV])
# plt.title("$t=$"+str(round(t[Ntfig], 3))+" s",fontsize=16)
# plt.tight_layout()


# SolV,SolS=analytics.analytics_Cauchy_1IModK(0.045,X,sourcef,modulation,frontiere,mat,NOde,dt,Nt)
# # plt.figure()
# # plt.plot(trace)
# # plt.plot(trace2)

# plt.figure()
# plt.plot(X, Vfilm[Ntfig,:],linewidth=3,linestyle=None,marker='o')
# plt.plot(X,SolV)


# if float(frontiere["Kn_0"])==0:
#     K0=np.infty
# else:
#     K0=float(frontiere["Kn_0"])
# frK=float(modulation["Freq"])
# delta=float(modulation["DeltaF"])

# M0=float(frontiere["Mn_0"])
# frM=float(modulation["Freq"])
# deltaM=float(modulation["DeltaM"])

# modulb=modulateP.fonctionModM(t,alpha[0],modulation,frontiere,0)
# modulM=M0*(1+modulb)
# dmodulb=np.zeros(int(Nt/pFilm))
# # for i in range(int(Nt/pFilm)):
# #     modulb[i]=M0*(1+deltaM*np.sin(2*np.pi*frM*t[pFilm*i]))
# #     dmodulb[i]=2*np.pi*frM*M0*(deltaM*np.cos(2*np.pi*frM*t[pFilm*i]))
    
# modula=modulateP.fonctionModK(t,alpha[0],modulation,frontiere,0)
# modulK=1/K0*(1+modula)
# dmodula=np.zeros(int(Nt/pFilm))

# energyV=np.zeros(int(Nt/pFilm))
# energyI=np.zeros(int(Nt/pFilm))
# for ti in range(int(Nt/pFilm)):
#     energyV[ti]=dx/2*np.sum(rhox*Vfilm[ti]*Vfilm[ti])+dx/2*np.sum(1/rhox/Celx/Celx*Sfilm[ti]*Sfilm[ti])
#     energyI[ti]=1/2*modulM[ti]*CSVM[ti,2]*CSVM[ti,2]+1/2*modulK[ti]*CSVM[ti,3]*CSVM[ti,3]     
    
# energyM=energyV+energyI
# plt.figure()
# plt.plot(t[:-1],energyV,label=r'$\mathcal{E}_b$',linewidth=3)
# plt.plot(t[:-1],energyI,label=r'$\mathcal{E}_i$',linewidth=3)
# plt.plot(t[:-1],energyM,label=r'$\mathcal{E}_m$',linewidth=3)
# plt.legend(loc="lower right",fontsize=16)
# plt.xlabel(r'$t$ (s)',fontsize=16)
# plt.ylabel(r'$\mathcal{E}$ (J)',fontsize=16)
# plt.imshow([1+modula], cmap='gray', extent=[0, t[-1], 0,max(energyM)], aspect='auto', alpha=0.35)
# plt.colorbar().set_label(label=r'$1+\varepsilon~\sin (\Omega t)$',size=16)
# plt.xlim([0,0.06])
# plt.tight_layout()
# plt.figure()

# def errL2(A,B):
#     err=0
#     B2=0
#     L=max(B.shape)
#     for l in range(L):
#         eps=A[l]-B[l]
#         err=err+eps**2
#         B2=B2+B[l]**2
#     return err,err/B2

# ndx=7
# # err=np.zeros([4,ndx])
# errIm=np.zeros([2,ndx])
# # errR=np.zeros([4,ndx])
# errRIm=np.zeros([2,ndx])
# data=np.zeros([4,ndx])
# k0=2**((ndx-1)/2)
# # mat=inputReading.readfile('Milieu.txt')
# # frontiere=inputReading.readfile('Frontiere.txt')     
# # sourcef=inputReading.readfile('Source.txt')   
# # sourceRev=inputReading.readfile('SourceRev.txt')
# for k in range(0,ndx):
#     nDem='Demarrer'+str(k)+'.txt'
#     configuration=inputReading.readfile(nDem)
#     mat=inputReading.readfile('Milieu.txt')
#     frontiere=inputReading.readfile('Frontiere.txt')     
#     sourcef=inputReading.readfile('Source.txt')   
#     modulation=inputReading.readfile('Modulation.txt')
#     ################ Parameters of the simulation   
#     X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(configuration)
#     fs,CFL,tshift=inputReading.source(sourcef)
#     Nmat,rho,cm=inputReading.material(mat)
#     alpha=inputReading.frontiere(frontiere)
#     rhox,Celx=inputReading.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
#     dt=min(scheme1D.timeStep(CFL,dx,cm))
# #%
# #% Initial conditions
#     U0=scheme1D.initCauchyProblem(configuration,frontiere,mat,sourcef)
# ###### Mise en oeuvre
#     Ntfig=int(np.floor(0.045/dt))
#     Tmax=Ntfig*dt
#     tmax=str(Tmax)+"s"
#     t=np.linspace(0,Tmax,Ntfig+1)
#     sce=np.zeros(Ntfig+1)
#     for i in range(Ntfig+1):
#         sce[i]=source.choice_timefct(sourcef,t[i])
#     THV0y=np.zeros(Nx)
#     THVTMy=np.zeros(Nx)
#     THS0y=np.zeros(Nx)
#     THSTMy=np.zeros(Nx)
#     Nxalpha=int(alpha[0]/dx)
#     WL=max(cm/fs)
#     data[0,k]=Tmax
#     data[1,k]=dx
#     data[2,k]=dt
#     data[3,k]=WL/dx 
#     SolV,SolS=analytics.analytics_Cauchy_1IModK(Tmax,X,sourcef,modulation,frontiere,mat,NOde,dt,Ntfig+1)
#     nSve='ResV'+str(k)+'.txt'
#     nSveS='ResS'+str(k)+'.txt'
#     # Res=np.zeros([2,Nx])
#     # Res[0,:]=np.loadtxt(nSve, delimiter='\t')
#     # Res[1,:]=np.loadtxt(nSveS, delimiter='\t')
#     Res,Vfilm,Sfilm,CSVM=scheme1D.FD_cauchyII(U0,Ntfig+1,Nx,K,rhox,Celx,dt,dx,configuration,mat,frontiere,modulation,pFilm,rCSVM="yes")
#     errIm[0,k],errRIm[0,k]=np.sqrt(dx)*np.sqrt(errL2(Res[0,:], SolV))
#     errIm[1,k],errRIm[1,k]=np.sqrt(dx)*np.sqrt(errL2(Res[1,:], SolS))
#     np.savetxt(nSve, Res[0,:], fmt='%15e', delimiter='\t')
#     np.savetxt(nSveS, Res[1,:], fmt='%15e', delimiter='\t')
#     plt.figure()
#     plt.plot(X,SolV,'*')
#     plt.plot(X,Res[0,:])
#     plt.figure()
#     plt.plot(X,SolS,'*')
#     plt.plot(X,Res[1,:])
    
# plt.figure()
# plt.plot(np.log(data[1,:]),np.log(errIm[1,:]),'*')
# plt.plot(np.log(data[1,:]),np.log(data[1,:]))
# plt.plot(np.log(data[1,:]),2*np.log(data[1,:]))
# plt.plot(np.log(data[1,:]),4*np.log(data[1,:]))
# # plt.figure()
# # plt.plot(np.log(data[1,:]),np.log(errRIm[0,:]),'*')
# # plt.plot(np.log(data[1,:]),np.log(data[1,:]))
# # plt.plot(np.log(data[1,:]),2*np.log(data[1,:]))
# # plt.plot(np.log(data[1,:]),4*np.log(data[1,:]))

# pentConv=(np.log(errIm[1,6])-np.log(errIm[1,0]))/(np.log(data[1,6])-np.log(data[1,0]))

# plt.figure()
# plt.plot(X,Res[0,:])
# plt.plot(X,THVTMy,"--")
# plt.xlabel("$x$ (m)")
# plt.ylabel("$v$ (m/s)")
# plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0006,0.0006], aspect='auto', alpha=0.35)
# plt.colorbar()
# plt.axvline(x=alpha[0],color='gray',linestyle='--')


################ REVERSE PROPAGATION ####################################
# time0=time.time()
# ResR,VfilmR,SfilmR=scheme1D.FD_sourceVII(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0R,sceR,configuration,mat,frontiere,modulation,pFilm)
# timeR=time.time()-time0
# print(timeR)  
# def updateR(frame):
#     plt.clf()
#     plt.plot(X, VfilmR[frame,:])#, label=f'Frame {frame}')
#     plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0008,0.0008], aspect='auto', alpha=0.35)
#     plt.axvline(x=alpha[0],color='gray',linestyle='--')
#     plt.colorbar()
#     # plt.plot(X, VfilmH[frame,:])
#     plt.xlabel("x (m)")
#     plt.ylabel("v (m/s)")
#     plt.ylim([-0.0008,0.0008])
#     plt.title(t[frame]*2)
#     # plt.legend()

# total_frames = VfilmR.shape[0]
# figR, axR = plt.subplots()

# # Cree l'animation
# animationR = FuncAnimation(figR, updateR, frames=total_frames, interval=100)
# # Affiche l'animation
# plt.show()
# # nomAnimR="Film_SMR_"+sourcef["Frequence"]+"Hz.mp4"
# # animationR.save(nomAnimR, writer='ffmpeg', fps=20)


# ###################### NO MODULATION ##########################################################################
# time0=time.time()
# ResNM,VfilmNM,SfilmNM=scheme1D.FD_sourceVII(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0,sce,configuration,mat,frontiere,nomodulation,pFilm)
# timeNM=time.time()-time0
# print(timeNM)  

# # Fonction pour mettre à jour le contenu du graphique à chaque image
# def updateNM(frame):
#     plt.clf()
#     plt.plot(X, VfilmNM[frame,:])#, label=f'Frame {frame}')
#     plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0008,0.0008], aspect='auto', alpha=0.35)
#     plt.axvline(x=alpha[0],color='gray',linestyle='--')
#     #plt.colorbar()
#     # plt.plot(X, VfilmH[frame,:])
#     plt.xlabel("x (m)")
#     plt.ylabel("v (m/s)")
#     plt.ylim([-maxV,maxV])
#     # plt.ylim([-0.0008,0.0008])
#     #plt.title(t[frame]*2)
#     # plt.legend()


# total_frames = VfilmNM.shape[0]
# figNM, axNM = plt.subplots()

# # Cree l'animation
# animationNM = FuncAnimation(figNM, updateNM, frames=total_frames, interval=100)
# # Affiche l'animation
# plt.show()
# nomAnimNM="Film_SMNM_"+sourcef["Frequence"]+"Hz.mp4"
# animationNM.save(nomAnimNM, writer='ffmpeg', fps=20)


# #######################################
# plt.figure()
# plt.clf()
# plt.plot(Vfilm[:,int(x_0R)])
# plt.plot(VfilmNM[:,int(x_0R)])
# plt.xlabel("t (s)")

# plt.figure()
# plt.clf()
# plt.plot(VfilmR[:,int(x_0)])
# plt.plot(Vfilm[:,int(x_0R)],"--")
# plt.xlabel("t (s)")

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
# plt.plot(X,b0y.T,'--')
# plt.xlabel("$x$ (m)")
# plt.ylabel("$\sigma$ (Pa)")\includegraphics[scale=0.52]{1Interface_t0} & 

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
