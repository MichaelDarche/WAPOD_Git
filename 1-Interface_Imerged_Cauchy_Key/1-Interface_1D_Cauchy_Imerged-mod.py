#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:50:43 2023

@author: michael
"""

import os
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/1-Interface_Imerged_Cauchy_Key")
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
################# Input files 
#%
configuration=inputReading.readfile('Demarrer.txt')
mat=inputReading.readfile('Milieu.txt')
frontiere=inputReading.readfile('Frontiere.txt')     
sourcef=inputReading.readfile('Source.txt')   
# sourceRev=inputReading.readfile('SourceRev.txt')
modulation=inputReading.readfile('Modulation.txt')
modulationR=inputReading.readfile('ModulationR.txt')
# nomodulation=inputReading.readfile('NoModulation.txt')
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




#% Initial conditions
U0=scheme1D.initCauchyProblem(configuration,frontiere,mat,sourcef)
# U0=scheme1D.initPointSourceProblem(configuration,frontiere,mat,sourcef)

###### Mise en oeuvre
Nt=1000
Tmax=Nt*dt
tmax=str(Tmax)+"s"
K=4

t=np.linspace(0,Tmax,Nt+1)

sce=np.zeros(Nt+1)
for i in range(Nt+1):
    sce[i]=source.choice_timefct(sourcef,t[i])
    
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
fftSig=1/Nt*np.fft.fftshift(np.fft.fft(sigcor,n=N))
freq = np.fft.fftshift(np.fft.fftfreq(N,d=dt))
plt.figure()
plt.plot(np.abs(freq),np.abs(fftSig),linewidth=3)
plt.xlim([0,300])
plt.ylim([0,0.25])
plt.xlabel(r'$f$ (Hz)',fontsize=16)
# plt.ylabel(r'$\overline{S}$',fontsize=16)

plt.ylabel(r'$\hat{S}$',fontsize=16)
plt.axvline(x=fs,color='black',linestyle='--')
plt.tight_layout()




    
    
# sceR=np.zeros(Nt+1)
# for i in range(Nt+1):
#     sceR[i]=source.choice_timefct(sourceRev,t[i])
    
    
from matplotlib.animation import FuncAnimation

pFilm=1
############### MODULATION DIRECT PROPAGATION ##########################################
import time
time0=time.time()
Res,Vfilm,Sfilm,CSVM=scheme1D.FD_cauchyII(U0,Nt,Nx,K,rhox,Celx,dt,dx,configuration,mat,frontiere,modulation,pFilm,rCSVM="yes")
#U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
timeM=time.time()-time0
maxV=np.max([np.max(Vfilm),-np.min(Vfilm)])
print(timeM)  

import time
time0=time.time()
ResR,VfilmR,SfilmR=scheme1D.FD_cauchyIIKey(Res,Nt,Nx,K,rhox,Celx,dt,dx,configuration,mat,frontiere,modulationR,pFilm,rCSVM="no")
#U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
timeM=time.time()-time0
maxV=np.max([np.max(Vfilm),-np.min(Vfilm)])
print(timeM)  

VfilmTot=np.zeros([2*Nt,Nx])
VfilmTot[0:Nt,:]=Vfilm
VfilmTot[Nt:2*Nt,:]=VfilmR
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

def update(frame):
    plt.clf()
    plt.plot(X, VfilmTot[frame,:])#, label=f'Frame {frame}')
    plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -maxV,maxV], aspect='auto', alpha=0.25)
    for i in range(alpha.size-1):
        if modulateP.fonctionModK(pFilm*dt*frame,alpha[i],modulation,frontiere,i)<0.:
            plt.axvline(x=alpha[i],color='red')#,linestyle='--')
        else:
            plt.axvline(x=alpha[i],color='lime')#,linestyle='--')
    #plt.axvline(x=alpha[0],color='gray',linestyle='--')
    # plt.colorbar()
    # plt.plot(X, VfilmH[frame,:])
    plt.xlabel("x (m)",fontsize=16)
    plt.ylabel("v (m/s)",fontsize=16)
    #plt.title(t[frame]*2)
    plt.ylim([-maxV,maxV])
    plt.tight_layout()
    # plt.legend()


total_frames = 2*Vfilm.shape[0]
fig, ax = plt.subplots()

# Cree l'animation
animation = FuncAnimation(fig, update, frames=total_frames, interval=2)
# Affiche l'animation
plt.show()
# nomAnim="Film_Key3_"+sourcef["Frequence"]+"Hz_"+modulation["Freq"]+"Hz.mp4"
# animation.save(nomAnim, writer='ffmpeg', fps=100)

plt.figure()
plt.plot(X, U0[0,:],linewidth=3,label="source")
plt.plot(X, VfilmTot[-1,:],"--",linewidth=2,label="time-reversal")
plt.plot(X, Res[0,:],linewidth=3,label="crypto")
plt.axvline(x=alpha[0],color='black')
plt.xlabel(r'$x$ (m)',fontsize=16)
plt.ylabel(r'$v$ (m/s)',fontsize=16)
plt.ylim([-maxV,maxV])
plt.legend()
# plt.title("$t=$0 s",fontsize=16)
plt.tight_layout()

Ntfig=int(np.floor(0.045/dt))
plt.figure()
plt.plot(X, Vfilm[Ntfig,:],linewidth=3)
plt.axvline(x=alpha[0],color='black')
plt.xlabel(r'$x$ (m)',fontsize=16)
plt.ylabel(r'$v$ (m/s)',fontsize=16)
plt.ylim([-maxV,maxV])
plt.title("$t=$"+str(round(t[Ntfig], 3))+" s",fontsize=16)
plt.tight_layout()

import analytics
# NOde=5000
# SolV,SolS=analytics.analytics_Cauchy_1IModK(0.045,X,sourcef,modulation,frontiere,mat,NOde,dt,Nt)
# plt.figure()
# plt.plot(trace)
# plt.plot(trace2)

plt.figure()
plt.axvline(x=alpha[0],color='black')
plt.plot(X, Vfilm[Ntfig,:],linestyle=None,marker='o', label='numerics')
# plt.plot(X,SolV,linewidth=2, label='analytics')
plt.xlabel(r'$x$ (m)',fontsize=16)
plt.ylabel(r'$v$ (m/s)',fontsize=16)
plt.legend(loc="lower left",fontsize=16)
plt.tight_layout()

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
dmodula=np.zeros(int(Nt/pFilm))

energyV=np.zeros(int(Nt/pFilm))
energyI=np.zeros(int(Nt/pFilm))
for ti in range(int(Nt/pFilm)):
    energyV[ti]=dx/2*np.sum(rhox*Vfilm[ti]*Vfilm[ti])+dx/2*np.sum(1/rhox/Celx/Celx*Sfilm[ti]*Sfilm[ti])
    energyI[ti]=1/2*modulM[ti]*CSVM[ti,2]*CSVM[ti,2]+1/2*modulK[ti]*CSVM[ti,3]*CSVM[ti,3]     
    
energyM=energyV+energyI
plt.figure()
plt.plot(t[:-1],energyV,label=r'$\mathcal{E}_b$',linewidth=3)
plt.plot(t[:-1],energyI,label=r'$\mathcal{E}_i$',linewidth=3)
plt.plot(t[:-1],energyM,label=r'$\mathcal{E}_m$',linewidth=3)
plt.legend(loc="lower right",fontsize=16)
plt.xlabel(r'$t$ (s)',fontsize=16)
plt.ylabel(r'$\mathcal{E}$ (J)',fontsize=16)
plt.imshow([1+modula], cmap='gray', extent=[0, t[-1], 0,max(energyM)], aspect='auto', alpha=0.35)
plt.colorbar().set_label(label=r'$1+\varepsilon~\sin (\Omega t)$',size=16)
plt.xlim([0,0.06])
plt.tight_layout()
plt.figure()



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
# plt.plot(X,LWS0y.T,'--')
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
