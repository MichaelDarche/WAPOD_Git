#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:50:43 2023

@author: michael
"""

import os
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/1-Interface_Imerged_Cauchy_modulated_errM")
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

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
################# Input files 
#%
# configuration=inputReading.readfile('Demarrer.txt')
# mat=inputReading.readfile('Milieu.txt')
# frontiere=inputReading.readfile('Frontiere.txt')     
# sourcef=inputReading.readfile('Source.txt')   
# # sourceRev=inputReading.readfile('SourceRev.txt')
# modulation=inputReading.readfile('Modulation.txt')
# # nomodulation=inputReading.readfile('NoModulation.txt')
# ################ Parameters of the simulation   
# X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(configuration)
# fs,CFL,x_0=inputReading.source(sourcef)
# # fsR,CFLR,x_0R=inputReading.source(sourceRev)
# Nmat,rho,cm=inputReading.material(mat)
# alpha=inputReading.frontiere(frontiere)
# rhox,Celx=inputReading.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
# dt=min(scheme1D.timeStep(CFL,dx,cm))
# # lambdaOnde=cm/fs
# xsmin=cm*tshift-lambdaOnde
# xsmax=cm*tshift
#%




# #% Initial conditions
# U0=scheme1D.initCauchyProblem(configuration,frontiere,mat,sourcef)

# ###### Mise en oeuvre
# Nt=4000
# Tmax=Nt*dt
# tmax=str(Tmax)+"s"
# K=4

# t=np.linspace(0,Tmax,Nt+1)

# sce=np.zeros(Nt+1)
# for i in range(Nt+1):
#     sce[i]=source.choice_timefct(sourcef,t[i])
    
# def cosine_taper(signal, alpha=0.00,zp=18):
#     zeropadding=2**zp
#     Size=signal.shape[0]
#     M=np.floor((Size*alpha)/2+0.5);
#     tapered=np.zeros(zeropadding)
#     for j in range(int(Size)):
#         if j<=M+1:
#             tapered[j]=signal[j] * (0.5 * ( 1-np.cos(j*np.pi/(M+1))))
#         elif j<Size - M-1:
#             tapered[j]=signal[j]
#         elif j<=Size:
#             tapered[j]=signal[j]* (0.5 * (1-np.cos((Size-j)*np.pi/(M+1))))
#     return tapered

# N=2**19
# sigcor = cosine_taper(sce)
# fftSig=1/Nt*np.fft.fftshift(np.fft.fft(sigcor,n=N))
# freq = np.fft.fftshift(np.fft.fftfreq(N,d=dt))
# plt.figure()
# plt.plot(np.abs(freq),np.abs(fftSig),linewidth=3)
# plt.xlim([0,300])
# plt.ylim([0,0.125])
# plt.xlabel(r'$f$ (Hz)',fontsize=16)
# # plt.ylabel(r'$\overline{S}$',fontsize=16)

# plt.ylabel(r'$\hat{S}$',fontsize=16)
# plt.axvline(x=fs,color='black',linestyle='--')
# plt.tight_layout()




    
    
# sceR=np.zeros(Nt+1)
# for i in range(Nt+1):
#     sceR[i]=source.choice_timefct(sourceRev,t[i])
    
    
from matplotlib.animation import FuncAnimation

pFilm=1
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

import analytics
NOde=20001
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

def errL2(A,B):
    err=0
    B2=0
    L=max(B.shape)
    for l in range(L):
        eps=A[l]-B[l]
        err=err+eps**2
        B2=B2+B[l]**2
    return err,err/B2

ndx=7
# err=np.zeros([4,ndx])
errIm=np.zeros([6,ndx])
# errR=np.zeros([4,ndx])
errRIm=np.zeros([6,ndx])
data=np.zeros([4,ndx])
# k0=2**((ndx-1)/2)
# mat=inputReading.readfile('Milieu.txt')
# frontiere=inputReading.readfile('Frontiere.txt')     
# sourcef=inputReading.readfile('Source.txt')   
# sourceRev=inputReading.readfile('SourceRev.txt')
for k in range(0,ndx):
    K=4
    nDem='Demarrer'+str(k)+'.txt'
    configuration=inputReading.readfile(nDem)
    mat=inputReading.readfile('Milieu.txt')
    frontiere=inputReading.readfile('Frontiere.txt')     
    sourcef=inputReading.readfile('Source.txt')   
    modulation0=inputReading.readfile('Modulation0.txt')
    modulation10=inputReading.readfile('Modulation10.txt')
    modulation100=inputReading.readfile('Modulation100.txt')
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
    Ntfig=int(np.floor(0.045/dt))
    Tmax=Ntfig*dt
    tmax=str(Tmax)+"s"
    t=np.linspace(0,Tmax,Ntfig+1)
    sce=np.zeros(Ntfig+1)
    for i in range(Ntfig+1):
        sce[i]=source.choice_timefct(sourcef,t[i])
    THVTMy=np.zeros(Nx)
    THSTMy=np.zeros(Nx)
    THVTMyF=np.zeros(Nx)
    THSTMyF=np.zeros(Nx)
    Nxalpha=int(alpha[0]/dx)
    WL=max(cm/fs)
    data[0,k]=Tmax
    data[1,k]=dx
    data[2,k]=dt
    data[3,k]=WL/dx 
    # NFourier=128
    # Kf=10
    # df=Kf*fs/NFourier
    # for pt in range(Nx):
    #     THVTMyF[pt],THSTMyF[pt]=analytics.analytics1InterfaceCauchy(tshift+Tmax,X[pt],sourcef,frontiere,NFourier,mat,Kf)
    # SolV,SolS=analytics.analytics_Cauchy_1IModM(Tmax,X,sourcef,modulation,frontiere,mat,NOde,dt,Ntfig+1)
    # R01=(rho[1]*cm[1]-rho[0]*cm[0])/(rho[1]*cm[1]+rho[0]*cm[0])
    # T01=2*rho[0]*cm[1]/(rho[1]*cm[1]+rho[0]*cm[0])
    # Nxalpha=int(alpha[0]/dx)
    # for x in range(Nxalpha):
    #     f0=source.choice_timefct(sourcef,tshift-x*dx/cm[0])
    #     fr=source.choice_timefct(sourcef,tshift+x*dx/cm[0]-2*alpha[0]/cm[0])
    #     fT=source.choice_timefct(sourcef,Tmax+tshift-x*dx/cm[0])
    #     frT=source.choice_timefct(sourcef,Tmax+tshift+x*dx/cm[0]-2*alpha[0]/cm[0])
    #     THVTMy[x]=1/cm[0]*fT-R01/cm[0]*frT
    #     THSTMy[x]=-rho[0]*fT-R01*rho[0]*frT
    # for x in range(Nxalpha,Nx):
    #     ft=source.choice_timefct(sourcef,tshift-x*dx/cm[1]+(1/cm[1]-1/cm[0])*alpha[0])
    #     ftT=source.choice_timefct(sourcef,Tmax+tshift-x*dx/cm[1]+(1/cm[1]-1/cm[0])*alpha[0])
    #     THVTMy[x]=1/cm[1]*T01*ftT
    #     THSTMy[x]=-rho[1]*T01*ftT
    nSve='ResV'+str(k)+str(100)+'.txt'
    nSveS='ResS'+str(k)+str(100)+'.txt'
    Res=np.zeros([2,Nx])
    Res[0,:]=np.loadtxt(nSve, delimiter='\t')
    Res[1,:]=np.loadtxt(nSveS, delimiter='\t')
    SolV,SolS=analytics.analytics_Cauchy_1IModM(Tmax,X,sourcef,modulation100,frontiere,mat,NOde,dt,Ntfig+1)
    Res,Vfilm,Sfilm,CSVM=scheme1D.FD_cauchyII(U0,Ntfig,Nx,K,rhox,Celx,dt,dx,configuration,mat,frontiere,modulation100,pFilm,rCSVM="yes")
    # np.savetxt(nSve, Res[0,:], fmt='%18e', delimiter='\t')
    # np.savetxt(nSveS, Res[1,:], fmt='%18e', delimiter='\t')
    errIm[0,k],errRIm[0,k]=np.sqrt(dx)*np.sqrt(errL2(Res[0,:], SolV))
    errIm[1,k],errRIm[1,k]=np.sqrt(dx)*np.sqrt(errL2(Res[1,:], SolS))
    ################################################
    nSve='ResV'+str(k)+str(10)+'.txt'
    nSveS='ResS'+str(k)+str(10)+'.txt'
    Res=np.zeros([2,Nx])
    Res[0,:]=np.loadtxt(nSve, delimiter='\t')
    Res[1,:]=np.loadtxt(nSveS, delimiter='\t')
    SolV,SolS=analytics.analytics_Cauchy_1IModM(Tmax,X,sourcef,modulation10,frontiere,mat,NOde,dt,Ntfig+1)
    Res,Vfilm,Sfilm,CSVM=scheme1D.FD_cauchyII(U0,Ntfig,Nx,K,rhox,Celx,dt,dx,configuration,mat,frontiere,modulation10,pFilm,rCSVM="yes")
    errIm[2,k],errRIm[2,k]=np.sqrt(dx)*np.sqrt(errL2(Res[0,:], SolV))
    errIm[3,k],errRIm[3,k]=np.sqrt(dx)*np.sqrt(errL2(Res[1,:], SolS))
    # np.savetxt(nSve, Res[0,:], fmt='%18e', delimiter='\t')
    # np.savetxt(nSveS, Res[1,:], fmt='%18e', delimiter='\t')
    ################################################
    nSve='ResV'+str(k)+str(0)+'.txt'
    nSveS='ResS'+str(k)+str(0)+'.txt'
    Res=np.zeros([2,Nx])
    Res[0,:]=np.loadtxt(nSve, delimiter='\t')
    Res[1,:]=np.loadtxt(nSveS, delimiter='\t')
    SolV,SolS=analytics.analytics_Cauchy_1IModM(Tmax,X,sourcef,modulation0,frontiere,mat,NOde,dt,Ntfig+1)
    # Res,Vfilm,Sfilm,CSVM=scheme1D.FD_cauchyII(U0,Ntfig,Nx,K,rhox,Celx,dt,dx,configuration,mat,frontiere,modulation0,pFilm,rCSVM="yes")
    errIm[4,k],errRIm[2,k]=np.sqrt(dx)*np.sqrt(errL2(Res[0,:], SolV))
    errIm[5,k],errRIm[3,k]=np.sqrt(dx)*np.sqrt(errL2(Res[1,:], SolS))
    # np.savetxt(nSve, Res[0,:], fmt='%18e', delimiter='\t')
    # np.savetxt(nSveS, Res[1,:], fmt='%18e', delimiter='\t')
    # plt.figure()
    # plt.plot(X,SolV,'*')
    # plt.plot(X,Res[0,:])
    # plt.figure()
    # plt.plot(X,SolS,'*')
    # plt.plot(X,Res[1,:])
    
plt.figure()
plt.loglog(data[1,1:-1],errIm[0,1:-1],'*-',label=r"$f_m=100$ Hz")
plt.loglog(data[1,1:-1],errIm[2,1:-1],'o--',label=r"$f_m=10$ Hz")
# plt.loglog(data[1,1:-1],errIm[4,1:-1],'*-')
# plt.plot(np.log(data[1,:]),np.log(errIm[1,:]),'*')
# plt.plot(np.log(data[1,:]),np.log(errIm[3,:]),'*')
# plt.plot(np.log(data[1,:]),np.log(errIm[5,:]),'*')
# plt.loglog(data[1,:],(data[1,:]))
# plt.loglog(data[1,:],(data[1,:])**2)
plt.loglog(data[1,1:-1],(data[1,1:-1])**4/200000,":",label=r"$\alpha\Delta_x^4$")
plt.xlabel(r'$\Delta_x$ (m)',fontsize=16)
plt.ylabel(r"$\varepsilon_v$ (m/s)",fontsize=16)
plt.legend(fontsize=16)
# plt.title("$\mathscr{M}(t)$,$\mathscr{C}(t)=0$")

# plt.figure()
# plt.plot(np.log(data[1,:]),np.log(errRIm[0,:]),'*')
# plt.plot(np.log(data[1,:]),np.log(data[1,:]))
# plt.plot(np.log(data[1,:]),2*np.log(data[1,:]))
# plt.plot(np.log(data[1,:]),4*np.log(data[1,:]))

pentConv=(np.log(errIm[1,6])-np.log(errIm[1,0]))/(np.log(data[1,6])-np.log(data[1,0]))

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
