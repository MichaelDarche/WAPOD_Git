#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:50:43 2023

@author: michael
"""

import os
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/TimeLaminated_Cauchy/")

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../WAPOD/")

import scheme1D_TDold as scheme1D
import source
import inputReading




################# Input files 
#%
configuration=inputReading.readfile('Demarrer.txt')
mat=inputReading.readfile('Milieu.txt')
frontiere=inputReading.readfile('Frontiere.txt')     
sourcef=inputReading.readfile('Source.txt')   
# sourceRev=inputReading.readfile('SourceRev.txt')
modulation=inputReading.readfile('Modulation.txt')
# nomodulation=inputReading.readfile('NoModulation.txt')
################ Parameters of the simulation   
X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(configuration)
fs,CFL,x_0=inputReading.source(sourcef)
# fs,CFL,x_0R=inputReading.source(sourceRev)
Nmat,rho,cm=inputReading.material(mat)
alpha=inputReading.frontiere(frontiere)
rhox,Celx=inputReading.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
dt=min(scheme1D.timeStep(CFL,dx,cm))
# lambdaOnde=cm/fs
# xsmin=cm*tshift-lambdaOnde
# xsmax=cm*tshift
#%

#% Initial conditions
U0=scheme1D.initCauchyProblemH(configuration,frontiere,mat,sourcef,modulation)

###### Mise en oeuvre
Nt=2**10
Tmax=Nt*dt
tmax=str(Tmax)+"s"
K=4

t=np.linspace(0,Tmax,Nt+1)


import time,homogenize,math
time0=time.time()
fCFLh=homogenize.updateCFL(Nt,dt,t,rho,cm,frontiere,modulation)
CFLH=CFL/fCFLh
dt2v=np.min(scheme1D.timeStep(CFLH,dx,cm))
kt=math.ceil(dt/dt2v)
dt2p=min(dt2v,dt)
Nt2=int(Tmax/dt2p)
Ntp=max(Nt2,Nt)
tmax=str(Tmax)+"s"

tp=np.linspace(0,Tmax,Ntp+1)


sce=np.zeros(Ntp+1)
for i in range(Ntp+1):
    sce[i]=source.choice_timefct(sourcef,t[i])
    
# sceR=np.zeros(Nt+1)
# for i in range(Nt+1):
#     sceR[i]=source.choice_timefct(sourceRev,t[i])
    
    
from matplotlib.animation import FuncAnimation

pFilm=1
############### MODULATION DIRECT PROPAGATION ##########################################
import time
time0=time.time()
Res,Vfilm,Sfilm=scheme1D.FD_CauchyVIIHomoGS_Period(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,configuration,mat,frontiere,modulation,pFilm)
#U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
timeM=time.time()-time0
maxV=np.max([np.max(Vfilm),-np.min(Vfilm)])
print(timeM)  

import homogenize
rhoT,CT=homogenize.timeProp(Nt, dt, t, rho, cm, frontiere, modulation)


# Fonction pour mettre à jour le contenu du graphique à chaque image
def update(frame):
    plt.clf()
    plt.plot(X, Vfilm[frame,:])#, label=f'Frame {frame}')
    plt.imshow([np.linspace(rhoT[frame],rhoT[frame],Nx)], cmap='seismic',vmax=max(rhoT),vmin=min(rhoT), extent=[xmin, xmax,-0,maxV], aspect='auto', alpha=0.35)
    plt.colorbar() 
    plt.imshow([np.linspace(CT[frame],CT[frame],Nx)], cmap='gray',vmax=max(CT),vmin=min(CT), extent=[xmin, xmax,-maxV,0], aspect='auto', alpha=0.35)
    plt.colorbar() 
    # plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -maxV,maxV], aspect='auto', alpha=0.35)
    # for i in range(alpha.size):
    #     plt.axvline(x=alpha[i],color='gray',linestyle='--')
    
    # plt.plot(X, VfilmH[frame,:])
    plt.xlabel("x (m)")
    plt.ylabel("v (m/s)")
    # plt.title(t[frame]*2)
    plt.ylim([-maxV,maxV])
    # plt.legend()

energyV=np.zeros(int(Nt/pFilm))
for ti in range(int(Nt/pFilm)):
   energyV[ti]=dx/2*np.sum(rhoT[ti*pFilm]*Vfilm[ti]*Vfilm[ti])+1/2*np.sum(1/rhoT[ti*pFilm]/CT[ti*pFilm]/CT[ti*pFilm]*Sfilm[ti]*Sfilm[ti])


total_frames = Vfilm.shape[0]
fig, ax = plt.subplots()

# Cree l'animation
animation = FuncAnimation(fig, update, frames=total_frames, interval=1)
# Affiche l'animation
plt.show()
# nomAnim="Film_SM_"+sourcef["Frequence"]+"Hz_"+modulation["FreqM"]+"Hz.mp4"
# animation.save(nomAnim, writer='ffmpeg', fps=100)
plt.figure()
plt.plot(tp[:-1],energyV)
plt.imshow([rhoT[::pFilm]], cmap='gray', extent=[t[0], t[-2], 0,np.max(energyV)], aspect='auto', alpha=0.35)
# plt.title("$$ (m/s)")
#############DISPERSION###########
h=alpha[1]-alpha[0]
K0=float(frontiere["Kn_0"])
DK=float(modulation["DeltaF"])
K1=K0/(1+DK)
K2=K0/(1-DK)
M0=float(frontiere["Mn_0"])
DM=float(modulation["DeltaM"])
M1=M0*(1+DM)
M2=M0*(1-DM)

cmet=cm[0]
rhoet=rho[0]
Eet=rhoet*cmet**2

K0et=K0*h/Eet
M0et=M0/h/rhoet
K1et=K1*h/Eet
M1et=M1/h/rhoet
K2et=K2*h/Eet
M2et=M2/h/rhoet

rhoef1=rhoet+M1/h
Eef1=K1*h/(K1*h/Eet+1)
cmef1=np.sqrt(Eef1/rhoef1)
rhoef2=rhoet+M2/h
Eef2=K2*h/(K2*h/Eet+1)
cmef2=np.sqrt(Eef2/rhoef2)

rhoeff0=M0/h+rhoet
deltarho=M0/h/rhoeff0*DM

Eeff0=1/(1/K0/h+1/Eet)
deltaE=Eeff0/h/K0*DK

RhoZero=1/2*(rhoef1+rhoef2)
DeltaRho=(rhoef1-rhoef2)/RhoZero/2
EZero=2/(1/Eef1+1/Eef2)
DeltaE=2*EZero*(1/Eef1-1/Eef2)
CZero=np.sqrt(EZero/RhoZero)
LambdaZero=CZero/fs
kZero=2*np.pi/LambdaZero

ratio=float(modulation["Ratio"])
c_r = cmef1*ratio+cmef2*(1-ratio)


fm=float(modulation["Freq"])
tau=1/fm

Z1 = np.sqrt(rhoef1*Eef1) 
Z2 = np.sqrt(rhoef2*Eef2)


L1 = ratio*tau 
L2 = (1-ratio)*tau 


def D_micro(k,omega):
    Dmicro=np.cos(omega*tau)-np.cos(cmef1*k*L1)*np.cos(cmef2*k*L2)+1/2*(Z1/Z2+Z2/Z1)*np.sin(cmef1*k*L1)*np.sin(cmef2*k*L2)
    return Dmicro

Nome=150
Nk=150

Omega_plot = np.linspace(0,2*np.pi,Nome)/tau 
k_map = np.linspace(0,4*np.pi,Nk)/tau/c_r;


DD_micro = np.zeros([Nome,Nk])
DD_det = np.zeros([Nome,Nk])


for ind in range(Nk):
    for ind2 in range(Nome):
        DD_micro[ind2,ind] = D_micro(k_map[ind],Omega_plot[ind2])
        #DD_det(ind2,ind) = Det_micro(k_map(ind),Omega_plot(ind2)



plt.figure()
plt.pcolor(k_map*c_r*tau,Omega_plot*tau,np.log10(np.abs(DD_micro)))
plt.axvline(x=kZero*c_r*tau,color='gray',linestyle='--')
plt.figure()
########################################



N=400
def cosine_taper(signal, alpha=0.12):
    Size=signal.shape[0]
    M=np.floor((Size*alpha)/2+0.5);
    tapered=np.zeros(Size)
    for j in range(int(Size)):
        if j<=M+1:
            tapered[j]=signal[j] * (0.5 * ( 1-np.cos(j*np.pi/(M+1))))
        elif j<Size - M-1:
            tapered[j]=signal[j]
        elif j<=Size:
            tapered[j]=signal[j] * (0.5 * (1-np.cos((Size-j)*np.pi/(M+1))))
    return tapered

# sigcor = cosine_taper(Vfilm[::,int(480)])
# plt.plot(Vfilm[::,int(480)])
# plt.plot(sigcor)
# fftSig=np.fft.fftshift(np.fft.fft(sigcor,n=N))
# freq = np.fft.fftshift(np.fft.fftfreq(N,d=dt))
# plt.figure()
# plt.plot(np.abs(freq),np.abs(fftSig))

# sigcor = cosine_taper(Vfilm[::,int(380)])
# plt.plot(Vfilm[::,int(380)])
# plt.plot(sigcor)
# fftSig=np.fft.fftshift(np.fft.fft(sigcor,n=N))
# freq = np.fft.fftshift(np.fft.fftfreq(N,d=dt))
# plt.figure()
# plt.plot(np.abs(freq),np.abs(fftSig))

# sigcor = cosine_taper(Vfilm[::,int(350)])
# plt.plot(Vfilm[::,int(350)])
# plt.plot(sigcor)
# fftSig=np.fft.fftshift(np.fft.fft(sigcor,n=N))
# freq = np.fft.fftshift(np.fft.fftfreq(N,d=dt))
# plt.figure()
# plt.plot(np.abs(freq),np.abs(fftSig))

def specfk(U,dx,dt,Nt,Nx):
    """
    FK spectrum of traces using the numpy.fft functions.
    """
    # Time padding
    zp=2**15
    Uzp=np.zeros([zp,Nx])
    Uzp[:Nt-1,:]=U
    #Fourier in time domain 
    FFTt = np.fft.rfft(Uzp, axis=0)
    #Fourier in sapce
    FFTx = np.fft.fft(FFTt, axis=1)
    FFT = np.flip(np.fft.fftshift(FFTx, axes=1), axis=1)
    FFK = np.absolute(FFT)

    # Get the frequency and K vectors
    frqv = np.fft.rfftfreq(zp, dt)
    wavv = np.fft.fftfreq(Nx, dx)
    # Centering
    wavv = np.fft.fftshift(wavv)
    return frqv,wavv,FFK

def specfkR(U,dx,dt,Nt,Nx):
    """
    FK spectrum of traces using the numpy.fft functions.
    """
    # Time padding
    zt=2**13
    zx=2**13
    Uzp=np.zeros([zt,zx])
    Uzp[:Nt-1,int(zx/2-Nx/2):int(zx/2+Nx/2)]=U
    #Fourier in time domain 
    FFTx = np.fft.fft(Uzp, axis=1)
    #Fourier in space
    FFTt = np.fft.rfft(FFTx, axis=0)
    FFT = np.flip(np.fft.fftshift(FFTt, axes=1), axis=1)
    FFK = np.absolute(FFT)

    # Get the frequency and K vectors
    frqv = np.fft.rfftfreq(zt, dt)
    wavv = np.fft.fftfreq(zx, dx)
    # Centering
    wavv = np.fft.fftshift(wavv)
    return frqv,wavv,FFK

# # frqv,wavv,FFK=specfk(Vfilm, dx, dt2p*pFilm, int(Ntp/pFilm)+1, Nx)
# # frqvS,wavvS,FFKS=specfkR(VfilmSS, dx, dt2p*pFilm, int(Ntp/pFilm)+1, Nx)
# frqvG,wavvG,FFKG=specfkR(Vfilm[::,::], dx, dt2p*pFilm, int(Ntp/pFilm)+1, int((Nx)))
# # plt.figure()
# # # plt.pcolor(wavv[::]*cmet*tau,frqv[::2]*tau,np.log(FFK[::2,::]),vmin=-12, vmax=2)
# # # plt.ylim([0,2])
# # # plt.xlim([0,12])
# # plt.colorbar()
# # plt.xlabel("k (1/m)")
# # plt.ylabel("f (Hz)")
# # plt.tight_layout()
# plt.figure()
# plt.pcolor(wavvG[::]*tau*c_r,frqvG[::]*tau,np.log(FFKG[::,::]),vmin=-6, vmax=2)
# # plt.ylim([0,8*np.pi])
# # plt.xlim([-0,8*np.pi])
# plt.colorbar()
# plt.xlabel("kc^*\tau")
# plt.ylabel("f*\tau")
# plt.tight_layout()
# nomFK="FK_JK_KM_"+sourcef["Frequence"]+"_"+modulation["Freq"]+"Hz.eps"
# plt.savefig(nomFK, format='eps')
# nomFK="FK_JK_KM_"+sourcef["Frequence"]+"_"+modulation["Freq"]+"Hz.png"
# plt.savefig(nomFK, format='png')
# left=FFK[:,0:500]
# right=FFK[:,501:1001]
# right=right[:,::-1]
# diff=np.abs(left-right)
# plt.figure()
# plt.pcolor(wavv[0:500],frqv,(diff))
# plt.colorbar()
# plt.ylim([0,200])
# plt.xlim([-np.pi/10,0])

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
#     plt.colorbar()
#     # plt.plot(X, VfilmH[frame,:])
#     plt.xlabel("x (m)")
#     plt.ylabel("v (m/s)")
#     plt.ylim([-0.0008,0.0008])
#     plt.title(t[frame]*2)
#     # plt.legend()


# total_frames = VfilmNM.shape[0]
# figNM, axNM = plt.subplots()

# # Cree l'animation
# animationNM = FuncAnimation(figNM, updateNM, frames=total_frames, interval=100)
# # Affiche l'animation
# plt.show()
# # nomAnimNM="Film_SMNM_"+sourcef["Frequence"]+"Hz.mp4"
# # animationNM.save(nomAnimNM, writer='ffmpeg', fps=20)


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
