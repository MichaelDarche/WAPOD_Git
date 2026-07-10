#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:50:43 2023

@author: michael
"""

import os
#os.chdir("/home/darche/Documents/ModulatedMetamaterials/WAPRO_TD/N-Interface_Imerged_PS_modulated/")
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/N-Interface_Imerged_PS_modulated_LFB2/")

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../WAPOD/")

import scheme1D_TDold as scheme1D
import source
import inputReading
import modulation as modul



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
Nt=2**13
Tmax=Nt*dt
tmax=str(Tmax)+"s"

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
K=4

tp=np.linspace(0,Tmax,Ntp+1)

sce=np.zeros(Ntp+1)
for i in range(Ntp+1):
    sce[i]=source.choice_timefct(sourcef,tp[i])
    
sceR=np.zeros(Ntp+1)
for i in range(Ntp+1):
    sceR[i]=source.choice_timefct(sourceRev,tp[i])
    
    
from matplotlib.animation import FuncAnimation

pFilm=4
############### MODULATION DIRECT PROPAGATION ##########################################
import time
time0=time.time()
# Res,Vfilm,Sfilm=scheme1D.FD_sourceVII(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
timef=time.time()
# maxV=np.max([np.max(Vfilm),-np.min(Vfilm)])
timeN_mod=timef-time0
print(timeN_mod)  
# Fonction pour mettre à jour le contenu du graphique à chaque image
# def colorK(fra,modulation,frontiere,i):
    
#     modulation

# def update(frame):
#     plt.clf()
#     plt.plot(X, Vfilm[frame,:])#, label=f'Frame {frame}')
#     plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -maxV,maxV], aspect='auto', alpha=0.35)
#     for i in range(alpha.size):
#         plt.axvline(x=alpha[i],color=colorK(frame),linestyle='--')
#     plt.colorbar()
#     # plt.plot(X, VfilmH[frame,:])
#     plt.xlabel("x (m)")
#     plt.ylabel("v (m/s)")    
#     plt.title(t[frame]*2)
#     plt.ylim([-maxV,maxV])
#     # plt.legend()

# plt.figure()
# plt.plot(Vfilm[:,500])

# total_frames = Vfilm.shape[0]
# fig, ax = plt.subplots()

# # Cree l'animationen certains pointspyten certains pointspythonen certains pointspythonen certains pointspythonen certains pointspythonen certains pointspythonen certains pointspythonhon
# animation = FuncAnimation(fig, update, frames=total_frames, interval=1)
# # Affiche l'animation
# plt.show()
# nomAnim="Film_SM_"+sourcef["Frequence"]+"Hz.mp4"
# animation.save(nomAnim, writer='ffmpeg', fps=20)



# Nt2=int(Tmax/dt)
# K=4
# timef=time.time()
# ResH,VfilmH,SfilmH=scheme1D.FD_sourceVIIHomoCara(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
# maxV=np.max([np.max(Vfilm),-np.min(Vfilm)])
# timeN_mod=timef-time0
# print(timeN_mod)  
# Fonction pour mettre à jour le contenu du graphique à chaque image
# def update(frame):
#     plt.clf()
#     plt.plot(X, Vfilm[frame,:])#, label=f'Frame {frame}')
#     plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -maxV,maxV], aspect='auto', alpha=0.35)
#     for i in range(alpha.size):
#         plt.axvline(x=alpha[i],color='gray',linestyle='--')
#     plt.colorbar()
#     # plt.plot(X, VfilmH[frame,:])
#     plt.xlabel("x (m)")
#     plt.ylabel("v (m/s)")    
#     plt.title(t[frame]*2)
#     plt.ylim([-maxV,max1V])
#     # plt.legend()

# plt.figure()
# plt.plot(Vfilm[:,0:1000],'b')
# plt.plot(VfilmH[:,0:1000],'k')


# ################ GODUNOV SPLITTING ####################################
# plt.grid()
# total_frames = Vfilm.shape[0]
# fig, ax = plt.subplots()
# timef=time.time()
ResGS,VfilmGS,SfilmGS=scheme1D.FD_sourceVII_GS_TL(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)
maxV=np.max([np.max(VfilmGS),-np.min(VfilmGS)])

nomtxtGS="WP_JK_KM_"+sourcef["Frequence"]+"_"+modulation["Freq"]+str(Nx)+"Hz.txt"
np.savetxt(nomtxtGS,VfilmGS)


energyM=np.zeros(int(Nt/pFilm))
valfct=np.zeros(int(Nt/pFilm))
for ti in range(int(Nt/pFilm)):
    energyM[ti]=dx/2*np.sum(rhox*VfilmGS[ti]*VfilmGS[ti])+1/2*np.sum(1/rhox/Celx/Celx*SfilmGS[ti]*SfilmGS[ti])
    valfct

    
    
plt.figure()
plt.plot(energyM)
plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0008,0.0008], aspect='auto', alpha=0.35)

# timeN_mod=timef-time0
# print(timeN_mod)  
# ################ STRANG SPLITTING ####################################
# plt.grid()
# total_frames = Vfilm.shape[0]
# fig, ax = plt.subplots()
# timef=time.time()
# ResSS,VfilmSS,SfilmSS=scheme1D.FD_sourceVIIHomoSS(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)
# maxV=np.max([np.max(VfilmSS),-np.min(VfilmSS)])
# timeN_mod=timef-time0
# print(timeN_mod)  
# # ###################### NO MODULATION #########25########################
# time0=time.time()
# ResNM,VfilmNM,SfilmNM=scheme1D.FD_sourceVII(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,nomodulation,pFilm)
# timeNM=time.time()-time0
# print(timeNM)  

# timef=time.time()
# ResSSC,VfilmSSC,SfilmSSC=scheme1D.FD_sourceVIIHomoSSC(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)
# timeN_mod=timef-time0
# print(timeN_mod)  


# # Cree l'animation
# animation = FuncAnimation(fig, update, frames=total_frames, interval=1)
# Affiche l'animation
# plt.show()
def updateComp(frame):
    plt.clf()
    for i in range(alpha.size-1):
        if modul.fonctionModK(pFilm*dt2p*frame,alpha[i],modulation,frontiere,i)<0.:
            plt.axvline(x=alpha[i],color='red')#,linestyle='--')
        else:
            plt.axvline(x=alpha[i],color='lime')#,linestyle='--')
    # plt.plot(X, Vfilm[frame,:],'o')#, label=f'Frame {frame}')
    # plt.plot(X, VfilmSS[frame,:],'-')
    plt.plot(X, VfilmGS[frame,:],"-")#, label=f'Frame {frame}')
    # # plt.plot(X, VfilmSS[frame,:]-Vfilm[frame,:],"-")#, label=f'Frame {frame}')
    # #plt.plot(X, VfilmSSC[frame,:]-Vfilm[frame,:],"-")
    # plt.plot(X, VfilmH[frame,:],"--")#, label=f'Frame {frame}')
    # plt.plot(X, VfilmNM[frame,:],"o")#, label=f'Frame {frame}')
    # plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0015,0.0015], aspect='auto', alpha=0.35)
    # plt.colorbar()
    # plt.plot(X, VfilmH[frame,:])
    plt.xlabel("x (m)")
    plt.ylabel("v (m/s)")
    plt.ylim([-maxV,maxV])
    plt.tight_layout()
    # plt.grid()
    # plt.legend()


plt.figure()
total_frames = VfilmGS.shape[0]
fig, ax = plt.subplots()#top=0.969,bottom=0.121,left=0.164,right=0.977,hspace=0.2,wspace=0.2)

# Cree l'animation25
animationComp = FuncAnimation(fig, updateComp, frames=total_frames, interval=2)
plt.show()
# nomAnim="WP_JK_KM_"+sourcef["Frequence"]+"_"+modulation["Freq"]+"Hz.mp4"
# animationComp.save(nomAnim, writer='ffmpeg', fps=20)

# plt.figure()
# plt.plot(Vfilm[700,:],"*")
# plt.plot(VfilmGS[700,:],"--")
# plt.plot(VfilmSS[700,:],"-")
# plt.plot(VfilmH[700,:],"--")
# plt.xlabel("x (m)")
# plt.ylabel("v (m/s)")
# plt.title("$nt=$500")
# # namefig="Figures/WP500_JK_K_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".eps"
# # plt.savefig(str(namefig), format='eps')
# plt.figure()
# plt.plot(X, Vfilm[400,:],'*')#, label=f'Frame {frame}')
# plt.plot(X, VfilmGS[400,:]-Vfilm[400,:],"-")#, label=f'Frame {frame}')
# plt.plot(X, VfilmSS[400,:]-Vfilm[400,:],"-")#, label=f'Frame {frame}')
# #plt.plot(X, VfilmSSC[800,:]-Vfilm[800,:],"-")

# print(np.max(np.abs(VfilmGS[400:900,:]-Vfilm[400:900,:])*1000))
# print(np.max(np.abs(VfilmSS[400:900,:]-Vfilm[400:900,:])*1000))
# #print(np.max(np.abs(VfilmSSC[200:900,:]-Vfilm[200:900,:])*1000))
# print(np.max(np.abs(Vfilm[400:900,:])*1000))

# plt.figure()
# plt.plot(Vfilm[900,:],"*")
# # plt.plot(VfilmGS[900,:],"-")
# # plt.plot(VfilmSS[1000,:],"--")
# # plt.plot(VfilmH[1000,:],"--")
# plt.xlabel("x (m)")
# plt.ylabel("v (m/s)")
# plt.tight_layout()
# # plt.title("$nt=$500")
# namefig="Figures/WP900_JK_K_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".eps"
# plt.savefig(str(namefig), format='eps')

def cosine_taper(signal, alpha=0.05,zp=18):
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
            tapered[j]=signal[j] * (0.5 * (1-np.cos((Size-j)*np.pi/(M+1))))
    return tapered

N=2**15
# FFTMaptot={}
# for i in range(nsou):
#     for j in range(nmod):
#         sigcor = cosine_1taper(results[f'Vfilm{i}{j}'][::,int(750)])
# # plt.figure()
# # plt.plot(Vfilm[::,int(750)])
# # plt.plot(sigcor)
# fftSig=np.fft.fftshift(np.fft.fft(sigcor,n=N))
# freq = np.fft.fftshift(np.fft.fftfreq(N,d=dt))
# plt.figure()
# plt.plot(np.abs(freq),np.abs(fftSig))



# plt.figure()
# plt.plot(Vfilm[2000,:])
# plt.plot(VfilmShft[2000,:])
#FFTMapSft=np.zeros([int(N),Nx])
freq = np.fft.fftshift(np.fft.fftfreq(N,d=dt2p))
# xpt=600
# sigsrc=cosine_taper(sce/rho[0]/2)
# FFTS=np.abs(np.fft.fftshift(np.fft.fft(sigsrc,n=int(N))))
# sigcor = cosine_taper(VfilmGS[::,int(xpt)])
# FFT0=np.abs(np.fft.fftshift(np.fft.fft(sigcor,n=int(N))))
#sigcorH = cosine_taper(VfilmH[::,int(xpt)])
#FFTH=np.abs(np.fft.fftshift(np.fft.fft(sigcorH,n=int(N))))
# sigcorGS = cosine_taper(VfilmGS[::,int(xpt)])
# FFTGS=np.abs(np.fft.fftshift(np.fft.fft(sigcorGS,n=int(N))))
#sigcorNM = cosine_taper(VfilmNM[::,int(xpt)])    
#FFTNM=np.abs(np.fft.fftshift(np.fft.fft(sigcorNM,n=int(N)))) 

# plt.figure()
# plt.plot(freq,FFT0)
# # plt.plot(freq,FFTGS)
# # plt.plot(freq,FFTH)
# # plt.plot(freq,FFTNM)
# plt.plot(freq,FFTS)
# for x in range(Nx+1):
#     sigcor = cosine_taper(Vfilm[::,x])
#     FFTMaptot[f'F{x}']=np.abs(np.fft.fftshift(np.fft.fft(sigcor,n=int(N))))
#     sigcorH = cosine_taper(VfilmH[::,x])
#     FFTMaptot[f'FH{x}']=np.abs(np.fft.fftshift(np.fft.fft(sigcorH,n=int(N))))
#     sigcorGS = cosine_taper(VfilmGS[::,x])
#     FFTMaptot[f'FGS{x}']=np.abs(np.fft.fftshift(np.fft.fft(sigcorGS,n=int(N))))
#     sigcorNM = cosine_taper(VfilmNM[::,x])  
#     FFTMaptot[f'FNM{x}']=np.abs(np.fft.fftshift(np.fft.fft(sigcorNM,n=int(N))))   
# ################ REVERSE PROPAGATION ####################################
# ResR,VfilmR,SfilmR=scheme1D.FD_sourceVII(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0R,sceR,configuration,mat,frontiere,modulation,pFilm)
# def updateR(frame):
#     plt.clf()
#     plt.plot(X, VfilmR[frame,:])#, label=f'Frame {frame}')
#     plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0008,0.0008], aspect='auto', alpha=0.35)
#     for i in range(alpha.size):
#         plt.axvline(x=alpha[i],color='gray',linestyle='--')
#     plt.colorbar()
#     # plt.plot(X, VfilmH[frame,:])
#     plt.xlabel("x (m)")
#     plt.ylabel("v (m/s)")    
#     plt.title(t[frame]*2)
#     plt.ylim([-maxV,maxV])
#     # plt.legend()

# total_frames = VfilmR.shape[0]
# figR, axR = plt.subplots()

# # Cree l'animation
# animationR = FuncAnimation(figR, updateR, frames=total_frames, interval=1)
# # Affiche l'animation
# plt.show()
# # nomAnimR="Film_SMR_"+sourcef["Frequence"]+"Hz.mp4"
# # animationR.save(nomAnimR, writer='ffmpeg', fps=20)

### Parametres adimensionnes
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

ratio=float(modulation["Ratio"])
c_r = cmef1*ratio+cmef2*(1-ratio)


fm=float(modulation["Freq"])
tau=1/fm

import numpy as np

VfilmGS=np.loadtxt(nomtxtGS)


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

# frqv,wavv,FFK=specfk(Vfilm, dx, dt2p*pFilm, int(Ntp/pFilm)+1, Nx)
# frqvS,wavvS,FFKS=specfkR(VfilmSS, dx, dt2p*pFilm, int(Ntp/pFilm)+1, Nx)
frqvG,wavvG,FFKG=specfk(Vfilm[::,::], dx, dt2p*pFilm, int(Ntp/pFilm)+1, int((Nx)))
# plt.figure()
# # plt.pcolor(wavv[::]*cmet*tau,frqv[::2]*tau,np.log(FFK[::2,::]),vmin=-12, vmax=2)
# # plt.ylim([0,2])
# # plt.xlim([0,12])
# plt.colorbar()
# plt.xlabel("k (1/m)")
# plt.ylabel("f (Hz)")
# plt.tight_layout()

plt.figure()
plt.pcolor(wavvG[::],frqvG[::],np.log(FFKG[::,::]),vmin=-12, vmax=4)
# plt.ylim([0,8*np.pi])
# plt.xlim([-0,8*np.pi])
plt.colorbar()
plt.xlabel("kc^*\tau")
plt.ylabel("f*\tau")
plt.tight_layout()
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
# # ###################### NO MODULATION ##########################################################################
# time0=time.time()
# ResNM,VfilmNM,SfilmNM=scheme1D.FD_sourceVII(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,nomodulation,pFilm)
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
# # nomAnimNM="Film_SMNM_"+sourcef["Frequence"]+"Hz.mp4"
# # animationNM.save(nomAnimNM, writer='ffmpeg', fps=20)


# #######################################
# # plt.figure()
# # plt.clf()
# # plt.plot(Vfilm[:,int(x_0R)])
# # plt.plot(VfilmNM[:,int(x_0R)])
# # plt.xlabel("t (s)")

# plt.figure()
# plt.clf()
# plt.plot(VfilmR[:,int(x_0/dx)])
# plt.plot(Vfilm[:,int(x_0R/dx)],"--")
# plt.xlabel("t (s)")
# # # plt.figure()
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
