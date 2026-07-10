#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:50:43 2023

@author: michael
"""

import os
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/1-Interface_Imerged_PS_modulated_complex/")

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
sourceSft=inputReading.readfile('SourceSft.txt')
sourceRev=inputReading.readfile('SourceRev.txt')
modulation=inputReading.readfile('Modulation.txt')
nomodulation=inputReading.readfile('NoModulation.txt')
################ Parameters of the simulation   
X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(configuration)
fs,CFL,x_0=inputReading.source(sourcef)
fs,CFL,x_0S=inputReading.source(sourceSft)
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
Nt=100#4096*2
Tmax=Nt*dt
tmax=str(Tmax)+"s"
K=4

tim=np.linspace(0,Tmax,Nt+1)

sce=np.zeros(Nt+1,dtype='complex')
for i in range(Nt+1):
    sce[i]=source.choice_timefct(sourcef,tim[i])

sceS=np.zeros(Nt+1,dtype='complex')
for i in range(Nt+1):
    sceS[i]=source.choice_timefct(sourceSft,tim[i])
    
    
sceR=np.zeros(Nt+1,dtype='complex')
for i in range(Nt+1):
    sceR[i]=source.choice_timefct(sourceRev,tim[i])
    
    
# from matplotlib.animation import FuncAnimation

pFilm=1
############### MODULATION DIRECT PROPAGATION ##########################################
import time
time0=time.time()
Res,Vfilm,Sfilm=scheme1D.FD_sourceVIIComplex(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
time1_mod=time.time()-time0
print(time1_mod)  
# Fonction pour mettre à jour le contenu du graphique à chaque image
def update(frame):
    plt.clf()
    plt.plot(X, Vfilm[frame,:])#, label=f'Frame {frame}')
    plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0008,0.0008], aspect='auto', alpha=0.35)
    plt.axvline(x=alpha[0],color='gray',linestyle='--')
    plt.colorbar()
    # plt.plot(X, VfilmH[frame,:])
    plt.xlabel("x (m)")
    plt.ylabel("v (m/s)")
    plt.title(ti[frame]*2)
    plt.ylim([-0.0008,0.0008])
    # plt.legend()

# energyM=np.zeros(int(Nt/pFilm))
# for ti in range(int(Nt/pFilm)):
#     energyM[ti]=1/2*np.sum(rhox*Vfilm[ti]*Vfilm[ti])+1/2*np.sum(1/rhox/Celx/Celx*Sfilm[ti]*Sfilm[ti])


# import time
# time0=time.time()
# ResShft,VfilmShft,SfilmShft=scheme1D.FD_sourceVII(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0S,sceS,configuration,mat,frontiere,nomodulation,pFilm)#(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0S,sce,configuration,mat,frontiere,nomodulation,pFilm)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
# timeSft_mod=time.time()-time0
# print(timeSft_mod)  


# total_frames = Vfilm.shape[0]
# fig, ax = plt.subplots()

# # Cree l'animation
# animation = FuncAnimation(fig, update, frames=total_frames, interval=100)
# # Affiche l'animation
# plt.show()
# nomAnim="Film_SM_"+sourcef["Frequence"]+"Hz.mp4"
# animation.save(nomAnim, writer='ffmpeg', fps=20)

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
plt.figure()
plt.plot(tim[0:-1],Vfilm[:,1500])
plt.xlabel("t (s)")
plt.ylabel("v (m/s)")

plt.figure()
plt.plot(X, Vfilm[100,:])
plt.plot(X, Vfilm[500,:])
plt.plot(X, Vfilm[2000,:])
# plt.plot(X, Vfilm[4000,:])

plt.figure()
plt.plot(X, Sfilm[100,:])
plt.plot(X, Sfilm[1000,:])
plt.plot(X, Sfilm[2000,:])
# plt.plot(X, Sfilm[4000,:])

# total_frames = VfilmNM.shape[0]
# figNM, axNM = plt.subplots()

# # Cree l'animation
# animationNM = FuncAnimation(figNM, updateNM, frames=total_frames, interval=100)
# # Affiche l'animation
# plt.show()
# nomAnimNM="Film_SMNM_"+sourcef["Frequence"]+"Hz.mp4"
# animationNM.save(nomAnimNM, writer='ffmpeg', fps=20)

#######################################
# with open("SMmodule-N1_Bruno/Vmesure-Iter2000.txt", "r") as fichier:
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

#######################################
# plt.figure()
# plt.clf()
# plt.plot(t[:-1:pFilm],Vfilm[:,int(x_0R)])
# plt.plot(t[:-1:pFilm],VfilmNM[:,int(x_0R)])
# plt.xlabel("t (s)")
# plt.ylabel("v (m/s)")

# plt.figure()
# plt.clf()
# plt.plot(t[:-1:pFilm],VfilmR[:,int(x_0)],label="Xsource<-Xrecepteur")
# plt.plot(t[:-1:pFilm],Vfilm[:,int(x_0R)],"--",label="Xsource->Xrecepteur")
# plt.xlabel("t (s)")
# plt.ylabel("v (m/s)")
# plt.legend()


# x_0B=301
# plt.figure()
# plt.clf()
# plt.plot(t[:-1:pFilm],Vfilm[:,int(x_0B)],label="Michaël")
# plt.plot(LWS0x,LWS0y,'*',label="Bruno")
# plt.xlabel("t (s)")
# plt.ylabel("v (m/s)")
# plt.legend()

# x_0B=301
# plt.figure()
# plt.clf()
# plt.plot(t[:-2:pFilm],100*(Vfilm[:-1,int(x_0B)]-LWS0y)/max(LWS0y))
# # plt.plot(LWS0x,,'*')
# plt.xlabel("t (s)")
# plt.ylabel("v (m/s)")
# plt.legend('Michaël','Bruno')

# def tFourier(signal,dt,Nfourier):
#     for n in range(Nfourier):
        

def cosine_taper(signal, alpha=0.05,zp=20):
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

N=2**20
sigcor = cosine_taper(Vfilm[::,int(750)])
plt.figure()
plt.plot(Vfilm[::,int(750)])
plt.plot(sigcor)
fftSig=np.fft.fftshift(np.fft.fft(sigcor,n=N))
freq = np.fft.fftshift(np.fft.fftfreq(N,d=dt))
plt.figure()
plt.plot(np.abs(freq),np.abs(fftSig))


omega=-8
Omega=5
sigtest=np.zeros(N)
for i in range(Nt):
    sigtest[i]=np.sin(omega*i*dt)*np.sin(Omega*i*dt)

# plt.figure()
# plt.plot(Vfilm[2000,:])
# plt.plot(VfilmShft[2000,:])
FFTMap=np.zeros([int(N),Nx])
#FFTMapSft=np.zeros([int(N),Nx])
for k in range(Nx):
    sigcor = cosine_taper(Vfilm[::,k])
    #sigcorS = cosine_taper(VfilmShft[::,k])
    FFTMap[:,k]=np.abs(np.fft.fftshift(np.fft.fft(sigcor,n=int(N))))
    #FFTMapSft[:,k]=np.abs(np.fft.fftshift(np.fft.fft(sigcorS,n=int(N))))
    
Nt=4096*4
Tmax=Nt*dt
tmax=str(Tmax)+"s"
K=4
n=8
sigSyn=np.zeros(N)
sigSyn2=np.zeros(N)

n=8
ome=42
Ome=27

frp=[ome]

for i in range(1,n+1):
    # frm.append(-ome+i*Ome)
    # frm.append(-ome-i*Ome)
    frp.append(ome-i*Ome)
    frp.append(ome+i*Ome)
N=2**21

sigSyn=np.zeros(N)
sigSyn2=np.zeros(N)
sigSyn3=np.zeros(N)
cRand=np.random.rand(2*n+1)
A=frp[0]
B=12/5*frp[5]
for nt in range(Nt):
    sigSyn[nt]=np.sin(2*np.pi*A*nt*dt)*np.sin(2*np.pi*B*nt*dt)
    sigSyn2[nt]=-1/4*(np.exp(2j*np.pi*(A+B)*nt*dt)+np.exp(-2j*np.pi*(A+B)*nt*dt)-(np.exp(2j*np.pi*(A-B)*nt*dt)+np.exp(-2j*np.pi*(A-B)*nt*dt)))
    sigSyn3[nt]=-1/2*(np.cos(2*np.pi*(A+B)*nt*dt)-np.cos(2*np.pi*(A-B)*nt*dt))
    # for i in range(2*n+1):
    #     sigSyn[nt]=sigSyn[nt]+cRand[i]*np.exp(-1j*2*np.pi*frp[i]*nt*dt)
# coefs=[float(Norm.evalf()),np.abs(float(syp.Abs(Normm.evalf())),float(Normp.evalf()))]
fftSyn=1/Nt*(np.fft.fftshift(np.abs(np.fft.fft(sigSyn,n=N))))
fftSyn2=1/Nt*(np.fft.fftshift(np.abs(np.fft.fft(sigSyn2,n=N))))
fftSyn3=1/Nt*(np.fft.fftshift(np.abs(np.fft.fft(sigSyn3,n=N))))

freq = np.fft.fftshift(np.fft.fftfreq(N,d=dt))
plt.figure()
plt.plot(tim,sigSyn[0:Nt+1],linewidth=3)
plt.plot(tim,sigSyn2[0:Nt+1],":",linewidth=2)
plt.plot(tim,sigSyn3[0:Nt+1],":",linewidth=1)
plt.figure()
plt.plot(freq,fftSyn,linewidth=2)
plt.plot(freq,fftSyn2,":",linewidth=2)
plt.plot(freq,fftSyn3,":",linewidth=1)
plt.scatter([-A-B,-A+B,A-B,A+B], [0.25,0.25,0.25,0.25], c='red')
plt.xlabel('f')
plt.ylabel('|F(u)|')

# freq = np.fft.fftshift(np.fft.fftfreq(int(N),d=dt))
# plt.figure()int(Nt)
# plt.pcolormesh(np.linspace(xmin,xmax,Nx),freq,np.real(FFTMap))
# plt.colorbar()

# plt.figure()
# plt.pcolormesh(np.linspace(xmin,xmax,Nx),np.linspace(0,Tmax,Nt),Vfilm)
# plt.colorbar()
ome=1000
Ome=20

frm=[-ome]
frp=[ome]

for i in range(1,n+1):
    frm.append(-ome+i*Ome)
    frm.append(-ome-i*Ome)
    frp.append(ome-i*Ome)
    frp.append(ome+i*Ome)

# coefs=[float(Norm.evalf()),np.abs(float(syp.Abs(Normm.evalf())),float(Normp.evalf()))]

plt.figure()
plt.plot(freq, np.abs(FFTMap[:, 1500])/np.max(np.abs(FFTMap[:, 1500])))
# plt.plot(freq,np.abs(FFTMapSft[:,750])/np.max(np.abs(FFTMapSft[:,750])))
#plt.scatter(frp, coefsTrNorm[::2], c='red')
#plt.scatter(frm, coefsTrNorm[1::2], c='red')
# plt.scatter(fr, coefsTiNorm[::2], c='green')
plt.scatter(frp, coefsTaNorm[::2], c='black')
plt.scatter(frm, coefsTaNorm[1::2], c='black')
omeS = 40
OmeS=30

frS=[omeS,omeS-OmeS,omeS+OmeS,ome-2*Ome,ome+2*Ome,ome-3*Ome,ome+3*Ome,ome-4*Ome,ome+4*Ome,ome-5*Ome,ome+5*Ome,ome-6*Ome,ome+6*Ome,ome-7*Ome,ome+7*Ome,ome-8*Ome,ome+8*Ome,ome-9*Ome,ome+9*Ome,ome-10*Ome,ome+10*Ome]
# coefs=[float(Norm.evalf()),np.abs(float(syp.Abs(Normm.evalf())),float(Normp.evalf()))]

plt.figure()
#plt.plot(freq,np.abs(FFTMap[:,750])/np.max(np.abs(FFTMap[:,750])))
plt.plot(freq,np.abs(FFTMapSft[:,750])/np.max(np.abs(FFTMapSft[:,750])))
plt.scatter(frS,coefsTrNormS,c='red')
plt.scatter(frS,coefsTiNormS,c='green')
plt.scatter(frS,coefsTaNormS,c='black')


plt.plot(freq,np.abs(FFTMap[:,600])/np.max(np.abs(FFTMap[:,600])))
plt.plot(freq,np.abs(FFTMap[:,625])/np.max(np.abs(FFTMap[:,625])))
# plt.scatter(fr,coefsT)

plt.figure()
plt.plot(freq,(np.abs(FFTMap[:,750])))#/np.max(np.abs(FFTMap[:,275])))
plt.plot(freq,(np.abs(FFTMap[:,600])))#/np.max(np.abs(FFTMap[:,300])))
plt.plot(freq,(np.abs(FFTMap[:,625])))#/np.max(np.abs(FFTMap[:,325])))
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
