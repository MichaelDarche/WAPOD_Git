 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:50:43 2023

@author: michael
"""

import os
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/1-Interface_Imerged_PS_modulated3/")

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../WAPOD/")

import scheme1D_TDold as scheme1D
import source
import inputReading
import modulation as modul
import analytics
NOde=2**15

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
################# Input files 
#%
configuration=inputReading.readfile('Demarrer.txt')
mat=inputReading.readfile('Milieu.txt')
frontiere=inputReading.readfile('Frontiere.txt')     
sourcef=inputReading.readfile('Source.txt')   
# sourceSft=inputReading.readfile('SourceSft.txt')
# sourceRev=inputReading.readfile('SourceRev.txt')
modulation=inputReading.readfile('Modulation.txt')
# nomodulation=inputReading.readfile('NoModulation.txt')
################ Parameters of the simulation   
X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(configuration)
fs,CFL,x_0=inputReading.source(sourcef)
# fs,CFL,x_0S=inputReading.source(sourceSft)
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
U0=scheme1D.initPointSourceProblem(configuration,frontiere,mat,sourcef)

###### Mise en oeuvre


K=4

Ntfig=int(np.floor(0.06/dt))
# Ntfig=600
Tmax=Ntfig*dt
tmax=str(Tmax)+"s"
t=np.linspace(0,Tmax,Ntfig+1)
sce=np.zeros(Ntfig+1)
for i in range(Ntfig+1):
    sce[i]=source.choice_timefct(sourcef,t[i])

# sceS=np.zeros(Nt+1)
# for i in range(Nt+1):
#     sceS[i]=source.choice_timefct(sourceSft,tim[i])
    
    
# sceR=np.zeros(Nt+1)
# for i in range(Nt+1):
#     sceR[i]=source.choice_timefct(sourceRev,tim[i])
    
    

pFilm=1
############### MODULATION DIRECT PROPAGATION ##########################################
import time
time0=time.time()
Res,Vfilm,Sfilm,CSVM,traces=scheme1D.FD_sourceVII(U0,Ntfig,Nx,K,rhox,Celx,dt,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm,rCSVM="yes",rtraces="yes")#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
time1_mod=time.time()-time0
print(time1_mod)  

SolV,SolS=analytics.analytics_1IMod(Tmax,X,sourcef,modulation,frontiere,mat,NOde,dt,Ntfig+1)

plt.figure()
plt.plot(X,Res[0,:],linestyle=None,marker='o', label='numerics')
plt.plot(X,SolV,linewidth=2,label='analytics')
plt.axvline(x=alpha[0],color='red')#,linestyle='--')
plt.figure()
plt.plot(X,Res[1,:],linestyle=None,marker='o', label='numerics')
plt.plot(X,SolS,linewidth=2,label='analytics')
plt.axvline(x=alpha[0],color='red')#,linestyle='--')



# Fonction pour mettre à jour le contenu du graphique à chaque image
# def update(frame):
#     plt.clf()
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
#     # plt.legend()

# energyM=np.zeros(int(Nt/pFilm))
# for ti in range(int(Nt/pFilm)):
#     energyM[ti]=1/2*np.sum(rhox*Vfilm[ti]*Vfilm[ti])+1/2*np.sum(1/rhox/Celx/Celx*Sfilm[ti]*Sfilm[ti])


# # import time
# # time0=time.time()
# # ResShft,VfilmShft,SfilmShft=scheme1D.FD_sourceVII(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0S,sceS,configuration,mat,frontiere,nomodulation,pFilm)#(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0S,sce,configuration,mat,frontiere,nomodulation,pFilm)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
# # timeSft_mod=time.time()-time0
# # print(timeSft_mod)  


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
# plt.figure()
# plt.plot(tim[0:-1],Vfilm[:,1500])
# plt.xlabel("t (s)")
# plt.ylabel("v (m/s)")

# plt.figure()
# plt.plot(X, Vfilm[100,:])
# plt.plot(X, Vfilm[1000,:])
# plt.plot(X, Vfilm[2000,:])
# # plt.plot(X, Vfilm[4000,:])

# plt.figure()
# plt.plot(X, Sfilm[100,:])
# plt.plot(X, Sfilm[1000,:])
# plt.plot(X, Sfilm[2000,:])
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
        

# def cosine_taper(signal, alpha=0.00,zp=16):
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
#             tapered[j]=signal[j] * (0.5 * (1-np.cos((Size-j)*np.pi/(M+1))))
#     return tapered

# N=2**16
# sigcor = cosine_taper(Vfilm[::,int(5050)])
# plt.figure()
# plt.plot(Vfilm[::,int(5050)])
# plt.plot(sigcor)
# fftSig=1/Nt*np.fft.fftshift(np.fft.fft(sigcor,n=N))
# freq = np.fft.fftshift(np.fft.fftfreq(N,d=dt))
# plt.figure()
# plt.plot(np.abs(freq),np.abs(fftSig))

# FFTMap=np.zeros([int(N),Nx])
# #FFTMapSft=np.zeros([int(N),Nx])
# for k in range(Nx):
#     sigcor = cosine_taper(Vfilm[::,k])
#     #sigcorS = cosine_taper(VfilmShft[::,k])
#     FFTMap[:,k]=np.abs(np.fft.fftshift(np.fft.fft(sigcor,n=int(N))))
    
    
# import HBM_Num as HBM

# n=4
# f=float(sourcef["Frequence"])
# F=float(modulation["FreqF"])
# dM=float(modulation["DeltaM"])
# dK=float(modulation["DeltaF"])
# K0=float(frontiere["Kn_0"])
# M0=float(frontiere["Mn_0"])
# rho0=float(mat["Rho_0"])
# c0=float(mat["Cel_0"])

# freqs,R,T=HBM.RTcoefs(n,dM,dK,rho0,c0,M0,K0,f,F)


# TN=np.abs(HBM.normV(T,n))
# TNm=np.abs(HBM.normV(T.conjugate(),n))


# plt.figure()
# plt.plot(Vfilm[2000,:])
# plt.plot(VfilmShft[2000,:])

    #FFTMapSft[:,k]=np.abs(np.fft.fftshift(np.fft.fft(sigcorS,n=int(N))))
    
# Nt=4096*2
# Tmax=Nt*dt
# tmax=str(Tmax)+"s"
# K=4
# n=10
# sigSyn=np.zeros(N)
# sigSyn2=np.zeros(N)

# plt.figure()
# plt.plot(freq, np.abs(FFTMap[:, 6001])/np.max(np.abs(FFTMap[:,6001])))
# plt.scatter(freqs, TN, c='red')
# plt.scatter(-freqs, TN, c='purple')

# plt.figure()
# plt.plot(freq, np.abs(FFTMap[:, 7500])/np.max(np.abs(FFTMap[:, 7500])))
# plt.scatter(freqs, TN, c='red')
# plt.scatter(-freqs, TNm, c='purple')

# # n=8
# ome=63
# Ome=115

# frp=[ome]

# for i in range(1,n+1):
#     # frm.append(-ome+i*Ome)
#     # frm.append(-ome-i*Ome)
#     frp.append(ome-i*Ome)
#     frp.append(ome+i*Ome)
# N=2**20



# freq = np.fft.fftshift(np.fft.fftfreq(int(N),d=dt))
# plt.figure()int(Nt)
# plt.pcolormesh(np.linspace(xmin,xmax,Nx),freq,np.real(FFTMap))
# plt.colorbar()

# plt.figure()
# plt.pcolormesh(np.linspace(xmin,xmax,Nx),np.linspace(0,Tmax,Nt),Vfilm)
# plt.colorbar()


# coefs=[float(Norm.evalf()),np.abs(float(syp.Abs(Normm.evalf())),float(Normp.evalf()))]

# plt.plot(freq,np.abs(FFTMapSft[:,750])/np.max(np.abs(FFTMapSft[:,750])))
#plt.scatter(frp, coefsTrNorm[::2], c='red')
#plt.scatter(frm, coefsTrNorm[1::2], c='red')
# plt.scatter(fr, coefsTiNorm[::2], c='green')


# plt.plot(freq,np.abs(FFTMap[:,600])/np.max(np.abs(FFTMap[:,600])))
# plt.plot(freq,np.abs(FFTMap[:,625])/np.max(np.abs(FFTMap[:,625])))
# # plt.scatter(fr,coefsT)

# plt.figure()
# plt.plot(freq,(np.abs(FFTMap[:,750])))#/np.max(np.abs(FFTMap[:,275])))
# plt.plot(freq,(np.abs(FFTMap[:,600])))#/np.max(np.abs(FFTMap[:,300])))
# plt.plot(freq,(np.abs(FFTMap[:,625])))#/np.max(np.abs(FFTMap[:,325])))
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
