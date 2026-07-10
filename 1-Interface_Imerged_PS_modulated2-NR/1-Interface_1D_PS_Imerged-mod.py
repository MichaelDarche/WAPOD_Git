#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:50:43 2023

@author: michael
"""

import os
#os.chdir("/home/darche/Documents/ModulatedMetamaterials/WAPRO_TD/1-Interface_Imerged_PS_modulated2/")
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/1-Interface_Imerged_PS_modulated2-NR/")

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../WAPOD/")

import scheme1D_TDold as scheme1D
import source
import inputReading
import modulation as modul
import figures as ffg

pFilm=1
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
################# Input files 
Ns=2
Nm=101
fmod=np.linspace(0,400,Nm)
err=np.zeros([Ns,Nm])
# nameSource="Source"+".txt"
# nameSourceR="SourceRev"+".txt"
for nsource in range(1,Ns):
    for nmodu in range(Nm):      
        nameSource="Source"+str(nsource)+".txt"
        nameSourceR="SourceRev"+str(nsource)+".txt"
        nameMod="Modulation.txt"#+str(nmodu)+".txt"
        configuration=inputReading.readfile('Demarrer.txt')
        # configuration2=inputReading.readfile('Demarrer2.txt')
        # configuration4=inputReading.readfile('Demarrer4.txt')
        mat=inputReading.readfile('Milieu.txt')
        frontiere=inputReading.readfile('Frontiere.txt')      
        sourcef=inputReading.readfile(nameSource)   
        # sourceSft=inputReading.readfile('SourceSft.txt')
        sourceRev=inputReading.readfile(nameSourceR)
        lineM="FreqM  :  "+str(fmod[nmodu])+"\n"
        lineK="FreqF  :  "+str(fmod[nmodu])+"\n"
        with open('Modulation.txt','r') as f:
            lines=f.readlines()
        lines[5]=lineM
        lines[6]=lineK
        with open('Modulation.txt','w') as f:
            f.writelines(lines)
        modulation=inputReading.readfile(nameMod)
        # nomodulation=inputReading.readfile('NoModulation.txt')
        ################ Parameters of the simulation   
        X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(configuration)
        fs,CFL,x_0=inputReading.source(sourcef)
        #fs,CFL,x_0S=inputReading.source(sourceSft)
        fs,CFL,x_0R=inputReading.source(sourceRev)
        Nmat,rho,cm=inputReading.material(mat)
        alpha=inputReading.frontiere(frontiere)
        rhox,Celx=inputReading.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
        dt=min(scheme1D.timeStep(CFL,dx,cm))
#%
        U0=scheme1D.initPointSourceProblem(configuration,frontiere,mat,sourcef)
        Tmax=0.07
        Nt=int(Tmax/dt)

        # Tmax=Nt*dt
        # tmax=str(Tmax)+"s"
        K=4

        tim=np.linspace(0,Tmax,Nt+1)

        sce=np.zeros(Nt+1)
        for i in range(Nt+1):
            sce[i]=source.choice_timefct(sourcef,tim[i])
        sceR=np.zeros(Nt+1)
        for i in range(Nt+1):
            sceR[i]=source.choice_timefct(sourceRev,tim[i])
            
        Res,Vfilm,Sfilm,CSVM=scheme1D.FD_sourceVII(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm,rCSVM="yes")#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
        ResR,VfilmR,SfilmR=scheme1D.FD_sourceVII(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0R,sceR,configuration,mat,frontiere,modulation,pFilm)

        xz=int(x_0/dx)
        xzR=int(x_0R/dx)
        
        diff=Vfilm[:,xzR]-VfilmR[:,xz]


        err[nsource-1,nmodu]=np.sum(diff**2)


        # plt.figure()
        # plt.plot(Vfilm[:,xzR],linewidth=3)
        # plt.plot(VfilmR[:,xz],"--",linewidth=3)
        # plt.xlabel("$t$ (s)",fontsize=16)
        # plt.ylabel("$v$ (m/s)",fontsize=16)
        # ffg.format_fig("$t$ (s)", "$v$ (m/s)")
        # plt.tight_layout()
        # nomFig="NR_"+sourcef["Frequence"]+"_"+modulation["FreqM"]+"Hz.pdf"
        # # animationComp.save(nomAnim, writer='ffmpeg', fps=20)
        # plt.savefig(nomFig, format="pdf")

fsou=65.

colorplot=['black','blue','red','purple','green']
plt.figure()
for i in range(nsource+1):
    for j in range(nmodu+1):
        plt.scatter(fmod[j],np.sqrt(err[i,j]),color=colorplot[i])
plt.figure()
for i in range(nsource+1):
    for j in range(nmodu+1):
        plt.scatter(fmod[j],err[i,j],color=colorplot[i])
    

plt.figure()
for i in range(Ns):
    # f=fsou[i]
    # fmodm=fmod/f
    plt.plot(fmod,np.sqrt(err[i]))
plt.xlabel(r'$f_m$ (Hz)')
plt.ylabel(r'$\varepsilon_v=\sqrt{\sum(v_d-v_r)^2}$ (m/s)')

# lambdaOnde=cm/fs
# xsmin=cm*tshift-lambdaOnde
# xsmax=cm*tshift
#%

#% Initial conditions


###### Mise en oeuvre 




# sceS=np.zeros(Nt+1)
# for i in range(Nt+1):
#     sceS[i]=source.choice_timefct(sourceSft,tim[i])
    
    





# Fonction pour mettre à jour le contenu du graphique à chaque image
# def update(frame):
#     plt.clf()
#     plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0008,0.0008], aspect='auto', alpha=0.35)
#     plt.plot(X[::], Vfilm[frame,:])#, label=f'Frame {frame}')
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


# # import time
# # time0=time.time()
# # ResShft,VfilmShft,SfilmShft=scheme1D.FD_sourceVII(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0S,sceS,configuration,mat,frontiere,nomodulation,pFilm)#(U0,Nt,Nx,K,rhox,Celx,dt,dx,x_0S,sce,configuration,mat,frontiere,nomodulation,pFilm)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
# # timeSft_mod=time.time()-time0
# # print(timeSft_mod)  


# total_frames = Vfilm.shape[0]
# fig, ax = plt.subplots()

# # Cree l'animation
# animation = FuncAnimation(fig, update, frames=total_frames, interval=10)
# # Affiche l'animation
# plt.show()





# nomAnim="Film_SM_"+sourcef["Frequence"]+"Hz.mp4"
# animation.save(nomAnim, writer='ffmpeg', fps=20)
# energyM=np.zeros(int(Nt/pFilm))
# valfct=np.zeros(int(Nt/pFilm))
# for ti in range(int(Nt/pFilm)):
#     energyM[ti]=dx/2*np.sum(rhox*Vfilm[ti]*Vfilm[ti])+1/2*np.sum(1/rhox/Celx/Celx*Sfilm[ti]*Sfilm[ti])
#     valfct

# U0=scheme1D.initPointSourceProblem(configuration2,frontiere,mat,sourcef)   
# X2,xmin,xmax,Nx2,dx2,ESIM,Npt2=inputReading.geometric(configuration2)
# rhox,Celx=inputReading.affectMaterials(alpha,rho,cm,xmin,xmax,Nx2,dx2)

# tim2=np.linspace(0,Tmax,2*Nt+1)

# sce=np.zeros(2*Nt+1)
# for i in range(2*Nt+1):
#     sce[i]=source.choice_timefct(sourcef,tim2[i])
# ############### MODULATION DIRECT PROPAGATION ##########################################
# import time
# time0=time.time()
# Res2,Vfilm2,Sfilm2,CSVM2=scheme1D.FD_sourceVII(U0,2*Nt,Nx2,K,rhox,Celx,dt/2,dx2,x_0,sce,configuration2,mat,frontiere,modulation,2*pFilm,rCSVM="yes")#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
# time1_mod=time.time()-time0
# print(time1_mod)  
# Fonction pour mettre à jour le contenu du graphique à chaque image


# U0=scheme1D.initPointSourceProblem(configuration4,frontiere,mat,sourcef) 
# X4,xmin,xmax,Nx4,dx4,ESIM,Npt4=inputReading.geometric(configuration4)
# rhox,Celx=inputReading.affectMaterials(alpha,rho,cm,xmin,xmax,Nx4,dx4)
# tim4=np.linspace(0,Tmax,4*Nt+1)

# sce=np.zeros(4*Nt+1)
# for i in range(4*Nt+1):
#     sce[i]=source.choice_timefct(sourcef,tim4[i])
############### MODULATION DIRECT PROPAGATION ##########################################
# import time
# time0=time.time()
# Res4,Vfilm4,Sfilm4,CSVM4=scheme1D.FD_sourceVII(U0,4*Nt,Nx4,K,rhox,Celx,dt/4,dx4,x_0,sce,configuration4,mat,frontiere,modulation,4*pFilm,rCSVM="yes")#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
# time1_mod=time.time()-time0
# print(time1_mod)   
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
# rhox,Celx=inputReading.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
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


# x_0B=301p
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
    

# Nalpha=int(alpha[0]/dx/4)
# intCSV=np.zeros(int(Nt/pFilm))
# # intCSV2=np.zeros(int(Nt/pFilm))
# # intCSV4=np.zeros(int(Nt/pFilm))
# for i in range(int(Nt/pFilm)):
#     intCSV[i]=intCSV[i-1]+Vfilm[i,Nalpha+1]-Vfilm[i,Nalpha-1]
#     # intCSV2[i]=intCSV2[i-1]+Vfilm2[i,2*(Nalpha+1)]-Vfilm2[i,2*(Nalpha-1)]
#     # intCSV4[i]=intCSV4[i-1]+Vfilm4[i,4*(Nalpha+1)]-Vfilm4[i,4*(Nalpha-1)]
# CSV=2*pFilm*dt*intCSV/(Sfilm[:,Nalpha+1]+Sfilm[:,Nalpha-1])
# # CSV2=2*pFilm*dt*intCSV2/(Sfilm2[:,2*(Nalpha+1)]+Sfilm2[:,2*(Nalpha-1)])
# # CSV4=2*pFilm*dt*intCSV4/(Sfilm4[:,4*(Nalpha+1)]+Sfilm4[:,4*(Nalpha-1)])

# if float(frontiere["Kn_0"])==0:
#     K0=np.infty
# else:
#     K0=float(frontiere["Kn_0"])
# frK=float(modulation["FreqF"])
# delta=float(modulation["DeltaF"])

# modula=np.zeros(int(Nt/pFilm))
# dmodula=np.zeros(int(Nt/pFilm))
# for i in range(int(Nt/pFilm)):
#     modula[i]=1/K0*(1+delta*np.sin(2*np.pi*frK*tim[pFilm*i]))
#     dmodula[i]=2*np.pi*frK/K0*(delta*np.cos(2*np.pi*frK*tim[pFilm*i]))
    



# # plt.figure()
# # plt.plot(tim[:-1:pFilm],CSV)
# # plt.plot(tim[:-1:pFilm],CSV2)
# # # plt.plot(tim[:-1:pFilm],CSV4)
# # plt.plot(tim[:-1:pFilm],modula)

# Nalpha=int(alpha[0]/dx/4)
# intCSS=np.zeros(int(Nt/pFilm))
# # intCSS2=np.zeros(int(Nt/pFilm))
# # intCSS4=np.zeros(int(Nt/pFilm))
# for i in range(int(Nt/pFilm)):
#     intCSS[i]=intCSS[i-1]+Sfilm[i,Nalpha+1]-Sfilm[i,Nalpha-1]
#     # intCSS2[i]=intCSS2[i-1]+Sfilm2[i,2*(Nalpha+1)]-Sfilm2[i,2*(Nalpha-1)]
#     # intCSS4[i]=intCSS4[i-1]+Sfilm4[i,4*(Nalpha+1)]-Sfilm4[i,4*(Nalpha-1)]
# CSS=2*pFilm*dt*intCSS/(Vfilm[:,Nalpha+1]+Vfilm[:,Nalpha-1])
# # CSS2=2*pFilm*dt*intCSS2/(Vfilm2[:,2*(Nalpha+1)]+Vfilm2[:,2*(Nalpha-1)])
# # CSS4=2*pFilm*dt*intCSS4/(Vfilm4[:,4*(Nalpha+1)]+Vfilm4[:,4*(Nalpha-1)])

# M0=float(frontiere["Mn_0"])
# frM=float(modulation["FreqM"])
# deltaM=float(modulation["DeltaM"])

# modulb=np.zeros(int(Nt/pFilm))
# dmodulb=np.zeros(int(Nt/pFilm))
# for i in range(int(Nt/pFilm)):
#     modulb[i]=M0*(1+deltaM*np.sin(2*np.pi*frM*tim[pFilm*i]))
#     dmodulb[i]=2*np.pi*frM*M0*(deltaM*np.cos(2*np.pi*frM*tim[pFilm*i]))



# # plt.figure()
# # plt.plot(tim[:-1:pFilm],CSS)
# # plt.plot(tim[:-1:pFilm],CSS2)
# # # plt.plot(tim[:-1:pFilm],CSS4)
# # plt.plot(tim[:-1:pFilm],modulb)      

  
# energyV=np.zeros(int(Nt/pFilm))
# energyI=np.zeros(int(Nt/pFilm))
# for ti in range(int(Nt/pFilm)):
#     energyV[ti]=dx/2*np.sum(rhox*Vfilm[ti]*Vfilm[ti])+dx/2*np.sum(1/rhox/Celx/Celx*Sfilm[ti]*Sfilm[ti])
#     energyI[ti]=1/2*modulb[ti]*CSVM[ti,2]*CSVM[ti,2]+1/2*modula[ti]*CSVM[ti,3]*CSVM[ti,3]    
    
# energyM=energyV+energyI
# plt.figure()
# plt.plot(energyV)
# plt.plot(max(energyV)-energyI)
# plt.plot(energyM)
# plt.imshow([modula], cmap='gray', extent=[0, int(Nt/pFilm), 0,max(energyM)], aspect='auto', alpha=0.35)

# energyV2=np.zeros(int(Nt/pFilm))
# energyI2=np.zeros(int(Nt/pFilm))
# for ti in range(int(Nt/pFilm)):
#     energyV2[ti]=1/2*np.sum(rhox*Vfilm2[ti,::2]*Vfilm2[ti,::2])+1/2*np.sum(1/rhox/Celx/Celx*Sfilm2[ti,::2]*Sfilm2[ti,::2])
#     energyI2[ti]=1/2*modulb[ti]*CSVM2[ti::2,2]*CSVM2[ti::2,2]+1/2*modula[ti]*CSVM2[ti::2,3]*CSVM2[ti::2,3]    
    
# energyM2=energyV2+energyI2
# plt.figure()
# plt.plot(energyV2)
# plt.plot(energyI2)
# plt.plot(energyM2)
# plt.imshow([modula], cmap='gray', extent=[0, int(Nt/pFilm), 0,max(energyM2)], aspect='auto', alpha=0.35)


# dVMV=np.zeros(Nt)
# dVMS=np.zeros(Nt)
# for i in range(1,Nt-1):
#     dVMV[i]=(CSVM[i+1,2]-CSVM[i-1,2])/dt/2
#     dVMS[i]=(CSVM[i+1,3]-CSVM[i-1,3])/dt/2


# plt.figure()
# plt.plot(CSVM[:,0])
# plt.plot(modula*dVMS[:]+dmodula*CSVM[:,3],'--')

# plt.figure()
# plt.plot(CSVM[:,1])
# plt.plot(modulb*dVMV[:]+dmodulb*CSVM[:,2],'--')


xmes=225
Ncap=int(xmes/dx)

xr=185
Nr=int(xr/dx)


# def cosine_taper(signal, alpha=0.0,zp=18):
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

# Zp=19
# Nzp=2**Zp
# sigcor = cosine_taper(Vfilm[::,Ncap],zp=Zp)
# # plt.figure()
# # plt.plot(Vfilm[::,int(2200)])
# # plt.plot(sigcor)
# fftSig=1/Nt*np.fft.fftshift(np.fft.fft(sigcor,n=Nzp))
# freq = np.fft.fftshift(np.fft.fftfreq(Nzp,d=dt*pFilm))
# # plt.figure()
# # plt.plot(np.abs(freq),np.abs(fftSig))

# matH=inputReading.readfile('MilieuH.txt')
# frontiereH=inputReading.readfile('FrontiereH.txt')  
# Nmat0,rho0,cm0=inputReading.material(matH)
# alpha0=inputReading.frontiere(frontiereH)
# rhox0,Celx0=inputReading.affectMaterials(alpha0,rho0,cm0,xmin,xmax,Nx,dx)
# dt=min(scheme1D.timeStep(CFL,dx,cm0))


# import time
# time0=time.time()
# # Res0,Vfilm0,Sfilm0,CSVM0=scheme1D.FD_sourceVII(U0,Nt,Nx,K,rhox0,Celx0,dt,dx,x_0,sce,configuration,matH,frontiereH,modulation,pFilm,rCSVM="yes")#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
# time1_mod=time.time()-time0
# print(time1_mod)  

# sigcorR = cosine_taper(Vfilm[::,Nr]-Vfilm0[::,Nr],zp=Zp)
# plt.figure()
# plt.plot(Vfilm[::,Nr])
# plt.plot(sigcorR)
# fftSigR=1/Nt*np.fft.fftshift(np.fft.fft(sigcorR,n=Nzp))
# freqR = np.fft.fftshift(np.fft.fftfreq(Nzp,d=dt*pFilm))
# plt.figure()
# plt.plot(np.abs(freqR),np.abs(fftSigR))

# # FFTMap=np.zeros([int(N),Nx])
# # #FFTMapSft=np.zeros([int(N),Nx])
# # for k in range(Nx):
# #     sigcor = cosine_taper(Vfilm[::,k])
# #     #sigcorS = cosine_taper(VfilmShft[::,k])
# #     FFTMap[:,k]=np.abs(np.fft.fftshift(np.fft.fft(sigcor,n=int(N))))
    

# import HBM_Num_Vform as HBM

# n=12
# f=float(sourcef["Frequence"])
# F=float(modulation["FreqF"])
# dM=float(modulation["DeltaM"])
# dK=float(modulation["DeltaF"])
# K0=float(frontiere["Kn_0"])
# M0=float(frontiere["Mn_0"])
# rho0=float(mat["Rho_0"])
# c0=float(mat["Cel_0"])


# freqs,R,T=HBM.RTcoefs(n,dM,dK,rho0,c0,M0,K0,f,F)

# n2=50
# APhi,BPhi=HBM.sysPhi(n2,dM,rho0,c0,M0,f,F)
# APsi,BPsi=HBM.sysPsi(n2,dK,rho0,c0,K0,f,F)
# freqs2,R2,T2=HBM.RTcoefs(n2,dM,dK,rho0,c0,M0,K0,f,F)

# # n3=30
# # freqs3,R3,T3=HBM.RTcoefs(n3,dM,dK,rho0,c0,M0,K0,f,F)
# def colorpoint(freq):
#     c=[]
#     for fr in freq:
#         if fr<0:
#             c.append('orchid')
#         else:
#             c.append('limegreen')
#     return c

# plt.figure()
# plt.axvline(f,linestyle="-",c='deeppink')
# plt.plot(np.abs(freq),np.abs(fftSig)/np.max(np.abs(fftSig)))
# plt.scatter(np.abs(freqs), np.abs(T)/np.max(np.abs(T)), c=colorpoint(freqs),s=60,label='$N=12$')
# plt.scatter(np.abs(freqs2), np.abs(T2)/np.max(np.abs(T2)), c='blue',marker='*',label='$N=5$',zorder=3)
# plt.scatter(-freqs, np.abs(T)/np.max(np.abs(T)), c='grey',label='$N=12$')
# plt.scatter(-freqs2, np.abs(T2)/np.max(np.abs(T2)), c='black',marker='*',label='$N=5$')
# plt.scatter(freqs3, np.abs(T3)/np.max(np.abs(T3)), c='black',marker='*',label='$N=30$')
# plt.xlim([0,700])
# plt.xlabel("$f$ (Hz)",fontsize=16)
# plt.ylabel("$|T_k|$",fontsize=16)
# #plt.legend(fontsize=16)
# plt.tight_layout()


# plt.figure()
# plt.plot(np.abs(freqR),np.abs(fftSigR)/np.max(np.abs(fftSigR)))
# plt.axvline(f,linestyle="--",c='black')
# plt.scatter(np.abs(freqs), np.abs(R)/np.max(np.abs(R)), c=colorpoint(freqs),label='$N=12$')
# # plt.scatter(np.abs(freqs2), np.abs(R2)/np.max(np.abs(R2)), c='blue',marker='*',label='$N=5$',zorder=3)
# # plt.scatter(-freqs, np.abs(T)/np.max(np.abs(T)), c='grey',label='$N=12$')
# # plt.scatter(-freqs2, np.abs(T2)/np.max(np.abs(T2)), c='black',marker='*',label='$N=5$')
# # plt.scatter(freqs3, np.abs(T3)/np.max(np.abs(T3)), c='black',marker='*',label='$N=30$')
# plt.xlim([0,700])
# plt.xlabel("$f$ (Hz)",fontsize=16)
# plt.ylabel("$|R_k|$",fontsize=16)
# #plt.legend(fontsize=16)
# plt.tight_layout()



    
# plt.figure()
# plt.axvline(f,linestyle="-",c="red")
# plt.plot(freq,np.abs(fftSig)*4*c0)
# plt.scatter(np.abs(freqs), np.abs(T), c=colorpoint(freqs),label='$N=12$',s=60)
# # plt.scatter(np.abs(freqs2), np.abs(T2), c='blue',marker='*',label='$N=5$')
# # plt.scatter(-freqs, np.abs(T), c='grey',label='$N=12$')
# # plt.scatter(-freqs2, np.abs(T2), c='black',marker='*',label='$N=5$')
# # plt.scatter(freqs3, np.abs(T3)/np.max(np.abs(T3)), c='black',marker='*',label='$N=30$')
# plt.xlim([0,700])
# plt.xlabel("$f$ (Hz)",fontsize=16)
# plt.ylabel("$|T_k|$",fontsize=16)
# # plt.legend(fontsize=16)
# plt.tight_layout()

# plt.figure()
# plt.axvline(f,linestyle="-",c="red")
# plt.plot(freq,np.abs(fftSigR)*4*c0)
# plt.scatter(np.abs(freqs), np.abs(R), c=colorpoint(freqs),label='$N=12$',s=60)
# # plt.scatter(np.abs(freqs2), np.abs(T2), c='blue',marker='*',label='$N=5$')
# # plt.scatter(-freqs, np.abs(T), c='grey',label='$N=12$')
# # plt.scatter(-freqs2, np.abs(T2), c='black',marker='*',label='$N=5$')
# # plt.scatter(freqs3, np.abs(T3)/np.max(np.abs(T3)), c='black',marker='*',label='$N=30$')
# plt.xlim([0,700])
# plt.xlabel("$f$ (Hz)",fontsize=16)
# plt.ylabel("$|R_k|$",fontsize=16)
# # plt.legend(fontsize=16)
# plt.tight_layout()

# plt.figure()
# # plt.semilogy()
# # plt.plot(freq,np.abs(fftSig)/np.max(np.abs(fftSig)))
# plt.axvline(f,linestyle="--")
# plt.scatter(freqs, np.abs(T), c='red',label='$N=12$')
# plt.scatter(freqs2, np.abs(T2), c='blue',marker='*',label='$N=5$')
# plt.xlabel("$f$ (Hz)",fontsize=16)
# plt.ylabel("$T_k$",fontsize=16)
# plt.legend(fontsize=16)

# plt.figure()
# # plt.semilogy()
# #plt.plot(freq,np.abs(fftSig)/np.max(np.abs(fftSig)))
# plt.axvline(f,linestyle="--")
# plt.scatter(freqs, np.abs(R), c='red',label='$N=12$')
# plt.scatter(freqs2, np.abs(R2), c='blue',marker='*',label='$N=5$')
# plt.xlabel("$f$ (Hz)",fontsize=16)
# plt.ylabel("$R_k$",fontsize=16)
# plt.legend(fontsize=16)

# TN=np.abs(HBM.normV(T,n))
# TNm=np.abs(HBM.normV(T.conjugate(),n))



# VTm=np.zeros(Nt,dtype=complex)
# VTp=np.zeros(Nt,dtype=complex)
# for j in range(2*n+1):
#     for tt in range(Nt):
#         VTm[tt]=VTm[tt]+np.conj(T[j])*np.exp(-1j*2*np.pi*freqs[j]*(tt*dt-1/c0*(xmes-x_0)))
#         VTp[tt]=VTp[tt]+T[j]*np.exp(1j*2*np.pi*freqs[j]*(tt*dt-1/c0*(xmes-x_0)))
# VT=1/2/2j/c0*(-VTm+VTp)#nO.subs(param).subs(syp.Indexed('tm', j+1),syp.exp(-syp.I*(omega-j*Omega)/c0*(x-x0))*coefsT[4*j+3]).subs(syp.Indexed('tp', j+1),syp.exp(-syp.I*(omega+j*Omega)/c0*(x-x0))*coefsT[4*j+5])

# plt.figure()
# plt.plot(tim[:-1],VT)
# plt.plot(tim[:-1:pFilm],Vfilm[:,Ncap])


# VRm=np.zeros(Nt,dtype=complex)
# VRp=np.zeros(Nt,dtype=complex)
# Vin=np.zeros(Nt,dtype=complex)
# for tt in range(Nt):
#     for j in range(2*n+1):
#         VRm[tt]=VRm[tt]+np.conj(R[j])*np.exp(-1j*2*np.pi*freqs[j]*(tt*dt+1/c0*(xc-x_0+2*alpha0)))
#         VRp[tt]=VRp[tt]+R[j]*np.exp(1j*2*np.pi*freqs[j]*(tt*dt+1/c0*(xc-x_0+2*alpha0)))
#     Vin[tt]=Vin[tt]-1/2j*np.exp(-1j*2*np.pi*freqs[n]*(tt*dt-1/c0*(xc-x_0)))+1/2j*np.exp(1j*2*np.pi*freqs[n]*(tt*dt-1/c0*(xc-x_0)))
# VR=1/2j/Nt*(+VRp-VRm)+Vin/2/cm0#nO.subs(param).subs(syp.Indexed('tm', j+1),syp.exp(-syp.I*(omega-j*Omega)/c0*(x-x0))*coefsT[4*j+3]).subs(syp.Indexed('tp', j+1),syp.exp(-syp.I*(omega+j*Omega)/c0*(x-x0))*coefsT[4*j+5])

# plt.figure()
# plt.plot(tim[:-1],VR)
# plt.plot(tim[:-1],Vin/2/cm0)
# [plt.plot(tim[:-1:pFilm],Vfilm[:,Ncap])]


# import figures as ffg

# Ucalc=ffg.VtoU(Vfilm, dt)
# plt.figure()
# plt.plot(Ucalc)

# plt.plot(Vfilm0[:,xmes])
# plt.plot(1/2j/cm0*(np.exp(1j*2*np.pi*freqs[n]*(tim[:]*dt-1/c0*(xmes-x_0)))-np.exp(-1j*2*np.pi*freqs[n]*(tim[:]*dt-1/c0*(xmes-x_0)))))

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
# plt.plot(freq, np.abs(FFTMap[:, 5001])/np.max(np.abs(FFTMap[:,5001])))
# plt.scatter(freqs, TN, c='red')
# plt.scatter(-freqs, TN, c='purple')

# plt.figure()
# plt.plot(freq, np.abs(FFTMap[:, 7500])/np.max(np.abs(FFTMap[:, 7500])))
# plt.scatter(freqs, TN, c='red')
# plt.scatter(-freqs, TNm, c='purple')

# plt.figure()
# plt.plot(freq, np.abs(FFTMap[:, 5001])/np.max(np.abs(FFTMap[:,5001])))
# plt.plot(freq, np.abs(FFTMap[:, 7500])/np.max(np.abs(FFTMap[:, 7500])))



# n=8
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
