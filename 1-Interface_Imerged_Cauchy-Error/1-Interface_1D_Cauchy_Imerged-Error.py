#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:50:43 2023

@author: michael
"""

import os
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/1-Interface_Imerged_Cauchy-Error/")

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../WAPOD/")

import scheme1D_TDold as scheme1D
import source
import inputReading



################# Input files 
#%
def errL2(A,B):
    err=0
    B2=0
    for l in range(max(B.shape)):
        eps=A[l]-B[l]
        err=err+eps**2
        B2=B2+B[l]**2
    return err,err/B2

ndx=7
err=np.zeros([4,ndx])
errIm=np.zeros([4,ndx])
errR=np.zeros([4,ndx])
errRIm=np.zeros([4,ndx])
data=np.zeros([4,ndx])
k0=2**((ndx-1)/2)
for k in range(ndx):
    nDem='Demarrer'+str(k)+'.txt'
    configuration=inputReading.readfile(nDem)
    mat=inputReading.readfile('Milieu.txt')
    frontiere=inputReading.readfile('Frontiere.txt')     
    sourcef=inputReading.readfile('Source.txt')   
    modulation=inputReading.readfile('Modulation.txt')
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
    Nt=int((Nx-1)/2)
    Tmax=Nt*dt
    tmax=str(Tmax)+"s"
    THV0y=np.zeros(Nx)
    THVTMy=np.zeros(Nx)
    THS0y=np.zeros(Nx)
    THSTMy=np.zeros(Nx)
    R01=(rho[1]*cm[1]-rho[0]*cm[0])/(rho[1]*cm[1]+rho[0]*cm[0])
    T01=2*rho[0]*cm[1]/(rho[1]*cm[1]+rho[0]*cm[0])
    Nxalpha=int(alpha[0]/dx)
    WL=max(cm/fs)
    data[0,k]=Tmax
    data[1,k]=dx
    data[2,k]=dt
    data[3,k]=WL/dx
    for x in range(Nxalpha):
        f0=source.choice_timefct(sourcef,tshift-x*dx/cm[0])
        fr=source.choice_timefct(sourcef,tshift+x*dx/cm[0]-2*alpha[0]/cm[0])
        THV0y[x]=1/cm[0]*f0-R01/cm[0]*fr
        THS0y[x]=-rho[0]*f0-R01*rho[0]*fr
        fT=source.choice_timefct(sourcef,Tmax+tshift-x*dx/cm[0])
        frT=source.choice_timefct(sourcef,Tmax+tshift+x*dx/cm[0]-2*alpha[0]/cm[0])
        THVTMy[x]=1/cm[0]*fT-R01/cm[0]*frT
        THSTMy[x]=-rho[0]*fT-R01*rho[0]*frT
    for x in range(Nxalpha,Nx):
        ft=source.choice_timefct(sourcef,tshift-x*dx/cm[1]+(1/cm[1]-1/cm[0])*alpha[0])
        ftT=source.choice_timefct(sourcef,Tmax+tshift-x*dx/cm[1]+(1/cm[1]-1/cm[0])*alpha[0])
        THV0y[x]=1/cm[1]*T01*ft
        THS0y[x]=-rho[1]*T01*ft
        THVTMy[x]=1/cm[1]*T01*ftT
        THSTMy[x]=-rho[1]*T01*ftT
    for i in range(2):
        K=2+2*i
        Res=scheme1D.FD_cauchy(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
        # plt.figure()
        # plt.plot(X,Res[0,:])
        # plt.plot(X,THVTMy,'--')
        # plt.xlabel("$x$ (m)")
        # plt.ylabel("$v_{th}$ (m/s)")
        err[0+2*i,k],errR[0+2*i,k]=np.sqrt(dx)*np.sqrt(errL2(Res[0,:], THVTMy))
        err[1+2*i,k],errR[1+2*i,k]=np.sqrt(dx)*np.sqrt(errL2(Res[1,:], THSTMy))
        ResIm=scheme1D.FD_cauchyII(U0,Nt,Nx,K,rhox,Celx,dt,dx,configuration,mat,frontiere,modulation,"no")
        # plt.figure()
        # plt.plot(X,Res[0,:])
        # plt.plot(X,THVTMy,'--')
        # plt.xlabel("$x$ (m)")
        # plt.ylabel("$v_{th}$ (m/s)")
        errIm[0+2*i,k],errRIm[0+2*i,k]=np.sqrt(dx)*np.sqrt(errL2(ResIm[0,:], THVTMy))
        errIm[1+2*i,k],errRIm[1+2*i,k]=np.sqrt(dx)*np.sqrt(errL2(ResIm[1,:], THSTMy))
        
plt.figure()
plt.loglog(data[1,:],err[0,:],'*--',label="Law-Wendroff")
plt.loglog(data[1,:],err[2,:],'*--',label="ADER4")
plt.loglog(data[1,:],errIm[0,:],'*-',label="Law-Wendroff+II")
plt.loglog(data[1,:],errIm[2,:],'*-',label="ADER4+II")
plt.loglog(data[1,:],data[1,:]**4*10e-8)
plt.xlabel("$\Delta x$ (m)")
plt.ylabel("$\epsilon_v$ (m/s)")
plt.title("$\epsilon=\sqrt{\sum (U-U_{th})^2}$")
# plt.savefig('Figures/ErrV_Npt3_ESIM3.eps', format='eps')
plt.legend()
plt.figure()
plt.loglog(data[1,:],err[1,:],'*--',label="Law-Wendroff")
plt.loglog(data[1,:],err[3,:],'*--',label="ADER4")
plt.loglog(data[1,:],errIm[1,:],'*-',label="Law-Wendroff+II")
plt.loglog(data[1,:],errIm[3,:],'*-',label="ADER4+II")
plt.xlabel("$\Delta x$ (m)")
plt.ylabel("$\epsilon_\sigma$ (Pa)")
plt.title("$\epsilon=\sqrt{\sum (U-U_{th})^2}$")
plt.legend()
# plt.savefig('Figures/ErrS_Npt3_ESIM3.eps', format='eps')


plt.figure()
plt.loglog(data[1,:],errR[0,:]*100,'*--',label="Law-Wendroff")
plt.loglog(data[1,:],errR[2,:]*100,'*--',label="ADER4")
plt.loglog(data[1,:],errRIm[0,:]*100,'*-',label="Law-Wendroff+II")
plt.loglog(data[1,:],errRIm[2,:]*100,'*-',label="ADER4+II")
plt.xlabel("$\Delta x$ (m)")
plt.ylabel("$\Delta\epsilon_v$ (%)")
plt.title("$\Delta\epsilon=\sqrt{\\frac{\sum(U-U_{th})^2}{\sum U_{th}^2}}$")
# plt.savefig('Figures/ErrRV_Npt3_ESIM3.eps', format='eps')
plt.legend()
plt.figure()
plt.loglog(data[1,:],errR[1,:]*100,'*--',label="Law-Wendroff")
plt.loglog(data[1,:],errR[3,:]*100,'*--',label="ADER4")
plt.loglog(data[1,:],errRIm[1,:]*100,'*-',label="Law-Wendroff+II")
plt.loglog(data[1,:],errRIm[3,:]*100,'*-',label="ADER4+II")
plt.xlabel("$\Delta x$ (m)")
plt.ylabel("$\Delta\epsilon_\sigma$ (%)")
plt.title("$\Delta\epsilon=\sqrt{\\frac{\sum(U-U_{th})^2}{\sum U_{th}^2}}$")
plt.legend()
# plt.savefig('Figures/ErrRS_Npt3_ESIM3.eps', format='eps')


pentLW=(np.log(errIm[0,6])-np.log(errIm[0,2]))/(np.log(data[1,6])-np.log(data[1,2]))
pentA4=(np.log(errIm[2,6])-np.log(errIm[2,2]))/(np.log(data[1,6])-np.log(data[1,2]))

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
