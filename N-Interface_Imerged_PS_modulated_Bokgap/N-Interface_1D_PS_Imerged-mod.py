#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:50:43 2023

@author: michael
"""

import os
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/N-Interface_Imerged_PS_modulated_Bokgap/")

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../WAPOD/")

import scheme1D_TDold as scheme1D
import scheme1D_TD_Homo as scheme1DH
import source
import inputReading
import figures as ffg
import generatePeriodic
import matplotlib as mpl

nbcell, nb_milieux, blocs =generatePeriodic.lire_fichier_milieu("MilieuPeriodique.txt")
nbcell, nb_frontieres, blocsF =generatePeriodic.lire_fichier_frontieres("FrontieresPeriodique.txt")

generatePeriodic.ecrire_milieu("Milieu.txt", nbcell, nb_milieux, blocs)
generatePeriodic.ecrire_frontieres_sortie("Frontiere.txt", nbcell, nb_frontieres, blocsF)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
################# Input files 
#%
configuration=inputReading.readfile('Demarrer.txt')
mat=inputReading.readfile('Milieu.txt')
frontiere=inputReading.readfile('Frontiere.txt')     
sourcef=inputReading.readfile('Source.txt')   
#sourceRev=inputReading.readfile('SourceRev.txt')
modulation=inputReading.readfile('Modulation.txt')
#nomodulation=inputReading.readfile('NoModulation.txt')
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
U0=scheme1D.initPointSourceProblem(configuration,frontiere,mat,sourcef)

###### Mise en oeuvre


Tmaxesp=0.25

Nt=int(Tmaxesp/dt)
Tmax=Nt*dt

t=np.linspace(0,Tmax,Nt+1)

import time,homogenize,math
time0=time.time()
fCFLh=homogenize.updateCFLP(Nt,dt,t,xmin,nb_frontieres,frontiere,mat)
CFLH=fCFLh*CFL
dt2v=np.min(scheme1D.timeStep(CFLH,dx,cm))
kt=math.ceil(dt/dt2v)
dt2p=min(dt2v,dt)

Nt=int(np.ceil(0.25/dt2p))
Tmax=Nt*dt
Nt2=int(Tmax/dt2p)
Ntp=max(Nt2,Nt)
tmax=str(Tmax)+"s"
K=4



tmax=str(Tmax)+"s"
tp=np.linspace(0,Tmax,Ntp+1)

sce=np.zeros(Ntp+1)
for i in range(Ntp+1):
    sce[i]=source.choice_timefct(sourcef,tp[i])
    
# sceR=np.zeros(Ntp+1)
# for i in range(Ntp+1):
#     sceR[i]=source.choice_timefct(sourceRev,tp[i])
    
    
from matplotlib.animation import FuncAnimation

pFilm=1
h=alpha[1+nb_frontieres]-alpha[1]
# h=alpha[1]-alpha[0]
# alpha0=alpha[0]

Clist=[]
DClist=[]
Mlist=[]
DMlist=[]
fmlist=[]
taulist=[]
for i in range(nb_frontieres-1):
    print(i)
    if float(frontiere["Kn_"+str(i)])!=0:
        Clist.append(1/float(frontiere["Kn_"+str(i)]))
    else:
        Clist.append(0)
    DClist.append(float(frontiere["DeltaF_"+str(i)]))
    Mlist.append(float(frontiere["Mn_"+str(i)]))
    DMlist.append(float(frontiere["DeltaM_"+str(i)]))
    fmlist.append(float(frontiere["Freq_"+str(i)]))
    taulist.append(1/float(frontiere["Freq_"+str(i)]))
# cel=float(modulation["ModC"])
# if cel !=0:
#     Nip=int(cel/h/fm)
# else:
#     Nip=1

rhoef=rho[0]+np.sum(Mlist)/h
unsE=1/(rho[0]*cm[0]*cm[0])+np.sum(Clist)/h
Eef=1/unsE

# K1=K0/(1+DK)
# K2=K0/(1-DK)

# M1=M0*(1+DM)
# M2=M0*(1-DM)

# cmet=cm[0]
# rhoet=rho[0]
# Eet=rhoet*cmet**2

# K0et=K0*h/Eet
# M0et=M0/h/rhoet
# K1et=K1*h/Eet
# M1et=M1/h/rhoet
# K2et=K2*h/Eet
# M2et=M2/h/rhoet

# if Nip%2==1:
#     rhoef1=rhoet+M0/h+DM*M0/h/Nip
#     Eef1=1/(1/Eet+1/K0/h+DK/K0/h/Nip)
#     rhoef2=rhoet+M0/h-DM*M0/h/Nip
#     Eef2=1/(1/Eet+1/K0/h-DK/K0/h/Nip)
# else:
#     rhoef1=rhoet+M0/h
#     Eef1=1/(1/Eet+1/K0/h)
#     rhoef2=rhoet+M0/h
#     Eef2=1/(1/Eet+1/K0/h)
# cmef1=np.sqrt(Eef1/rhoef1)
# cmef2=np.sqrt(Eef2/rhoef2)
# cmefmin=min(cmef1,cmef2)

# ratio=float(modulation["Ratio"])
# Eef=1/(ratio/Eef1+(1-ratio)/Eef2)
# rhoef=ratio*rhoef1+(1-ratio)*rhoef2


# c_r = cmef1*ratio+cmef2*(1-ratio)
c_0 = np.sqrt(Eef/rhoef)


eta0=2*np.pi*fs/c_0*h
eta1=2*np.pi*(fs+max(fmlist))/c_0*h




############### MODULATION DIRECT PROPAGATION ##########################################
import time
time0=time.time()
U,dt2=scheme1DH.temporalSchemePS0MKP(Ntp,dt2p,configuration,sourcef,nb_frontieres,frontiere,mat)
Res,Vfilm,Sfilm,CSVM=scheme1D.FD_sourceVII(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm,rCSVM="yes")#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
timef=time.time()
maxV=np.max([np.max(Vfilm),-np.min(Vfilm)])
timeN_mod=timef-time0
print(timeN_mod)  


# nVfilm='Vfilm'+str(fs)+'fm'+str(fmlist[0])+'.txt'
# np.savetxt(nVfilm, Vfilm,fmt='%18e', delimiter='\t')
# nSveCS='CSVM'+str(fs)+'fm'+str(fmlist[0])+'.txt'
# np.savetxt(nSveCS, CSVM,fmt='%18e', delimiter='\t')
# nU='Uhomo'+str(fs)+'fm'+str(fmlist[0])+'.txt'
# np.savetxt(nU, U,fmt='%18e', delimiter='\t')
# nSfilm='Sfilm'+str(fs)+'fm'+str(fmlist[0])+'.txt'
# np.savetxt(nSfilm, Sfilm,fmt='%18e', delimiter='\t')
# nRes='Rfin'+str(fs)+'fm'+str(fmlist[0])+'.txt'
# np.savetxt(nRes, Res,fmt='%18e', delimiter='\t')
# ntp='tp'+str(fs)+'fm'+str(fmlist[0])+'.txt'
# np.savetxt(ntp, tp,fmt='%18e', delimiter='\t')

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

# energyM=np.zeros(int(Nt/pFilm))
# for ti in range(int(Nt/pFilm)):
#     energyM[ti]=1/2*np.sum(rhox*Vfilm[ti]*Vfilm[ti])+1/2*np.sum(1/rhox/Celx/Celx*Sfilm[ti]*Sfilm[ti])


# Nt2=int(Tmax/dt)
K=4
# timef=time.time()
# ResH,VfilmH,SfilmH=scheme1D.FD_sourceVIIHomo1(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
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
# plt.plot(Vfilm[:,0:1000],'b')4
# plt.plot(VfilmH[:,0:1000],'k')


# ################ GODUNOV SPLITTING ####################################
# plt.grid()
# total_frames = Vfilm.shape[0]
# fig, ax = plt.subplots()
# timef=time.time()
# ResGS,VfilmGS,SfilmGS=scheme1D.FD_sourceVIIHomoGS(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)
# timeN_mod=timef-time0
# print(timeN_mod)  
# ################ STRANG SPLITTING ####################################
# plt.grid()
# total_frames = Vfilm.shape[0]
# fig, ax = plt.subplots()
timef=time.time()
# ResSS,VfilmSS,SfilmSS=scheme1D.FD_sourceVIIHomoSS2(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)
timeN_mod=timef-time0
print(timeN_mod)  
# # # ###################### NO MODULATION #################################
# time0=time.time()
# ResNM,VfilmNM,SfilmNM=scheme1D.FD_sourceVII(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,nomodulation,pFilm)
# timeNM=time.time()-time0
# print(timeNM)  


import modulation as modul

Ufilm=ffg.VtoU(Vfilm,dt2p)
V=ffg.UtoV(U,dt2p)
S=ffg.UtoS(U,dx)
maxU=np.max(Ufilm)

Ntfig=int(np.floor(0.2/dt2p))
maxVt=np.max([np.max(Vfilm[Ntfig,:]),-np.min(Vfilm[Ntfig,:])])

psmax=np.where(Ufilm[Ntfig-1,:Nx//2]==max(Ufilm[Ntfig-1,:Nx//2]))
xzoommin=X[psmax[0][0]-200]
xzoommax=X[psmax[0][0]+600]
zoomlim=[xzoommin,xzoommax]

vals = []
for i in range(alpha.size-1):
    if frontiere["Contact_"+str(i)] != "parfait":
        v = modul.fonctionModMP(Ntfig*dt2p, alpha[i], frontiere, i)
        vals.append(v)
    else:
        vals.append(0)
vals = np.array(vals)

cmap = plt.cm.PiYG   # à ajuster si tu veux
norm = mpl.colors.Normalize(vmin=2*vals.min(), vmax=2*vals.max())



# # Cree l'animation
# # animation = FuncAnimation(fig, update, frames=total_frames, interval=1)
# # Affiche l'animation
# plt.show()
def updateComp(frame):
    plt.clf()
    # for i in range(alpha.size-1):
    #     if modul.fonctionModMP(pFilm*dt2p*frame,alpha[i],frontiere,i)<0.:
    #         plt.axvline(x=alpha[i],color='palegreen',alpha=-modul.fonctionModMP(pFilm*dt2p*frame,alpha[i],frontiere,i))#,linestyle='--')
    #     elif modul.fonctionModMP(pFilm*dt2p*frame,alpha[i],frontiere,i)>0.:
    #         plt.axvline(x=alpha[i],color='plum',alpha=modul.fonctionModMP(pFilm*dt2p*frame,alpha[i],frontiere,i))#,linestyle='--')
    #     else:
    #         plt.axvline(x=alpha[i],color='lavender')
    vals = []
    for i in range(alpha.size-1):
        if frontiere["Contact_"+str(i)] != "parfait":
            v = modul.fonctionModMP(frame*dt2p, alpha[i], frontiere, i)
            vals.append(v)
        else:
            vals.append(0)
    vals = np.array(vals)
    for i in range(alpha.size):
        if frontiere["Contact_"+str(i)]!="parfait":
            plt.axvline(x=alpha[i],ymin=0,ymax=1,color=cmap(vals[i]),alpha=np.abs(vals[i]/2),linewidth=0.5)#,linestyle='--')
            # elif modul.fonctionModMP(Ntfig*dt2p,alpha[i],frontiere,i)>0.:
            #     plt.axvline(x=alpha[i],ymin=0,ymax=1,color='orchid',alpha=modul.fonctionModMP(Ntfig*dt2p,alpha[i]/2,frontiere,i))
            #plt.colorbar()#,linestyle='--')
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.scatter(X[::1],Ufilm[frame,::1],marker=r'$\bigcirc$',label=r'$U_h$',color='red',s=30)
    # plt.plot(X[::1],Ufilm[Ntfig,::1],linestyle=None,marker='o',color='white',markersize=2)
    plt.plot(X,U[frame,:],linewidth=3,label=r'$U_0$')
    # plt.plot(X, Ufilm[frame,:],'*')#, label=f'Frame {frame}')
    # plt.plot(X, U[frame,:],'-')#, label=f'Frame {frame}')
    # plt.plot(X, VfilmGS[frame,:],"-")#, label=f'Frame {frame}')
    # plt.plot(X, VfilmSS[frame,:],"-")#, label=f'Frame {frame}')
    # plt.plot(X, VfilmH[frame,:],"--")#, label=f'Frame {frame}')
    # plt.plot(X, VfilmNM[frame,:],"o")#, label=f'Frame {frame}')
    # plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0000020,0.000003], aspect='auto', alpha=0.35)
    # plt.colorbar()
    # plt.plot(X, VfilmH[frame,:])
    plt.xlabel(r'$X$ (m)',fontsize=16)
    plt.ylabel(r'$U$ (m)',fontsize=16)
    plt.ylim([-0.0000050,0.000005])
    plt.tight_layout()
    # plt.legend()


total_frames = Vfilm.shape[0]
fig, ax = plt.subplots()

# Cree l'animation
animationComp = FuncAnimation(fig, updateComp, frames=total_frames, interval=1)
plt.show()
# nomAnim="Film_MicroH_"+sourcef["Frequence"]+frontiere["Freq_0"]+"Hz.mp4"
# animationComp.save(nomAnim, writer='ffmpeg', fps=20)


# def updateComp(frame):
#     plt.clf()
#     for i in range(alpha.size-1):
#         if modul.fonctionModM(pFilm*dt2p*frame,alpha[i],modulation,frontiere,i)<0.:
#             plt.axvline(x=alpha[i],color='red')#,linestyle='--')
#         else:
#             plt.axvline(x=alpha[i],color='lime')#,linestyle='--')
#     plt.plot(X, Ufilm[frame,:],'*')#, label=f'Frame {frame}')
#     # plt.plot(X, UfilmGS[frame,:],"-")#, label=f'Frame {frame}')
#     plt.plot(X, VfilmSS[frame,:],"-")#, label=f'Frame {frame}')
#     # plt.plot(X, VfilmH[frame,:],"--")#, label=f'Frame {frame}')
#     # plt.plot(X, VfilmNM[frame,:],"o")#, label=f'Frame {frame}')
#     plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -maxU,maxU], aspect='auto', alpha=0.35)
#     plt.colorbar()
#     # plt.plot(X, VfilmH[frame,:])
#     plt.xlabel("$x$ (m)",fontsize=20)
#     plt.ylabel("$v$ (m/s)",fontsize=20)
#     plt.ylim([-maxU,maxU])
#     plt.tight_layout()
    # plt.legend()


# total_frames = VfilmSS.shape[0]
# fig, ax = plt.subplots()

# # Cree l'animation
# animationComp = FuncAnimation(fig, updateComp, frames=total_frames, interval=10)
# plt.show()



if float(frontiere["Kn_0"])==0:
    K0=np.infty
else:
    K0=float(frontiere["Kn_0"])
frK=float(modulation["Freq"])
delta=float(modulation["DeltaF"])

M0=float(frontiere["Mn_0"])
frM=float(modulation["Freq"])
deltaM=float(modulation["DeltaM"])

nfront=frontiere["nombre de frontieres "]
import modulation as modulateP
import homogenize as homoP
# for i in range(nfront):
#     modula=modulateP.fonctionModK(t,alpha[i],modulation,frontiere,i)
#     modulb=modulateP.fonctionModM(t,alpha[i],modulation,frontiere,i)
#     print(nfront)

import sympy as sp
rhoefft,Eefft,Cefft,QCefft,QMefft=homoP.propHomoFormP(xmin,mat,frontiere,nb_frontieres)
t=sp.symbols('t')
Ceff=sp.lambdify(t,Cefft)
rhoeff=sp.lambdify(t,rhoefft)
modulb=modulateP.fonctionModMP(tp,alpha[0],frontiere,0)
modulM=M0*(1+modulb)
# dmodulb=np.zeros(int(Nt/pFilm))
# for i in range(int(Nt/pFilm)):
#     modulb[i]=M0*(1+deltaM*np.sin(2*np.pi*frM*t[pFilm*i]))
# #     dmodulb[i]=2*np.pi*frM*M0*(deltaM*np.cos(2*np.pi*frM*t[pFilm*i]))


modula=modulateP.fonctionModKP(tp,alpha[0],frontiere,0)
modulK=1/K0*(1+modula)
# dmodula=np.zeros(int(Ntfig/pFilm))

energyV=np.zeros(int(Nt2/pFilm))
energyI=np.zeros(int(Nt2/pFilm))
energyhomo=np.zeros(int(Nt2/pFilm))

for ti in range(int(Nt2/pFilm)):
    energyV[ti]=dx/2*np.sum(rhox*Vfilm[ti]*Vfilm[ti])+dx/2*np.sum(1/rhox/Celx/Celx*Sfilm[ti]*Sfilm[ti])
    energyhomo[ti]=dx/2*np.sum(rhoeff(ti*dt2p)*V[ti]*V[ti])+dx/2*np.sum(rhoeff(ti*dt2p)*Ceff(ti*dt2p)*Ceff(ti*dt2p)*S[ti]*S[ti])
    for interi in range(int(nfront)):
        if frontiere["Contact_"+str(interi)]!="parfait":
            energyI[ti]=energyI[ti]+1/2*modulM[ti]*CSVM[ti,2+4*interi]*CSVM[ti,2+4*interi]+1/2*modulK[ti]*CSVM[ti,3+4*interi]*CSVM[ti,3+4*interi]     
energymicro=energyV+energyI

plt.figure()
plt.semilogy(tp[:-1],energymicro,label=r'$\mathcal{E}_{h}$',linewidth=3)
plt.semilogy(tp[:-1],energyhomo,label=r'$\mathcal{E}_{0}$',linewidth=3)
plt.xlabel(r'$T$ (s)',fontsize=16)
plt.ylabel(r'$\mathcal{E}$ (J)',fontsize=16)
plt.xlim([0,tp[-1]])
plt.ylim([2e-3,max(energyhomo)])
plt.imshow([1+modula], cmap='gray', extent=[0, tp[-1], 0,max(energyhomo)], aspect='auto', alpha=0.35)
# plt.title("$t= $"+str(round(tp[Nt2], 3))+" s",fontsize=18)
plt.legend(fontsize=16)#loc="lower right",
plt.colorbar().set_label(label=r'$1+\varepsilon_{C,1}~\sin (\Omega_m T)$',size=16)
plt.tight_layout()

ffg.sismo(X, tp, U,12,pFilm,Nt=2400)

Ntfig=int(np.floor(0.2/dt2p))
maxVt=np.max([np.max(Vfilm[Ntfig,:]),-np.min(Vfilm[Ntfig,:])])

plt.figure()
psmax=np.where(Ufilm[Ntfig-1,:Nx//2]==max(Ufilm[Ntfig-1,:Nx//2]))
xzoommin=X[psmax[0][0]-200]
xzoommax=X[psmax[0][0]+600]
zoomlim=[xzoommin,xzoommax]
# for i in range(alpha.size):
#     plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(X[::1],Ufilm[Ntfig,::1],marker=r'$\bigcirc$',label=r'$U_h$',color='blue',s=40)
#plt.plot(X[::1],Ufilm[Ntfig,::1],linestyle='None',marker='o',color='white',markersize=3)
#plt.plot(X,Uz[Ntfig-1,:],'grey',linewidth=3,label=r'$U_0$')
plt.plot(X,U[Ntfig,:],color='red',linestyle='None',marker='o',markersize=3,label=r'$U_0$')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.title("$T=$"+str(round(tp[Ntfig], 3))+" s",fontsize=16)
plt.legend(fontsize=14)#loc="lower right",)
plt.tight_layout()


psmax=np.where(Ufilm[Ntfig-1,:Nx//2]==max(Ufilm[Ntfig-1,:Nx//2]))
xzoommin=X[psmax[0][0]-200]
xzoommax=X[psmax[0][0]+600]
zoomlim=[xzoommin,xzoommax]

plt.figure()
plt.axvspan(xzoommin, xzoommax, facecolor='g', alpha=0.3)
# for i in range(alpha.size):
#     plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(X[::1],Ufilm[Ntfig,::1],marker=r'$\bigcirc$',label=r'$U_h$',color='red',s=30)
# plt.plot(X[::1],Ufilm[Ntfig,::1],linestyle=None,marker='o',color='white',markersize=2)
plt.plot(X,U[Ntfig,:],linewidth=3,label=r'$U_0$')
# plt.axvline(x=alpha[0],color='black')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.title("$t= $"+str(round(tp[Ntfig], 3))+" s",fontsize=18)
plt.legend(fontsize=16)#loc="lower right",
plt.tight_layout()
namefig="Figures/PRSA_U_MK_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".pdf"
plt.savefig(str(namefig), format='pdf')
namefig="Figures/PRSA_U_MK_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".png"
plt.savefig(str(namefig), format='png')

plt.figure()
plt.axvspan(xzoommin, xzoommax, facecolor='g', alpha=0.3)
# for i in range(alpha.size):
#     plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(X[::1],Vfilm[Ntfig,::1],marker=r'$\bigcirc$',label=r'$V_h$',color='red',s=30)
# plt.plot(X[::1],Vfilm[Ntfig,::1],linestyle=None,marker='o',color='white',markersize=2)
plt.plot(X,V[Ntfig,:],linewidth=3,label=r'$V_0$')
# plt.axvline(x=alpha[0],color='black')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$V$ (m/s)',fontsize=16)
plt.title("$t= $"+str(round(tp[Ntfig], 3))+" s",fontsize=18)
plt.legend(fontsize=16)#loc="lower right",
plt.tight_layout()
# plt.plot(VfilmSS[500,:])
# plt.plot(VfilmH[500,:])
# plt.title("$nt=$500")
namefig="Figures/PRSA_V_MK_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".pdf"
plt.savefig(str(namefig), format='pdf')
namefig="Figures/PRSA_V_MK_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".png"
plt.savefig(str(namefig), format='png')


plt.figure()
for i in range(alpha.size):
    plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(X[::1],Ufilm[Ntfig,::1],marker=r'$\bigcirc$',label=r'$U_h$',color='red',s=30)
# plt.plot(X[::1],Ufilm[Ntfig,::1],linestyle=None,marker='o',color='white',markersize=2)
plt.plot(X,U[Ntfig,:],linewidth=3,label=r'$U_0$')
# plt.axvline(x=alpha[0],color='black')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.title("$t= $"+str(round(tp[Ntfig], 3))+" s",fontsize=18)
plt.xlim(zoomlim)
plt.legend(fontsize=16)#loc="lower right",
plt.tight_layout()
namefig="Figures/PRSAzoom_U_MK_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".pdf"
plt.savefig(str(namefig), format='pdf')
namefig="Figures/PRSAzoom_U_MK_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".png"
plt.savefig(str(namefig), format='png')

plt.figure()
for i in range(alpha.size):
    plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(X[::1],Vfilm[Ntfig,::1],marker=r'$\bigcirc$',label=r'$V_h$',color='red',s=30)
plt.plot(X[::1],Vfilm[Ntfig,::1],linestyle=None,marker='o',color='white',markersize=2)
plt.plot(X,V[Ntfig,:],linewidth=3,label=r'$V_0$')
# plt.axvline(x=alpha[0],color='black')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$V$ (m/s)',fontsize=16)
plt.title("$t= $"+str(round(tp[Ntfig], 3))+" s",fontsize=18)
plt.xlim(zoomlim)
plt.legend(fontsize=16)#loc="lower right",
plt.tight_layout()
# plt.plot(VfilmSS[500,:])
# plt.plot(VfilmH[500,:])
# plt.title("$nt=$500")
namefig="Figures/PRSAzoom_V_MK_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".pdf"
plt.savefig(str(namefig), format='pdf')
namefig="Figures/PRSAzoom_V_MK_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".png"
plt.savefig(str(namefig), format='png')

plt.figure()
plt.plot(Vfilm[1000,:])
# plt.plot(VfilmGS[1000,:])
plt.plot(V[1000,:])
# plt.plot(VfilmH[1000,:])
plt.xlabel(r"$x$ (m)")
plt.ylabel(r"$v$ (m/s)")
plt.title("$nt=$500")
# namefig="Figures/WP1000_sin1_MK_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".eps"
# plt.savefig(str(namefig), format='eps')

def cosine_taper(signal, alpha=0.05,zp=16):
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


import funcPWE
fctModE="np.sin"
fctModR="np.sin"
NE=8
tau=1/fmlist[0]
Emicro=rho[0]*cm[0]*cm[0]
EpsR=Mlist[0]*DMlist[0]/(Mlist[0]+rho[0]*h)
EpsE=Clist[0]*Emicro*DClist[0]/(Clist[0]*Emicro+h)
if tau!=0:
    omegavect,VP=funcPWE.fillMatrix(tau,Eef,rhoef,EpsE,EpsR,fctModE,fctModR,NE,N=16)
    kmax=max(np.real(VP[NE-1,:])/2/np.pi)
    ffg.dispFKB(Ufilm[:Ntfig,:],dt2p,dx)
    for i in range(NE):
        plt.plot((VP[i,:])/2/np.pi,omegavect/2/np.pi,'--r',linewidth=3)
#     ffg.dispFK_tperiodic2(U[:Ntfig,:],dt2p,dx,fm,kmax)
#     for i in range(NE):
#         plt.plot((VP[i,:])/2/np.pi,omegavect/2/np.pi,'--r',linewidth=3)
# # N=2**15
# # FFTMaptot={}
# # for i in range(nsou):
# #     for j in range(nmod):
# #         sigcor = cosine_taper(results[f'Vfilm{i}{j}'][::,int(750)])
# # # plt.figure()
# # # plt.plot(Vfilm[::,int(750)])
# # # plt.plot(sigcor)
# # fftSig=np.fft.fftshift(np.fft.fft(sigcor,n=N))
# # freq = np.fft.fftshift(np.fft.fftfreq(N,d=dt))
# # plt.figure()
# # plt.plot(np.abs(freq),np.abs(fftSig))


