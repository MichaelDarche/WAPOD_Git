#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:50:43 2023

@author: michael
"""

import os
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/N-Interface_Imerged_PS_modulated_homophi_ImpMatchP/")

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
nbcell, nb_frontieres, blocsF =generatePeriodic.lire_fichier_frontieres("FrontieresPeriodique2.txt")

generatePeriodic.ecrire_milieu("Milieu2.txt", nbcell, nb_milieux, blocs)
generatePeriodic.ecrire_frontieres_sortie("Frontiere2.txt", nbcell, nb_frontieres, blocsF)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
################# Input files 
#%
configuration=inputReading.readfile('Demarrer.txt')
mat=inputReading.readfile('Milieu2.txt')
frontiere=inputReading.readfile('Frontiere2.txt')     
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


Tmaxesp=0.3

Nt=int(Tmaxesp/dt)
Tmax=Nt*dt

t=np.linspace(0,Tmax,Nt+1)

import time,homogenize,math
time0=time.time()
fCFLh=homogenize.updateCFL(Nt,dt,t,rho,cm,frontiere,modulation)
CFLH=fCFLh*CFL
dt2v=np.min(scheme1D.timeStep(CFLH,dx,cm))
kt=math.ceil(dt/dt2v)
dt2p=min(dt2v,dt)

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



c_0 = np.sqrt(Eef/rhoef)


eta0=2*np.pi*fs/c_0*h
eta1=2*np.pi*(fs+max(fmlist))/c_0*h

# rhoef0=rhoet+M0/h
# EpsR=DM*M0/h/Nip/rhoef0
# Eef0=Eet*K0*h/(K0*h+Eet)
# EpsE=DK*Eef0/K0/h/Nip

timef=time.time()
#ResSS,VfilmSS,SfilmSS=scheme1D.FD_sourceVIIHomoSS2(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)
# U,dt2=scheme1DH.temporalSchemePS0MK(Ntp,dt2p,configuration,sourcef,modulation,frontiere,mat)
timeN_mod=timef-time0
print(timeN_mod)  

############### MODULATION DIRECT PROPAGATION #############################
import time
time0=time.time()
U,dt2=scheme1DH.temporalSchemePS0MKP(Ntp,dt2p,configuration,sourcef,nb_frontieres,frontiere,mat)
Res,Vfilm,Sfilm,CSVM=scheme1D.FD_sourceVII(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm,rCSVM="yes")#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
timef=time.time()
maxV=np.max([np.max(Vfilm),-np.min(Vfilm)])
timeN_mod=timef-time0
print(timeN_mod)  
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

energyM=np.zeros(int(Nt/pFilm))
for ti in range(int(Nt/pFilm)):
    energyM[ti]=1/2*np.sum(rhox*Vfilm[ti]*Vfilm[ti])+1/2*np.sum(1/rhox/Celx/Celx*Sfilm[ti]*Sfilm[ti])


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
#ResSS,VfilmSS,SfilmSS=scheme1D.FD_sourceVIIHomoSS2(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)
# U,dt=scheme1DH.temporalScheePS0MK(Nt,dt,configuration,source,modulation,frontiere,mat)
timeN_mod=timef-time0
print(timeN_mod)  
# # # ###################### NO MODULATION #################################
# time0=time.time()
# ResNM,VfilmNM,SfilmNM=scheme1D.FD_sourceVII(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,nomodulation,pFilm)
# timeNM=time.time()-time0
# print(timeNM)  
# nSveVmicro='V_MK-IM_fc'+str(fs)+'fm'+str(fm)+'.txt'
# nSveVh='Vh_MK-IM_fc'+str(fs)+'fm'+str(fm)+'.txt'
# np.savetxt(nSveVmicro, Vfilm, fmt='%18e', delimiter='\t')
# np.savetxt(nSveVh, VfilmSS, fmt='%18e', delimiter='\t')

# nSveSmicro='S_MK-IM_fc'+str(fs)+'fm'+str(fm)+'.txt'
# nSveSh='Sh_MK-IM_fc'+str(fs)+'fm'+str(fm)+'.txt'
# np.savetxt(nSveSmicro, Sfilm, fmt='%18e', delimiter='\t')
# np.savetxt(nSveSh, SfilmSS, fmt='%18e', delimiter='\t')

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


nVfilm='Vfilm'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nVfilm, Vfilm,fmt='%18e', delimiter='\t')
nU='Uhomo'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nU, U,fmt='%18e', delimiter='\t')
nSfilm='Sfilm'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nSfilm, Sfilm,fmt='%18e', delimiter='\t')
nRes='Rfin'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nRes, Res,fmt='%18e', delimiter='\t')
ntp='tp'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(ntp, tp,fmt='%18e', delimiter='\t')


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
    plt.ylim([-0.0000020,0.000008])
    plt.tight_layout()
    # plt.legend()


total_frames = U.shape[0]
fig, ax = plt.subplots()

# Cree l'animation
animationComp = FuncAnimation(fig, updateComp, frames=total_frames, interval=1)
plt.show()
nomAnim="Film_1I_"+sourcef["Frequence"]+frontiere["Freq_0"]+"Hz.mp4"
animationComp.save(nomAnim, writer='ffmpeg', fps=20)



# plt.figure()
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


#ffg.sismo(X, t, VfilmSS,20,pFilm,Nt=2000)

Ntfig=int(np.floor(0.2/dt2p))
maxVt=np.max([np.max(Vfilm[Ntfig,:]),-np.min(Vfilm[Ntfig,:])])

psmax=np.where(Ufilm[Ntfig-1,:Nx//2]==max(Ufilm[Ntfig-1,:Nx//2]))
xzoommin=X[psmax[0][0]-200]
xzoommax=X[psmax[0][0]+600]
zoomlim=[xzoommin,xzoommax]


plt.figure()
plt.axvspan(xzoommin, xzoommax, facecolor='g', alpha=0.3)
# for i in range(alpha.size):
#     plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(X[::1],Ufilm[Ntfig,::1],marker=r'$\bigcirc$',label=r'$U_h$',color='blue',s=40)
#plt.plot(X[::1],Ufilm[Ntfig,::1],linestyle='None',marker='o',color='white',markersize=3)
#plt.plot(X,Uz[Ntfig-1,:],'grey',linewidth=3,label=r'$U_0$')
plt.plot(X,U[Ntfig,:],color='red',linestyle='None',marker='o',markersize=3,label=r'$U_0$')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.title("$T=$"+str(round(tp[Ntfig], 3))+" s",fontsize=16)
plt.legend(fontsize=14)#loc="lower right",)<
plt.tight_layout()


plt.figure()
plt.axvspan(xzoommin, xzoommax, facecolor='g', alpha=0.3)
# for i in range(alpha.size):
#     plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(X[::1],Ufilm[Ntfig,::1],marker=r'$\bigcirc$',label=r'$U_h$',color='red',s=30)
# plt.plot(X[::1],Ufilm[Ntfig,::1],linestyle=None,marker='o',color='white',markersize=2)
# plt.plot(X,UfilmSS[Ntfig,:],linewidth=3,label=r'$U_0$')
plt.plot(X,U[Ntfig,:],linewidth=3,label=r'$U_0$')
# plt.axvline(x=alpha[0],color='black')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.title("$t= $"+str(round(tp[Ntfig], 3))+" s",fontsize=18)
plt.legend(fontsize=16)#loc="lower right",
plt.tight_layout()
# namefig="Figures/PRSA_IIM_U_MK_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".pdf"
# plt.savefig(str(namefig), format='pdf')
# namefig="Figures/PRSA_IIM_U_MK_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".png"
# plt.savefig(str(namefig), format='png')

# plt.figure()
# plt.axvspan(xzoommin, xzoommax, facecolor='g', alpha=0.3)
# # for i in range(alpha.size):
# #     plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
# plt.scatter(X[::1],Vfilm[Ntfig,::1],marker=r'$\bigcirc$',label=r'$V_h$',color='red',s=30)
# # plt.plot(X[::1],Vfilm[Ntfig,::1],linestyle=None,marker='o',color='white',markersize=2)
# plt.plot(X,VfilmSS[Ntfig,:],linewidth=3,label=r'$V_0$')
# # plt.axvline(x=alpha[0],color='black')
# plt.xlabel(r'$X$ (m)',fontsize=16)
# plt.ylabel(r'$V$ (m/s)',fontsize=16)
# plt.title("$t= $"+str(round(tp[Ntfig], 3))+" s",fontsize=18)
# plt.legend(fontsize=16)#loc="lower right",
# plt.tight_layout()
# # plt.plot(VfilmSS[500,:])
# # plt.plot(VfilmH[500,:])
# # plt.title("$nt=$500")
# namefig="Figures/PRSA_IIM_V_MK_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".pdf"
# plt.savefig(str(namefig), format='pdf')
# namefig="Figures/PRSA_IIM_V_MK_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".png"
# plt.savefig(str(namefig), format='png')


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
# namefig="Figures/PRSAzoom_IIM_U_MK_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".pdf"
# plt.savefig(str(namefig), format='pdf')
# namefig="Figures/PRSAzoom_IIM_U_MK_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".png"
# plt.savefig(str(namefig), format='png')

# plt.figure()
# for i in range(alpha.size):
#     plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
# plt.scatter(X[::1],Vfilm[Ntfig,::1],marker=r'$\bigcirc$',label=r'$V_h$',color='red',s=30)
# # plt.plot(X[::1],Vfilm[Ntfig,::1],linestyle=None,marker='o',color='white',markersize=2)
# plt.plot(X,VfilmSS[Ntfig,:],linewidth=3,label=r'$V_0$')
# # plt.axvline(x=alpha[0],color='black')
# plt.xlabel(r'$X$ (m)',fontsize=16)
# plt.ylabel(r'$V$ (m/s)',fontsize=16)
# plt.xlim(zoomlim)
# plt.title("$t= $"+str(round(tp[Ntfig], 3))+" s",fontsize=18)
# plt.legend(fontsize=16)#loc="lower right",
# plt.tight_layout()
# # plt.plot(VfilmSS[500,:])
# # plt.plot(VfilmH[500,:])
# # plt.title("$nt=$500")
# namefig="Figures/PRSAzoom_IIM_V_MK_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".pdf"
# plt.savefig(str(namefig), format='pdf')
# namefig="Figures/PRSAzoom_IIM_V_MK_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".png"
# plt.savefig(str(namefig), format='png')

# plt.figure()
# plt.plot(Vfilm[1000,:])
# # plt.plot(VfilmGS[1000,:])
# plt.plot(VfilmSS[1000,:])
# # plt.plot(VfilmH[1000,:])
# plt.xlabel(r"$x$ (m)")
# plt.ylabel(r"$v$ (m/s)")
# plt.title("$nt=$500")
# namefig="Figures/WP1000_sin1_MK_"+sourcef["Frequence"]+"_"+modulation["Freq"]+".eps"
# plt.savefig(str(namefig), format='eps')

import funcPWE
fctModE="np.sin"
fctModR="np.sin"
NE=8
if tau!=0:
    omegavect,VP=funcPWE.fillMatrix(tau,Eef0,rhoef0,EpsE,EpsR,fctModE,fctModR,NE,N=16)
    kmax=max(np.real(VP[NE-1,:])/2/np.pi)
    ffg.dispFK_tperiodic2(Ufilm[:Ntfig,:],dt2p,dx,fm,kmax)
    for i in range(NE):
        plt.plot((VP[i,:])/2/np.pi,omegavect/2/np.pi,'--r',linewidth=3)
    ffg.dispFK_tperiodic2(Vfilm[:Ntfig,:],dt2p,dx,fm,kmax)
    for i in range(NE):
        plt.plot((VP[i,:])/2/np.pi,omegavect/2/np.pi,'--r',linewidth=3)
# N=2**15
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



# # plt.figure()
# # plt.plot(Vfilm[2000,:])
# # plt.plot(VfilmShft[2000,:])
# #FFTMapSft=np.zeros([int(N),Nx])
# freq = np.fft.fftshift(np.fft.fftfreq(Nx,d=dx))
# xpt=635
# # sigsrc=cosine_taper(sce/rho[0]/2)
# # FFTS=np.abs(np.fft.fftshift(np.fft.fft(sigsrc,n=int(N))))
# sigcor = cosine_taper(Vfilm[int(xpt),::])
# FFT0=np.abs(np.fft.fftshift(np.fft.fft(sigcor,n=int(Nx))))

# xpt=230
# # sigsrc=cosine_taper(sce/rho[0]/2)
# # FFTS=np.abs(np.fft.fftshift(np.fft.fft(sigsrc,n=int(N))))
# sigcor = cosine_taper(Vfilm[int(xpt),::])
# FFT1=np.abs(np.fft.fftshift(np.fft.fft(sigcor,n=int(Nx))))

# xpt=1000
# # sigsrc=cosine_taper(sce/rho[0]/2)
# # FFTS=np.abs(np.fft.fftshift(np.fft.fft(sigsrc,n=int(N))))
# sigcor = cosine_taper(Vfilm[int(xpt),::])
# FFT2=np.abs(np.fft.fftshift(np.fft.fft(sigcor,n=int(Nx))))
# # sigcorH = cosine_taper(VfilmH[::,int(xpt)])
# # FFTH=np.abs(np.fft.fftshift(np.fft.fft(sigcorH,n=int(N))))
# # sigcorGS = cosine_taper(VfilmGS[::,int(xpt)])
# # FFTGS=np.abs(np.fft.fftshift(np.fft.fft(sigcorGS,n=int(N))))
# # sigcorNM = cosine_taper(VfilmNM[::,int(xpt)])    
# # FFTNM=np.abs(np.fft.fftshift(np.fft.fft(sigcorNM,n=int(N)))) 
# # sigcorSS = cosine_taper(VfilmSS[::,int(xpt)])
# # FFTSS=np.abs(np.fft.fftshift(np.fft.fft(sigcorSS,n=int(N))))

# plt.figure()
# plt.plot(freq,FFT0)
# plt.plot(freq,FFT1)
# plt.plot(freq,FFT2)
# plt.plot(freq,FFTSS)
# plt.plot(freq,FFTGS)
# plt.plot(freq,FFTH)
# plt.plot(freq,FFTNM)
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
