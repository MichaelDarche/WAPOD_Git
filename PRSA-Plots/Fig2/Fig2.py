#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 16:53:54 2026

@author: michael
"""

import sys
import os

sys.path.append(os.path.abspath(".."))


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


plt.rc('text', usetex=True)
plt.rc('font', family='serif')



import functions as func


##################
## Load parameters
##################

configuration=func.readfile('Demarrer.txt')
mat=func.readfile('Milieu.txt')
frontiere=func.readfile('Frontiere.txt')     
sourcef=func.readfile('Source.txt')   

X,xmin,xmax,Nx,dx,ESIM,Npt=func.geometric(configuration)
fs,CFL,x_0=func.source(sourcef)
Nmat,rho,cm=func.material(mat)
alpha=func.frontiere(frontiere)
rhox,Celx=func.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
nbFront=int(frontiere["nombre de frontieres "])

nbFrontieresPerCell=2
pFilm=1
h=alpha[1+nbFrontieresPerCell]-alpha[1]


Clist=[]
DClist=[]
Mlist=[]
DMlist=[]
fmlist=[]
taulist=[]
for i in range(nbFrontieresPerCell-1): # -1 because the last one is a perfect interface
    if float(frontiere["Kn_"+str(i)])!=0:
        Clist.append(1/float(frontiere["Kn_"+str(i)]))
    else:
        Clist.append(0)
    DClist.append(float(frontiere["DeltaF_"+str(i)]))
    Mlist.append(float(frontiere["Mn_"+str(i)]))
    DMlist.append(float(frontiere["DeltaM_"+str(i)]))
    fmlist.append(float(frontiere["Freq_"+str(i)]))
    taulist.append(1/float(frontiere["Freq_"+str(i)]))    
fm=fmlist[0]


##################
## Load datas
##################
# Raw datas
nSveResV='Vfilm'+str(fs)+'fm'+str(fm)+'.txt'
nSveResS='Sfilm'+str(fs)+'fm'+str(fm)+'.txt'
nSveCS='CS0'+str(fs)+'fm'+str(fm)+'.txt'
nSveHomo='Uhomo'+str(fs)+'fm'+str(fm)+'.txt'
nSvetp='tp'+str(fs)+'fm'+str(fm)+'.txt'
 
Vfilm=np.loadtxt(nSveResV, delimiter='\t')#Velocity (micro)
Sfilm=np.loadtxt(nSveResS, delimiter='\t')#Stress field (micro)
CSVM=np.loadtxt(nSveCS, delimiter='\t')#Jump conditions and Mean values at the interfaces
U=np.loadtxt(nSveHomo, delimiter='\t')#Displacement field (homogenised)
tp=np.loadtxt(nSvetp,delimiter='\t')#Time vector


# Postprocessing
dt=tp[1]-tp[0]#Time step
Ufilm=func.VtoU(Vfilm,dt)#Displacement field (micro)
V=func.UtoV(U,dt)#Velocity (homogenised)
S=func.UtoS(U,dx)#Space gradient of the displacement field (homogenised)

##################
## Modulated parameters
##################
# Interface parameters
modula=func.fonctionModKP(tp,frontiere,i)
modulC=Clist[0]*(1+modula)#Time dependance of the compliance C
modulb=func.fonctionModMP(tp,frontiere,i)
modulM=Mlist[0]*(1+modulb)#Time dependance of the inertia M

#Homogenized parameters
rhoef=rho[0]+np.sum(Mlist)/h
unsE=1/(rho[0]*cm[0]*cm[0])+np.sum(Clist)/h
Eef=1/unsE
c_0 = np.sqrt(Eef/rhoef)
rhoeff,Eeff=func.propEff(rho[0],cm[0],modulC,modulM,h)

##################
## Energies
##################
Ntmax=int(np.floor(0.25/dt))
# Microstructured medium
energymicro=func.energyMicro(Vfilm,Sfilm,CSVM,Ntmax,dx,rhox,Celx,nbFront,modulC,modulM,frontiere)

# Homogenised medium
energyhomo=func.energyHomogenisedU(V,S,Ntmax,dx,dt,rhoeff,Eeff)

##################
## Make plots
##################

Ntfig=int(np.floor(0.2/dt))
skipSpacePoints=2

psmax=np.where(Ufilm[Ntfig-1,:Nx//2]==max(Ufilm[Ntfig-1,:Nx//2]))
xzoommin=X[psmax[0][0]-200]
xzoommax=X[psmax[0][0]+600]
zoomlim=[xzoommin,xzoommax]

## Displacement fields
plt.figure()
plt.axvspan(xzoommin, xzoommax, facecolor='g', alpha=0.3)
plt.scatter(X[::skipSpacePoints],Ufilm[Ntfig,::skipSpacePoints],marker=r'$\bigcirc$',label=r'$U_h$',color='red',s=30)
plt.plot(X[::skipSpacePoints],U[Ntfig,::skipSpacePoints],linewidth=3,label=r'$U_0$')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.title("$t= $"+str(round(tp[Ntfig], 3))+" s",fontsize=18)
plt.legend(fontsize=16)
plt.tight_layout()
plt.show()

## Zoom 
plt.figure()
for i in range(alpha.size):
    if frontiere["Contact_"+str(i)]!="parfait":
        plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(X[::skipSpacePoints],Ufilm[Ntfig,::skipSpacePoints],marker=r'$\bigcirc$',label=r'$U_h$',color='red',s=30)
plt.plot(X[::skipSpacePoints],U[Ntfig,::skipSpacePoints],linewidth=3,label=r'$U_0$')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.title("$t= $"+str(round(tp[Ntfig], 3))+" s",fontsize=18)
plt.xlim(zoomlim)
plt.legend(fontsize=16)#loc="lower right",
plt.tight_layout()
plt.show()

## Space-time representation
func.sismo(X,tp,Ufilm,U,12,pFilm=1,Nt=2400)
plt.show()


## Dispersion diagram
fctModE="np.sin"
fctModR="np.sin"
NE=8
tau=1/fmlist[0]
Emicro=rho[0]*cm[0]*cm[0]
EpsR=Mlist[0]*DMlist[0]/(Mlist[0]+rho[0]*h)
EpsE=Clist[0]*Emicro*DClist[0]/(Clist[0]*Emicro+h)
if tau!=0:
    omegavect,VP=func.fillMatrix(tau,Eef,rhoef,EpsE,EpsR,fctModE,fctModR,NE,N=16)
    kmax=max(np.real(VP[NE-1,:])/2/np.pi)
    func.dispFK_tperiodic2(Ufilm[:Ntfig,:],dt,dx,fmlist[0],kmax)
    for i in range(NE):
        plt.plot((VP[i,:])/2/np.pi,omegavect/2/np.pi,'--r',linewidth=3)
plt.xlim([0,0.05])
plt.ylim([-fm/2,fm/2])
plt.clim(-9.75,-3.75) 
plt.show()

## Energy evolution
# Without dissipation
plt.figure()
plt.semilogy(tp[:-2],energymicro,label=r'$\mathcal{E}_{h}$',linewidth=3)
plt.semilogy(tp[:-2],energyhomo,label=r'$\mathcal{E}_{0}$',linewidth=3)
plt.xlabel(r'$T$ (s)',fontsize=16)
plt.ylabel(r'$\mathcal{E}$ (J)',fontsize=16)
plt.xlim([0,0.25])
plt.ylim([2e-3,max(energyhomo)])
plt.imshow([1+modula], cmap='gray', extent=[0, tp[-2], 0,max(energyhomo)], aspect='auto', alpha=0.35)
plt.legend(fontsize=16)#loc="lower right",
plt.colorbar().set_label(label=r'$1+\varepsilon_{C,1}~\sin (\Omega_m t)$',size=16)
plt.tight_layout()
plt.show()


# With dissipation
##################
## Load datas
##################
# Raw datas
nSvetp='tp'+str(fs)+'fm'+str(fm)+'.txt'
tp=np.loadtxt(nSvetp,delimiter='\t')
#%No dissipation
nSveResV0='Vfilm0'+str(fs)+'fm'+str(fm)+'.txt'
nSveResS0='Sfilm0'+str(fs)+'fm'+str(fm)+'.txt'
nSveCS0='CS0'+str(fs)+'fm'+str(fm)+'.txt'
nSveResVH0='VfilmSSH0'+str(fs)+'fm'+str(fm)+'.txt'
nSveResSH0='SfilmSSH0'+str(fs)+'fm'+str(fm)+'.txt'

Vfilm0=np.loadtxt(nSveResV0, delimiter='\t')
Sfilm0=np.loadtxt(nSveResS0, delimiter='\t')
CSVM0=np.loadtxt(nSveCS0, delimiter='\t')
VfilmSS0=np.loadtxt(nSveResVH0, delimiter='\t')
SfilmSS0=np.loadtxt(nSveResSH0, delimiter='\t')
#%Low
nSveResVL='VfilmLow'+str(fs)+'fm'+str(fm)+'.txt'
nSveResSL='SfilmLow'+str(fs)+'fm'+str(fm)+'.txt'
nSveCSL='CSLow'+str(fs)+'fm'+str(fm)+'.txt'
nSveResVHL='VfilmSSHLow'+str(fs)+'fm'+str(fm)+'.txt'
nSveResSHL='SfilmSSHLow'+str(fs)+'fm'+str(fm)+'.txt'

VfilmL=np.loadtxt(nSveResVL, delimiter='\t')
SfilmL=np.loadtxt(nSveResSL, delimiter='\t')
CSVML=np.loadtxt(nSveCSL, delimiter='\t')
VfilmSSL=np.loadtxt(nSveResVHL, delimiter='\t')
SfilmSSL=np.loadtxt(nSveResSHL, delimiter='\t')
#%Mid
nSveResVM='VfilmMid'+str(fs)+'fm'+str(fm)+'.txt'
nSveResSM='SfilmMid'+str(fs)+'fm'+str(fm)+'.txt'
nSveCSM='CSMid'+str(fs)+'fm'+str(fm)+'.txt'
nSveResVHM='VfilmSSHMid'+str(fs)+'fm'+str(fm)+'.txt'
nSveResSHM='SfilmSSHMid'+str(fs)+'fm'+str(fm)+'.txt'

VfilmM=np.loadtxt(nSveResVM, delimiter='\t')
SfilmM=np.loadtxt(nSveResSM, delimiter='\t')
CSVMM=np.loadtxt(nSveCSM, delimiter='\t')
VfilmSSM=np.loadtxt(nSveResVHM, delimiter='\t')
SfilmSSM=np.loadtxt(nSveResSHM, delimiter='\t')
#%High
nSveResVH='VfilmHigh'+str(fs)+'fm'+str(fm)+'.txt'
nSveResSH='SfilmHigh'+str(fs)+'fm'+str(fm)+'.txt'
nSveCSH='CSHigh'+str(fs)+'fm'+str(fm)+'.txt'
nSveResVHH='VfilmSSHHigh'+str(fs)+'fm'+str(fm)+'.txt'
nSveResSHH='SfilmSSHHigh'+str(fs)+'fm'+str(fm)+'.txt'

VfilmH=np.loadtxt(nSveResVH, delimiter='\t')
SfilmH=np.loadtxt(nSveResSH, delimiter='\t')
CSVMH=np.loadtxt(nSveCSH, delimiter='\t')
VfilmSSH=np.loadtxt(nSveResVHH, delimiter='\t')
SfilmSSH=np.loadtxt(nSveResSHH, delimiter='\t')
#



# Postprocessing
dt=tp[1]-tp[0]
Ufilm=func.VtoU(Vfilm,dt)
V=func.UtoV(U,dt)
S=func.UtoS(U,dt)

## No dissipation
# Microstructured medium
energymicro0=func.energyMicro(Vfilm0,Sfilm0,CSVM0,Ntmax,dx,rhox,Celx,nbFront,modulC,modulM,frontiere)
# Homogenised medium
energyhomo0=func.energyHomogenisedVS(VfilmSS0,SfilmSS0,Ntmax,dx,dt,rhoeff,Eeff)

## Low dissipation
# Microstructured medium
energymicroL=func.energyMicro(VfilmL,SfilmL,CSVML,Ntmax,dx,rhox,Celx,nbFront,modulC,modulM,frontiere)
# Homogenised medium
energyhomoL=func.energyHomogenisedVS(VfilmSSL,SfilmSSL,Ntmax,dx,dt,rhoeff,Eeff)

## Mid dissipation
# Microstructured medium
energymicroM=func.energyMicro(VfilmM,SfilmM,CSVMM,Ntmax,dx,rhox,Celx,nbFront,modulC,modulM,frontiere)
# Homogenised medium
energyhomoM=func.energyHomogenisedVS(VfilmSSM,SfilmSSM,Ntmax,dx,dt,rhoeff,Eeff)

## Low dissipation
# Microstructured medium
energymicroH=func.energyMicro(VfilmH,SfilmH,CSVMH,Ntmax,dx,rhox,Celx,nbFront,modulC,modulM,frontiere)
# Homogenised medium
energyhomoH=func.energyHomogenisedVS(VfilmSSH,SfilmSSH,Ntmax,dx,dt,rhoeff,Eeff)

plt.figure()
plt.semilogy(tp[:-2],energymicro0,color="grey",label=r'$\mathcal{E}_{h}$, without dissipation',linewidth=3)
plt.semilogy(tp[:-2],energyhomo0,":",color="black",label=r'$\mathcal{E}_{0}$, without dissipation',linewidth=3)
plt.semilogy(tp[:-2],energymicroL,color="salmon",label=r'$\mathcal{E}_{h}$, low dissipation',linewidth=3)
plt.semilogy(tp[:-2],energyhomoL,":",color="firebrick",label=r'$\mathcal{E}_{0}$, low dissipation',linewidth=3)
plt.semilogy(tp[:-2],energymicroM,color="chartreuse",label=r'$\mathcal{E}_{h}$, medium dissipation',linewidth=3)
plt.semilogy(tp[:-2],energyhomoM,":",color="forestgreen",label=r'$\mathcal{E}_{0}$, medium dissipation',linewidth=3)
plt.semilogy(tp[:-2],energymicroH,color="violet",label=r'$\mathcal{E}_{h}$, high dissipation',linewidth=3)
plt.semilogy(tp[:-2],energyhomoH,":",color="darkviolet",label=r'$\mathcal{E}_{0}$, high dissipation',linewidth=3)
plt.xlabel(r'$T$ (s)',fontsize=16)
plt.ylabel(r'$\mathcal{E}$ (J)',fontsize=16)
plt.xlim([0,0.25])
plt.ylim([2e-6,2])
plt.imshow([1+modulC], cmap='gray', extent=[0, tp[-1], 0,max(energyhomo0)], aspect='auto', alpha=0.35)
plt.legend(fontsize=10)#loc="lower right",
plt.colorbar().set_label(label=r'$1+\varepsilon_{C,1}~\sin (\Omega_m t)$',size=16)
plt.tight_layout()
plt.show()
