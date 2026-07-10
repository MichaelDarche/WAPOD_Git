#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:50:43 2023

@author: michael
"""

import os
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/N-Interface_Imerged_PS_modulated_homophiPeriod-Ni-amor/")

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../WAPOD/")

import scheme1D_TDold as scheme1D
import scheme1D_TD_Homo as scheme1DH
import source
import inputReading
import figures as ffg
import matplotlib as mpl


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
################# Input files 
#%
import generatePeriodic
nbcell, nb_milieux, blocs =generatePeriodic.lire_fichier_milieu("MilieuPeriodique.txt")
nbcell, nb_frontieres, blocsF =generatePeriodic.lire_fichier_frontieres("FrontieresPeriodique0.txt")

generatePeriodic.ecrire_milieu("Milieu.txt", nbcell, nb_milieux, blocs)
generatePeriodic.ecrire_frontieres_sortie("Frontiere.txt", nbcell, nb_frontieres, blocsF)


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

nfront=frontiere["nombre de frontieres "]

#% Initial conditions
U0=scheme1D.initPointSourceProblem(configuration,frontiere,mat,sourcef)

###### Mise en oeuvre


Tmaxesp=0.3

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

Nt=int(np.ceil(0.3/dt2p))
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
    if float(frontiere["Kn_"+str(i)])==0:
        Clist.append(0.)
    else:
        Clist.append(1/float(frontiere["Kn_"+str(i)]))
    DClist.append(float(frontiere["DeltaF_"+str(i)]))
    Mlist.append(float(frontiere["Mn_"+str(i)]))
    DMlist.append(float(frontiere["DeltaM_"+str(i)]))
    fmlist.append(float(frontiere["Freq_"+str(i)]))
    if float(frontiere["Freq_"+str(i)])==0:
        taulist.append(0.)
    else:
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

# rhoef0=rhoet+M0/h
# EpsR=DM*M0/h/Nip/rhoef0
# Eef0=Eet*K0*h/(K0*h+Eet)
# EpsE=DK*Eef0/K0/h/Nip




############### MODULATION DIRECT PROPAGATION ##########################################
import time

timef=time.time()
#U,dt2=scheme1DH.temporalSchemePS0MK(Ntp,dt2p,configuration,sourcef,modulation,frontiere,mat)
ResSS0,VfilmSS0,SfilmSS0=scheme1D.FD_sourceVIIHomoSSP(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,xmin,nb_frontieres,pFilm)
timeN_mod=timef-time0
print(timeN_mod)  

time0=time.time()
Res0,Vfilm0,Sfilm0,CSVM0=scheme1D.FD_sourceVII(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm,rCSVM="yes")#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
timef=time.time()
maxV=np.max([np.max(Vfilm0),-np.min(Vfilm0)])
timeN_mod=timef-time0
print(timeN_mod)  
import modulation as modulateP
import homogenize as homoP
# for i in range(nfront):
#     modula=modulateP.fonctionModK(t,alpha[i],modulation,frontiere,i)
#     modulb=modulateP.fonctionModM(t,alpha[i],modulation,frontiere,i)
#     print(nfront)
import sympy as sp

rhoeff0,Eeff0,Ceff0,QCeff0,QMeff0=homoP.propHomoFormP(xmin,mat,frontiere,nb_frontieres)
t = sp.symbols('t')
rhoeffnum0=sp.lambdify(t,rhoeff0)
Ceffnum0=sp.lambdify(t,Ceff0)
# modulb=modulateP.fonctionModMP(t,alpha[0],frontiere,0)
# modulM=M0*(1+modulb)
# dmodulb=np.zeros(int(Nt/pFilm))
# for i in range(int(Nt/pFilm)):
#     modulb[i]=M0*(1+deltaM*np.sin(2*np.pi*frM*t[pFilm*i]))
# #     dmodulb[i]=2*np.pi*frM*M0*(deltaM*np.cos(2*np.pi*frM*t[pFilm*i]))


# modula=modulateP.fonctionModKP(t,alpha[0],frontiere,0)
# modulK=1/K0*(1+modula)
# dmodula=np.zeros(int(Ntfig/pFilm))

energyV=np.zeros(int(Nt2/pFilm))
energyI=np.zeros(int(Nt2/pFilm))
energyhomo0=np.zeros(int(Nt2/pFilm))


Clist=[]
DClist=[]
Mlist=[]
DMlist=[]
fmlist=[]
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


for ti in range(int(Nt2/pFilm)):
    energyV[ti]=dx/2*np.sum(rhox*Vfilm0[ti]*Vfilm0[ti])+dx/2*np.sum(1/rhox/Celx/Celx*Sfilm0[ti]*Sfilm0[ti])
    energyhomo0[ti]=dx/2*np.sum(rhoeffnum0(ti*dt2p)*VfilmSS0[ti]*VfilmSS0[ti])+dx/2*np.sum(1/(rhoeffnum0(ti*dt2p)*Ceffnum0(ti*dt2p)*Ceffnum0(ti*dt2p))*SfilmSS0[ti]*SfilmSS0[ti])
    for interi in range(int(nfront)):
        ind=interi%nb_frontieres
        if ind!=nb_frontieres-1:
            modulM=Mlist[ind]*(1+modulateP.fonctionModMP(ti*dt2p,alpha[0],frontiere,interi))
            modulC=Clist[ind]*(1+modulateP.fonctionModKP(ti*dt2p,alpha[0],frontiere,interi))
            energyI[ti]=energyI[ti]+1/2*modulM*CSVM0[ti,2+4*interi]*CSVM0[ti,2+4*interi]+1/2*modulC*CSVM0[ti,3+4*interi]*CSVM0[ti,3+4*interi]     
energymicro0=energyV+energyI

nbcell, nb_milieux, blocs =generatePeriodic.lire_fichier_milieu("MilieuPeriodique.txt")
nbcell, nb_frontieres, blocsF =generatePeriodic.lire_fichier_frontieres("FrontieresPeriodique1.txt")

generatePeriodic.ecrire_milieu("Milieu.txt", nbcell, nb_milieux, blocs)
generatePeriodic.ecrire_frontieres_sortie("Frontiere.txt", nbcell, nb_frontieres, blocsF)


configuration=inputReading.readfile('Demarrer.txt')
mat=inputReading.readfile('Milieu.txt')
frontiere=inputReading.readfile('Frontiere.txt')     
sourcef=inputReading.readfile('Source.txt')   
#sourceRev=inputReading.readfile('SourceRev.txt')
modulation=inputReading.readfile('Modulation.txt')

import time

timef=time.time()
#U,dt2=scheme1DH.temporalSchemePS0MK(Ntp,dt2p,configuration,sourcef,modulation,frontiere,mat)
ResSS1,VfilmSS1,SfilmSS1=scheme1D.FD_sourceVIIHomoSSP(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,xmin,nb_frontieres,pFilm)
timeN_mod=timef-time0
print(timeN_mod)  

time0=time.time()
Res1,Vfilm1,Sfilm1,CSVM1=scheme1D.FD_sourceVII(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm,rCSVM="yes")#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
timef=time.time()
maxV=np.max([np.max(Vfilm1),-np.min(Vfilm1)])
timeN_mod=timef-time0
print(timeN_mod)  

rhoeff1,Eeff1,Ceff1,QCeff1,QMeff1=homoP.propHomoFormP(xmin,mat,frontiere,nb_frontieres)
t = sp.symbols('t')
rhoeffnum1=sp.lambdify(t,rhoeff1)
Ceffnum1=sp.lambdify(t,Ceff1)

energyV=np.zeros(int(Nt2/pFilm))
energyI=np.zeros(int(Nt2/pFilm))
energyhomo1=np.zeros(int(Nt2/pFilm))


Clist=[]
DClist=[]
Mlist=[]
DMlist=[]
fmlist=[]
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


for ti in range(int(Nt2/pFilm)):
    energyV[ti]=dx/2*np.sum(rhox*Vfilm1[ti]*Vfilm1[ti])+dx/2*np.sum(1/rhox/Celx/Celx*Sfilm1[ti]*Sfilm1[ti])
    energyhomo1[ti]=dx/2*np.sum(rhoeffnum1(ti*dt2p)*VfilmSS1[ti]*VfilmSS1[ti])+dx/2*np.sum(1/(rhoeffnum1(ti*dt2p)*Ceffnum1(ti*dt2p)*Ceffnum1(ti*dt2p))*SfilmSS1[ti]*SfilmSS1[ti])
    for interi in range(int(nfront)):
        ind=interi%nb_frontieres
        if ind!=nb_frontieres-1:
            modulM=Mlist[ind]*(1+DMlist[ind]*np.sin(2*np.pi*fmlist[ind]*ti*dt2p))
            modulK=Clist[ind]*(1+DClist[ind]*np.sin(2*np.pi*fmlist[ind]*ti*dt2p))
            energyI[ti]=energyI[ti]+1/2*modulM*CSVM1[ti,2+4*interi]*CSVM1[ti,2+4*interi]+1/2*modulK*CSVM1[ti,3+4*interi]*CSVM1[ti,3+4*interi]     
energymicro1=energyV+energyI



nbcell, nb_milieux, blocs =generatePeriodic.lire_fichier_milieu("MilieuPeriodique.txt")
nbcell, nb_frontieres, blocsF =generatePeriodic.lire_fichier_frontieres("FrontieresPeriodique2.txt")

generatePeriodic.ecrire_milieu("Milieu.txt", nbcell, nb_milieux, blocs)
generatePeriodic.ecrire_frontieres_sortie("Frontiere.txt", nbcell, nb_frontieres, blocsF)


configuration=inputReading.readfile('Demarrer.txt')
mat=inputReading.readfile('Milieu.txt')
frontiere=inputReading.readfile('Frontiere.txt')     
sourcef=inputReading.readfile('Source.txt')   
#sourceRev=inputReading.readfile('SourceRev.txt')
modulation=inputReading.readfile('Modulation.txt')

timef=time.time()
#U,dt2=scheme1DH.temporalSchemePS0MK(Ntp,dt2p,configuration,sourcef,modulation,frontiere,mat)
ResSS2,VfilmSS2,SfilmSS2=scheme1D.FD_sourceVIIHomoSSP(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,xmin,nb_frontieres,pFilm)
timeN_mod=timef-time0
print(timeN_mod)  

time0=time.time()
Res2,Vfilm2,Sfilm2,CSVM2=scheme1D.FD_sourceVII(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm,rCSVM="yes")#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
timef=time.time()
maxV=np.max([np.max(Vfilm2),-np.min(Vfilm2)])
timeN_mod=timef-time0
print(timeN_mod)  

rhoeff2,Eeff2,Ceff2,QCeff2,QMeff2=homoP.propHomoFormP(xmin,mat,frontiere,nb_frontieres)
t = sp.symbols('t')
rhoeffnum2=sp.lambdify(t,rhoeff2)
Ceffnum2=sp.lambdify(t,Ceff2)



energyV=np.zeros(int(Nt2/pFilm))
energyI=np.zeros(int(Nt2/pFilm))
energyhomo2=np.zeros(int(Nt2/pFilm))


Clist=[]
DClist=[]
Mlist=[]
DMlist=[]
fmlist=[]
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


for ti in range(int(Nt2/pFilm)):
    energyV[ti]=dx/2*np.sum(rhox*Vfilm2[ti]*Vfilm2[ti])+dx/2*np.sum(1/rhox/Celx/Celx*Sfilm2[ti]*Sfilm2[ti])
    energyhomo2[ti]=dx/2*np.sum(rhoeffnum2(ti*dt2p)*VfilmSS2[ti]*VfilmSS2[ti])+dx/2*np.sum(1/(rhoeffnum2(ti*dt2p)*Ceffnum2(ti*dt2p)*Ceffnum2(ti*dt2p))*SfilmSS2[ti]*SfilmSS2[ti])
    for interi in range(int(nfront)):
        ind=interi%nb_frontieres
        if ind!=nb_frontieres-1:
            modulM=Mlist[ind]*(1+DMlist[ind]*np.sin(2*np.pi*fmlist[ind]*ti*dt2p))
            modulK=Clist[ind]*(1+DClist[ind]*np.sin(2*np.pi*fmlist[ind]*ti*dt2p))
            energyI[ti]=energyI[ti]+1/2*modulM*CSVM2[ti,2+4*interi]*CSVM2[ti,2+4*interi]+1/2*modulK*CSVM2[ti,3+4*interi]*CSVM2[ti,3+4*interi]     
energymicro2=energyV+energyI

nbcell, nb_milieux, blocs =generatePeriodic.lire_fichier_milieu("MilieuPeriodique.txt")
nbcell, nb_frontieres, blocsF =generatePeriodic.lire_fichier_frontieres("FrontieresPeriodique3.txt")

generatePeriodic.ecrire_milieu("Milieu.txt", nbcell, nb_milieux, blocs)
generatePeriodic.ecrire_frontieres_sortie("Frontiere.txt", nbcell, nb_frontieres, blocsF)


configuration=inputReading.readfile('Demarrer.txt')
mat=inputReading.readfile('Milieu.txt')
frontiere=inputReading.readfile('Frontiere.txt')     
sourcef=inputReading.readfile('Source.txt')   
#sourceRev=inputReading.readfile('SourceRev.txt')
modulation=inputReading.readfile('Modulation.txt')

timef=time.time()
#U,dt2=scheme1DH.temporalSchemePS0MK(Ntp,dt2p,configuration,sourcef,modulation,frontiere,mat)
ResSS3,VfilmSS3,SfilmSS3=scheme1D.FD_sourceVIIHomoSSP(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,xmin,nb_frontieres,pFilm)
timeN_mod=timef-time0
print(timeN_mod)  

time0=time.time()
Res3,Vfilm3,Sfilm3,CSVM3=scheme1D.FD_sourceVII(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm,rCSVM="yes")#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
timef=time.time()
maxV=np.max([np.max(Vfilm3),-np.min(Vfilm3)])
timeN_mod=timef-time0
print(timeN_mod)  

rhoeff3,Eeff3,Ceff3,QCeff3,QMeff3=homoP.propHomoFormP(xmin,mat,frontiere,nb_frontieres)
t = sp.symbols('t')
rhoeffnum3=sp.lambdify(t,rhoeff3)
Ceffnum3=sp.lambdify(t,Ceff3)



energyV=np.zeros(int(Nt2/pFilm))
energyI=np.zeros(int(Nt2/pFilm))
energyhomo3=np.zeros(int(Nt2/pFilm))


Clist=[]
DClist=[]
Mlist=[]
DMlist=[]
fmlist=[]
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


for ti in range(int(Nt2/pFilm)):
    energyV[ti]=dx/2*np.sum(rhox*Vfilm3[ti]*Vfilm3[ti])+dx/2*np.sum(1/rhox/Celx/Celx*Sfilm3[ti]*Sfilm3[ti])
    energyhomo3[ti]=dx/2*np.sum(rhoeffnum3(ti*dt2p)*VfilmSS3[ti]*VfilmSS3[ti])+dx/2*np.sum(1/(rhoeffnum3(ti*dt2p)*Ceffnum3(ti*dt2p)*Ceffnum3(ti*dt2p))*SfilmSS3[ti]*SfilmSS3[ti])
    for interi in range(int(nfront)):
        ind=interi%nb_frontieres
        if ind!=nb_frontieres-1:
            modulM=Mlist[ind]*(1+DMlist[ind]*np.sin(2*np.pi*fmlist[ind]*ti*dt2p))
            modulK=Clist[ind]*(1+DClist[ind]*np.sin(2*np.pi*fmlist[ind]*ti*dt2p))
            energyI[ti]=energyI[ti]+1/2*modulM*CSVM3[ti,2+4*interi]*CSVM3[ti,2+4*interi]+1/2*modulK*CSVM3[ti,3+4*interi]*CSVM3[ti,3+4*interi]     
energymicro3=energyV+energyI

nbcell, nb_milieux, blocs =generatePeriodic.lire_fichier_milieu("MilieuPeriodique.txt")
nbcell, nb_frontieres, blocsF =generatePeriodic.lire_fichier_frontieres("FrontieresPeriodique4.txt")

generatePeriodic.ecrire_milieu("Milieu.txt", nbcell, nb_milieux, blocs)
generatePeriodic.ecrire_frontieres_sortie("Frontiere.txt", nbcell, nb_frontieres, blocsF)


configuration=inputReading.readfile('Demarrer.txt')
mat=inputReading.readfile('Milieu.txt')
frontiere=inputReading.readfile('Frontiere.txt')     
sourcef=inputReading.readfile('Source.txt')   
#sourceRev=inputReading.readfile('SourceRev.txt')
modulation=inputReading.readfile('Modulation.txt')

timef=time.time()
#U,dt2=scheme1DH.temporalSchemePS0MK(Ntp,dt2p,configuration,sourcef,modulation,frontiere,mat)
ResSS4,VfilmSS4,SfilmSS4=scheme1D.FD_sourceVIIHomoSSP(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,xmin,nb_frontieres,pFilm)
timeN_mod=timef-time0
print(timeN_mod)  

time0=time.time()
Res4,Vfilm4,Sfilm4,CSVM4=scheme1D.FD_sourceVII(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm,rCSVM="yes")#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
timef=time.time()
maxV=np.max([np.max(Vfilm4),-np.min(Vfilm4)])
timeN_mod=timef-time0
print(timeN_mod)  

rhoeff4,Eeff4,Ceff4,QCeff4,QMeff4=homoP.propHomoFormP(xmin,mat,frontiere,nb_frontieres)
t = sp.symbols('t')
rhoeffnum4=sp.lambdify(t,rhoeff4)
Ceffnum4=sp.lambdify(t,Ceff4)



energyV=np.zeros(int(Nt2/pFilm))
energyI=np.zeros(int(Nt2/pFilm))
energyhomo4=np.zeros(int(Nt2/pFilm))


Clist=[]
DClist=[]
Mlist=[]
DMlist=[]
fmlist=[]
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


for ti in range(int(Nt2/pFilm)):
    energyV[ti]=dx/2*np.sum(rhox*Vfilm4[ti]*Vfilm4[ti])+dx/2*np.sum(1/rhox/Celx/Celx*Sfilm4[ti]*Sfilm4[ti])
    energyhomo4[ti]=dx/2*np.sum(rhoeffnum4(ti*dt2p)*VfilmSS4[ti]*VfilmSS4[ti])+dx/2*np.sum(1/(rhoeffnum4(ti*dt2p)*Ceffnum4(ti*dt2p)*Ceffnum4(ti*dt2p))*SfilmSS4[ti]*SfilmSS4[ti])
    for interi in range(int(nfront)):
        ind=interi%nb_frontieres
        if ind!=nb_frontieres-1:
            modulM=Mlist[ind]*(1+DMlist[ind]*np.sin(2*np.pi*fmlist[ind]*ti*dt2p))
            modulK=Clist[ind]*(1+DClist[ind]*np.sin(2*np.pi*fmlist[ind]*ti*dt2p))
            energyI[ti]=energyI[ti]+1/2*modulM*CSVM4[ti,2+4*interi]*CSVM4[ti,2+4*interi]+1/2*modulK*CSVM4[ti,3+4*interi]*CSVM4[ti,3+4*interi]     
energymicro4=energyV+energyI

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

total_frames = Vfilm2.shape[0]
fig, ax = plt.subplots()

# # Cree l'animationen certains pointspyten certains pointspythonen certains pointspythonen certains pointspythonen certains pointspythonen certains pointspythonen certains pointspythonhon
# animation = FuncAnimation(fig, update, frames=total_frames, interval=1)
# # Affiche l'animation
# plt.show()
# nomAnim="Film_SM_"+sourcef["Frequence"]+"Hz.mp4"
# animation.save(nomAnim, writer='ffmpeg', fps=20)

# energyM=np.zeros(int(Nt/pFilm))
# for ti in range(int(Nt/pFilm)):
#     energyM[ti]=1/2*np.sum(rhox*Vfilm2[ti]*Vfilm[ti])+1/2*np.sum(1/rhox/Celx/Celx*Sfilm[ti]*Sfilm[ti])


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
# timef=time.time()
# U,dt2=scheme1DH.temporalSchemePS0MK(Ntp,dt2p,configuration,sourcef,modulation,frontiere,mat)
# #ResSS,VfilmSS,SfilmSS=scheme1D.FD_sourceVIIHomoSSP(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,xmin,nb_frontieres,pFilm)
# timeN_mod=timef-time0
# print(timeN_mod)  
# # # ###################### NO MODULATION #################################
# time0=time.time()
# ResNM,VfilmNM,SfilmNM=scheme1D.FD_sourceVII(U0,Ntp,Nx,K,rhox,Celx,dt2p,dx,x_0,sce,configuration,mat,frontiere,nomodulation,pFilm)
# timeNM=time.time()-time0
# print(timeNM)  
Ufilm2=ffg.VtoU(Vfilm2,dt2p)
UfilmSS2=ffg.VtoU(VfilmSS2,dt2p)
maxU=np.max(Ufilm2)


Ufilm0=ffg.VtoU(Vfilm0,dt2p)
UfilmSS0=ffg.VtoU(VfilmSS0,dt2p)
maxU=np.max(Ufilm2)

Ufilm1=ffg.VtoU(Vfilm1,dt2p)
UfilmSS1=ffg.VtoU(VfilmSS1,dt2p)
maxU=np.max(Ufilm2)


Ufilm1b=ffg.VtoU2(Vfilm1,dt2p)
UfilmSS1b=ffg.VtoU2(VfilmSS1,dt2p)
maxU=np.max(Ufilm2)

Ufilm1t=ffg.VtoU3(Vfilm1,dt2p)
UfilmSS1t=ffg.VtoU3(VfilmSS1,dt2p)
maxU=np.max(Ufilm2)

Ufilm1q=ffg.VtoU4(Vfilm1,dt2p)
UfilmSS1q=ffg.VtoU4(VfilmSS1,dt2p)
maxU=np.max(Ufilm2)




Ufilm3=ffg.VtoU(Vfilm3,dt2p)
UfilmSS3=ffg.VtoU(VfilmSS3,dt2p)
maxU=np.max(Ufilm2)

Ufilm4=ffg.VtoU(Vfilm4,dt2p)
UfilmSS4=ffg.VtoU(VfilmSS4,dt2p)
maxU=np.max(Ufilm2)


nVfilm0='Vfilm0'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nVfilm0, Vfilm0,fmt='%18e', delimiter='\t')
nCS0='CS0'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nCS0, CSVM0,fmt='%18e', delimiter='\t')
nVfilmSS0='VfilmSSH0'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nVfilmSS0, VfilmSS0,fmt='%18e', delimiter='\t')

nVfilmL='VfilmLow'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nVfilmL, Vfilm4,fmt='%18e', delimiter='\t')
nCSL='CSLow'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nCSL, CSVM4,fmt='%18e', delimiter='\t')
nVfilmSSL='VfilmSSHLow'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nVfilmSSL, VfilmSS4,fmt='%18e', delimiter='\t')

nVfilmM='VfilmMid'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nVfilmM, Vfilm1,fmt='%18e', delimiter='\t')
nCSM='CSMid'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nCSM, CSVM1,fmt='%18e', delimiter='\t')
nVfilmSSM='VfilmSSHMid'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nVfilmSSM, VfilmSS1,fmt='%18e', delimiter='\t')

nVfilmH='VfilmHigh'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nVfilmH, Vfilm2,fmt='%18e', delimiter='\t')
nCSH='CSHigh'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nCSH, CSVM2,fmt='%18e', delimiter='\t')
nVfilmSSH='VfilmSSHHigh'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nVfilmSSH, VfilmSS2,fmt='%18e', delimiter='\t')

nSfilm0='Sfilm0'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nSfilm0, Sfilm0,fmt='%18e', delimiter='\t')
nSfilmSS0='SfilmSSH0'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nSfilmSS0, SfilmSS0,fmt='%18e', delimiter='\t')

nSfilmL='SfilmLow'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nSfilmL, Sfilm4,fmt='%18e', delimiter='\t')
nSfilmSSL='SfilmSSHLow'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nSfilmSSL, SfilmSS4,fmt='%18e', delimiter='\t')

nSfilmM='SfilmMid'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nSfilmM, Sfilm1,fmt='%18e', delimiter='\t')
nSfilmSSM='SfilmSSHMid'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nSfilmSSM, SfilmSS1,fmt='%18e', delimiter='\t')

nSfilmH='SfilmHigh'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nSfilmH, Sfilm2,fmt='%18e', delimiter='\t')
nSfilmSSH='SfilmSSHHigh'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(nSfilmSSH, SfilmSS2,fmt='%18e', delimiter='\t')


dt2='dt2'+str(fs)+'fm'+str(fmlist[0])+'.txt'
np.savetxt(dt2, dt2, fmt='%18e', delimiter='\t')


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
    plt.scatter(X[::1],Ufilm1t[frame,::1],marker=r'$\bigcirc$',label=r'$U_h$',color='red',s=30)
    # plt.plot(X[::1],Ufilm[Ntfig,::1],linestyle=None,marker='o',color='white',markersize=2)
    plt.plot(X,UfilmSS1t[frame,:],linewidth=3,label=r'$U_0$')
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
    plt.ylabel(r'$U$ (m/s)',fontsize=16)
    plt.ylim([-0.0000050,0.000005])
    plt.tight_layout()
    # plt.legend()


total_frames = VfilmSS1.shape[0]
fig, ax = plt.subplots()

# Cree l'animation
animationComp = FuncAnimation(fig, updateComp, frames=total_frames, interval=1)
plt.show()
nomAnim="Film_HighDissp_"+sourcef["Frequence"]+frontiere["Freq_0"]+"Hz.mp4"
animationComp.save(nomAnim, writer='ffmpeg', fps=20)

import modulation as modul


Ntfig=int(np.floor(0.2/dt2p))
maxVt=np.max([np.max(Vfilm2[Ntfig,:]),-np.min(Vfilm2[Ntfig,:])])

psmax=np.where(Ufilm2[Ntfig-1,:Nx//2]==max(Ufilm2[Ntfig-1,:Nx//2]))
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

# Cree l'animation
#animation = FuncAnimation(fig, update, frames=total_frames, interval=1)
# Affiche l'animation
# plt.show()
# def updateComp(frame):
#     plt.clf()
#     vals = []
#     for i in range(alpha.size-1):
#         if frontiere["Contact_"+str(i)] != "parfait":
#             v = modul.fonctionModMP(frame*dt2p, alpha[i], frontiere, i)
#             vals.append(v)
#         else:
#             vals.append(0)
#     vals = np.array(vals)
#     for i in range(alpha.size):
#         if frontiere["Contact_"+str(i)]!="parfait":
#             plt.axvline(x=alpha[i],ymin=0,ymax=1,color=cmap(vals[i]),alpha=np.abs(vals[i]/2),linewidth=0.5)#,linestyle='--')
#             # elif modul.fonctionModMP(Ntfig*dt2p,alpha[i],frontiere,i)>0.:
#             #     plt.axvline(x=alpha[i],ymin=0,ymax=1,color='orchid',alpha=modul.fonctionModMP(Ntfig*dt2p,alpha[i]/2,frontiere,i))
#             #plt.colorbar()#,linestyle='--')
#     sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
#     sm.set_array([])
#     plt.scatter(X[::1],Ufilm2[frame,::1],marker=r'$\bigcirc$',label=r'$U_h$',color='red',s=30)
#     # plt.plot(X[::1],Ufilm[Ntfig,::1],linestyle=None,marker='o',color='white',markersize=2)
#     plt.plot(X,UfilmSS2[frame,:],linewidth=3,label=r'$U_0$')
#     # plt.plot(X, VfilmH[frame,:],"--")#, label=f'Frame {frame}')
#     # plt.plot(X, VfilmNM[frame,:],"o")#, label=f'Frame {frame}')
#     # plt.imshow([rhox], cmap='gray', extent=[xmin, xmax, -0.0015,0.0015], aspect='auto', alpha=0.35)
#     # plt.colorbar()
#     # plt.plot(X, VfilmH[frame,:])
#     plt.xlabel("X (m)",fontsize=16)
#     plt.ylabel("U (m)",fontsize=16)
#     plt.ylim([-0.000005,0.000005])
#     plt.tight_layout()
#     # plt.legend()


# total_frames = VfilmSS2.shape[0]
# fig, ax = plt.subplots()

# # Cree l'animation
# animationComp = FuncAnimation(fig, updateComp, frames=total_frames, interval=1)
# plt.show()
# nomAnim="Film_1IFB2_"+sourcef["Frequence"]+modulation["Freq"]+"Hz.mp4"
# animationComp.save(nomAnim, writer='ffmpeg', fps=20)



plt.figure()
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
import modulation as modulateP
import homogenize as homoP
# for i in range(nfront):
#     modula=modulateP.fonctionModK(t,alpha[i],modulation,frontiere,i)
#     modulb=modulateP.fonctionModM(t,alpha[i],modulation,frontiere,i)
#     print(nfront)
# import sympy as sp
# rhoeff,Eeff,Ceff,QCeff,QMeff=homoP.propHomoFormP(xmin,mat,frontiere,nb_frontieres)
# t = sp.symbols('t')
# rhoeffnum=sp.lambdify(t,rhoeff)
# Ceffnum=sp.lambdify(t,Ceff)
# # modulb=modulateP.fonctionModMP(t,alpha[0],frontiere,0)
# # modulM=M0*(1+modulb)
# # dmodulb=np.zeros(int(Nt/pFilm))
# # for i in range(int(Nt/pFilm)):
# #     modulb[i]=M0*(1+deltaM*np.sin(2*np.pi*frM*t[pFilm*i]))
# # #     dmodulb[i]=2*np.pi*frM*M0*(deltaM*np.cos(2*np.pi*frM*t[pFilm*i]))


# # modula=modulateP.fonctionModKP(t,alpha[0],frontiere,0)
# # modulK=1/K0*(1+modula)
# # dmodula=np.zeros(int(Ntfig/pFilm))

# # energyV=np.zeros(int(Nt2/pFilm))
# # energyI=np.zeros(int(Nt2/pFilm))
# energyhomo2=np.zeros(int(Nt2/pFilm))

# for ti in range(int(Nt2/pFilm)):
#     # energyV[ti]=dx/2*np.sum(rhox*Vfilm[ti]*Vfilm[ti])+dx/2*np.sum(1/rhox/Celx/Celx*Sfilm[ti]*Sfilm[ti])
#     energyhomo2[ti]=dx/2*np.sum(rhoeffnum(ti*dt2p)*VfilmSS2[ti]*VfilmSS2[ti])+dx/2*np.sum(1/(rhoeffnum(ti*dt2p)*Ceffnum(ti*dt2p)*Ceffnum(ti*dt2p))*SfilmSS2[ti]*SfilmSS2[ti])
#     # for interi in range(int(nfront)):
#         # energyI[ti]=energyI[ti]+1/2*modulM[ti]*CSVM[ti,2+4*interi]*CSVM[ti,2+4*interi]+1/2*modulK[ti]*CSVM[ti,3+4*interi]*CSVM[ti,3+4*interi]     
# # energymicro=energyV+energyI

fct=modulateP.fonctionModKP(tp[:-1],0,frontiere,0)

plt.figure()
# plt.semilogy(tp[:-1],energymicro,label=r'$\mathcal{E}_{h}$',linewidth=3)
plt.semilogy(tp[:-1],energymicro0,color="grey",label=r'$\mathcal{E}_{h}$, without dissipation',linewidth=3)
plt.semilogy(tp[:-1],energyhomo0,":",color="black",label=r'$\mathcal{E}_{0}$, without dissipation',linewidth=3)
plt.semilogy(tp[:-1],energymicro4,color="salmon",label=r'$\mathcal{E}_{h}$, low dissipation',linewidth=3)
plt.semilogy(tp[:-1],energyhomo4,":",color="firebrick",label=r'$\mathcal{E}_{0}$, low dissipation',linewidth=3)
plt.semilogy(tp[:-1],energymicro1,color="chartreuse",label=r'$\mathcal{E}_{h}$, medium dissipation',linewidth=3)
plt.semilogy(tp[:-1],energyhomo1,":",color="forestgreen",label=r'$\mathcal{E}_{0}$, medium dissipation',linewidth=3)
plt.semilogy(tp[:-1],energymicro2,color="violet",label=r'$\mathcal{E}_{h}$, high dissipation',linewidth=3)
plt.semilogy(tp[:-1],energyhomo2,":",color="darkviolet",label=r'$\mathcal{E}_{0}$, high dissipation',linewidth=3)
# plt.semilogy(tp[:-1],energyhomo3,label=r'$\mathcal{E}_{0}$, high dissipation',linewidth=3)
plt.xlabel(r'$T$ (s)',fontsize=16)
plt.ylabel(r'$\mathcal{E}$ (J)',fontsize=16)
plt.xlim([0,0.25])
# plt.xlim([0,tp[-1]])
plt.ylim([2e-6,2])
plt.imshow([1+fct], cmap='gray', extent=[0, tp[-1], 0,max(energyhomo0)], aspect='auto', alpha=0.35)
# plt.title("$t= $"+str(round(tp[Nt2], 3))+" s",fontsize=18)
plt.legend(fontsize=10)#loc="lower right",
plt.colorbar().set_label(label=r'$1+\varepsilon_{C,1}~\sin (\Omega_m t)$',size=16)
plt.tight_layout()

ffg.sismo(X, t, VfilmSS,20,pFilm,Nt=500)

Ntfig=int(np.floor(0.15/dt2p))
maxVt=np.max([np.max(Vfilm[Ntfig,:]),-np.min(Vfilm[Ntfig,:])])

psmax=np.where(Ufilm0[Ntfig-1,:Nx//2]==max(Ufilm0[Ntfig-1,:Nx//2]))
xzoommin=X[psmax[0][0]-200]
xzoommax=X[psmax[0][0]+600]
zoomlim=[xzoommin,xzoommax]

plt.figure()
for i in range(alpha.size):
    plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(X[::1],Ufilm0[Ntfig,::1],marker=r'$\bigcirc$',label=r'$U_h$',color='blue',s=40)
#plt.plot(X[::1],Ufilm[Ntfig,::1],linestyle='None',marker='o',color='white',markersize=3)
#plt.plot(X,Uz[Ntfig-1,:],'grey',linewidth=3,label=r'$U_0$')
plt.plot(X,UfilmSS0[Ntfig,:],color='red',linestyle='None',marker='o',markersize=3,label=r'$U_0$')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.title("$T=$"+str(round(tp[Ntfig], 3))+" s",fontsize=16)
plt.legend(fontsize=14)#loc="lower right",)
plt.tight_layout()


plt.figure()
plt.axvspan(xzoommin, xzoommax, facecolor='g', alpha=0.3)
# for i in range(alpha.size):
#     plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(X[::1],Ufilm2[Ntfig,::1],marker=r'$\bigcirc$',label=r'$U_h$',color='red',s=30)
plt.plot(X[::1],Ufilm2[Ntfig,::1],linestyle=None,marker='o',color='white',markersize=2)
plt.plot(X,UfilmSS2[Ntfig,:],linewidth=3,label=r'$U_0$')
# plt.axvline(x=alpha[0],color='black')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.title("$t= $"+str(round(tp[Ntfig], 3))+" s",fontsize=18)
plt.legend(fontsize=16)#loc="lower right",
plt.tight_layout()
namefig="Figures/PRSA_U_MK-amorKMQCQM_"+sourcef["Frequence"]+"_"+frontiere["Freq_0"]+".pdf"
plt.savefig(str(namefig), format='pdf')
namefig="Figures/PRSA_U_MK-amorKMQCQM_"+sourcef["Frequence"]+"_"+frontiere["Freq_0"]+".png"
plt.savefig(str(namefig), format='png')

plt.figure()
plt.axvspan(xzoommin, xzoommax, facecolor='g', alpha=0.3)
# for i in range(alpha.size):
#     plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(X[::1],Vfilm[Ntfig,::1],marker=r'$\bigcirc$',label=r'$V_h$',color='red',s=30)
plt.plot(X[::1],Vfilm[Ntfig,::1],linestyle=None,marker='o',color='white',markersize=2)
plt.plot(X,VfilmSS[Ntfig,:],linewidth=3,label=r'$V_0$')
# plt.axvline(x=alpha[0],color='black')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$V$ (m/s)',fontsize=16)
plt.title("$t= $"+str(round(tp[Ntfig], 3))+" s",fontsize=18)
plt.legend(fontsize=16)#loc="lower right",
plt.tight_layout()
# plt.plot(VfilmSS[500,:])
# plt.plot(VfilmH[500,:])
# plt.title("$nt=$500")
namefig="Figures/PRSA_V_MK-amorKMQCQM_"+sourcef["Frequence"]+"_"+frontiere["Freq_0"]+".pdf"
plt.savefig(str(namefig), format='pdf')
namefig="Figures/PRSA_V_MK-amorKMQCQM_"+sourcef["Frequence"]+"_"+frontiere["Freq_0"]+".png"
plt.savefig(str(namefig), format='png')



plt.figure()
for i in range(alpha.size):
    plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(X[::1],Ufilm[Ntfig,::1],marker=r'$\bigcirc$',label=r'$U_h$',color='red',s=30)
plt.plot(X[::1],Ufilm[Ntfig,::1],linestyle=None,marker='o',color='white',markersize=2)
plt.plot(X,UfilmSS[Ntfig,:],linewidth=3,label=r'$U_0$')
# plt.axvline(x=alpha[0],color='black')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.title("$t= $"+str(round(tp[Ntfig], 3))+" s",fontsize=18)
plt.xlim(zoomlim)
plt.legend(fontsize=16)#loc="lower right",
plt.tight_layout()
namefig="Figures/PRSAzoom_U_MK-amorKMQCQM_"+sourcef["Frequence"]+"_"+frontiere["Freq_0"]+".pdf"
plt.savefig(str(namefig), format='pdf')
namefig="Figures/PRSAzoom_U_MK-amorKMQCQM_"+sourcef["Frequence"]+"_"+frontiere["Freq_0"]+".png"
plt.savefig(str(namefig), format='png')

plt.figure()
for i in range(alpha.size):
    plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(X[::1],Vfilm[Ntfig,::1],marker=r'$\bigcirc$',label=r'$V_h$',color='red',s=30)
plt.plot(X[::1],Vfilm[Ntfig,::1],linestyle=None,marker='o',color='white',markersize=2)
plt.plot(X,VfilmSS[Ntfig,:],linewidth=3,label=r'$V_0$')
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
namefig="Figures/PRSAzoom_V_MK-amorKMQCQM_"+sourcef["Frequence"]+"_"+frontiere["Freq_0"]+".pdf"
plt.savefig(str(namefig), format='pdf')
namefig="Figures/PRSAzoom_V_MK-amorKMQCQM_"+sourcef["Frequence"]+"_"+frontiere["Freq_0"]+".png"
plt.savefig(str(namefig), format='png')

plt.figure()
plt.plot(Vfilm[1000,:])
# plt.plot(VfilmGS[1000,:])
plt.plot(VfilmSS[1000,:])
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
