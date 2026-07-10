#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:50:43 2023

@author: michael
"""
import os
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/Homogeneous_Cauchy_Leo/")

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
################# Parameters of the simulation   
X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(configuration)
fs,CFL,tshift=inputReading.source(sourcef)
Nmat,rho,cm=inputReading.material(mat)
alpha=inputReading.frontiere(frontiere)
rhox,Celx=inputReading.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
dt=min(scheme1D.timeStep(CFL,dx,cm))
lambdaOnde=cm/fs
xsmin=cm*tshift-lambdaOnde
xsmax=cm*tshift

#%

#% Initial conditions
U0=scheme1D.initCauchyProblem(configuration,frontiere,mat,sourcef)

###### Mise en oeuvre
def errL2(A,B):
    err=0
    B2=0
    L=max(B.shape)
    for l in range(L):
        eps=A[l]-B[l]
        err=err+eps**2
        B2=B2+B[l]**2
    return err,err/B2

ndx=7
# err=np.zeros([4,ndx])
errIm=np.zeros([6,ndx])
# errR=np.zeros([4,ndx])
errRIm=np.zeros([6,ndx])
data=np.zeros([4,ndx])
# k0=2**((ndx-1)/2)
# mat=inputReading.readfile('Milieu.txt')
# frontiere=inputReading.readfile('Frontiere.txt')     
# sourcef=inputReading.readfile('Source.txt')   
# sourceRev=inputReading.readfile('SourceRev.txt')
for k in range(0,ndx):
    K=4
    nDem='Demarrer'+str(k)+'.txt'
    configuration=inputReading.readfile(nDem)
    mat=inputReading.readfile('Milieu.txt')
    frontiere=inputReading.readfile('Frontiere.txt')     
    sourcef=inputReading.readfile('Source.txt')   
    modulation0=inputReading.readfile('Modulation.txt')
    # modulation100=inputReading.readfile('Modulation100.txt')
    # modulation500=inputReading.readfile('Modulation500.txt')
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
    Ntfig=int(np.floor(0.045/dt))
    Tmax=Ntfig*dt
    tmax=str(Tmax)+"s"
    t=np.linspace(0,Tmax,Ntfig+1)
    sce=np.zeros(Ntfig+1)
    for i in range(Ntfig+1):
        sce[i]=source.choice_timefct(sourcef,t[i])
    THVTMy=np.zeros(Nx)
    THSTMy=np.zeros(Nx)
    THVTMyF=np.zeros(Nx)
    THSTMyF=np.zeros(Nx)
    Nxalpha=int(alpha[0]/dx)
    WL=max(cm/fs)
    data[0,k]=Tmax
    data[1,k]=dx
    data[2,k]=dt
    data[3,k]=WL/dx 
    # NFourier=128
    # Kf=10
    # df=Kf*fs/NFourier
    # for pt in range(Nx):
    #     THVTMyF[pt],THSTMyF[pt]=analytics.analytics1InterfaceCauchy(tshift+Tmax,X[pt],sourcef,frontiere,NFourier,mat,Kf)
    # SolV,SolS=analytics.analytics_Cauchy_1IModM(Tmax,X,sourcef,modulation,frontiere,mat,NOde,dt,Ntfig+1)
    # R01=(rho[1]*cm[1]-rho[0]*cm[0])/(rho[1]*cm[1]+rho[0]*cm[0])
    # T01=2*rho[0]*cm[1]/(rho[1]*cm[1]+rho[0]*cm[0])
    # Nxalpha=int(alpha[0]/dx)
    # for x in range(Nxalpha):
    #     f0=source.choice_timefct(sourcef,tshift-x*dx/cm[0])
    #     fr=source.choice_timefct(sourcef,tshift+x*dx/cm[0]-2*alpha[0]/cm[0])
    #     fT=source.choice_timefct(sourcef,Tmax+tshift-x*dx/cm[0])
    #     frT=source.choice_timefct(sourcef,Tmax+tshift+x*dx/cm[0]-2*alpha[0]/cm[0])
    #     THVTMy[x]=1/cm[0]*fT-R01/cm[0]*frT
    #     THSTMy[x]=-rho[0]*fT-R01*rho[0]*frT
    # for x in range(Nxalpha,Nx):
    #     ft=source.choice_timefct(sourcef,tshift-x*dx/cm[1]+(1/cm[1]-1/cm[0])*alpha[0])
    #     ftT=source.choice_timefct(sourcef,Tmax+tshift-x*dx/cm[1]+(1/cm[1]-1/cm[0])*alpha[0])
    #     THVTMy[x]=1/cm[1]*T01*ftT
    #     THSTMy[x]=-rho[1]*T01*ftT
    # nSve='ResV'+str(k)+str(100)+'.txt'
    # nSveS='ResS'+str(k)+str(100)+'.txt'
    # # Res=np.zeros([2,Nx])
    # # Res[0,:]=np.loadtxt(nSve, delimiter='\t')
    # # Res[1,:]=np.loadtxt(nSveS, delimiter='\t')
    # SolV,SolS=analytics.analytics_Cauchy_1IModK(Tmax,X,sourcef,modulation100,frontiere,mat,NOde,dt,Ntfig+1)
    # Res,Vfilm,Sfilm,CSVM=scheme1D.FD_cauchyII(U0,Ntfig,Nx,K,rhox,Celx,dt,dx,configuration,mat,frontiere,modulation100,pFilm,rCSVM="yes")
    # # np.savetxt(nSve, Res[0,:], fmt='%18e', delimiter='\t')
    # # np.savetxt(nSveS, Res[1,:], fmt='%18e', delimiter='\t')
    # errIm[0,k],errRIm[0,k]=np.sqrt(dx)*np.sqrt(errL2(Res[0,:], SolV))
    # errIm[1,k],errRIm[1,k]=np.sqrt(dx)*np.sqrt(errL2(Res[1,:], SolS))
    # ################################################
    # nSve='ResV'+str(k)+str(10)+'.txt'
    # nSveS='ResS'+str(k)+str(10)+'.txt'
    # # Res=np.zeros([2,Nx])
    # # Res[0,:]=np.loadtxt(nSve, delimiter='\t')
    # # Res[1,:]=np.loadtxt(nSveS, delimiter='\t')
    # SolV,SolS=analytics.analytics_Cauchy_1IModK(Tmax,X,sourcef,modulation10,frontiere,mat,NOde,dt,Ntfig+1)
    # Res,Vfilm,Sfilm,CSVM=scheme1D.FD_cauchyII(U0,Ntfig,Nx,K,rhox,Celx,dt,dx,configuration,mat,frontiere,modulation10,pFilm,rCSVM="yes")
    # errIm[2,k],errRIm[2,k]=np.sqrt(dx)*np.sqrt(errL2(Res[0,:], SolV))
    # errIm[3,k],errRIm[3,k]=np.sqrt(dx)*np.sqrt(errL2(Res[1,:], SolS))
    # np.savetxt(nSve, Res[0,:], fmt='%18e', delimiter='\t')
    # np.savetxt(nSveS, Res[1,:], fmt='%18e', delimiter='\t')
    # ################################################
    # nSve='ResV'+str(k)+str(0)+'.txt'
    # nSveS='ResS'+str(k)+str(0)+'.txt'
    Res=np.zeros([2,Nx])
    # Res[0,:]=np.loadtxt(nSve, delimiter='\t')
    # Res[1,:]=np.loadtxt(nSveS, delimiter='\t')
    SolV,SolS=analytics.analytics_Cauchy_1IModK(Tmax,X,sourcef,modulation0,frontiere,mat,NOde,dt,Ntfig+1)
    # Res,Vfilm,Sfilm,CSVM=scheme1D.FD_cauchyII(U0,Ntfig,Nx,K,rhox,Celx,dt,dx,configuration,mat,frontiere,modulation0,pFilm,rCSVM="yes")
    errIm[4,k],errRIm[2,k]=np.sqrt(dx)*np.sqrt(errL2(Res[0,:], SolV))
    errIm[5,k],errRIm[3,k]=np.sqrt(dx)*np.sqrt(errL2(Res[1,:], SolS))
    # # np.savetxt(nSve, Res[0,:], fmt='%18e', delimiter='\t')
    # # np.savetxt(nSveS, Res[1,:], fmt='%18e', delimiter='\t')
    # plt.figure()
    # plt.plot(X,SolV,'*')
    # plt.plot(X,Res[0,:])
    # plt.figure()
    # plt.plot(X,SolS,'*')
    # plt.plot(X,Res[1,:])
    
plt.figure()
plt.loglog(data[1,1:-1],errIm[0,1:-1],'*-',label=r"$f_m=100$ Hz")
plt.loglog(data[1,1:-1],errIm[2,1:-1],'o--',label=r"$f_m=10$ Hz")

###############################

with open("DataBruno/LW-S-Iter0.txt", "r") as fichier:
    lignes = fichier.readlines()
colonne1 = []
colonne2 = []
for ligne in lignes:
    mots = ligne.split()
    if len(mots) >= 2:
        try:
            valeur_colonne1 = float(mots[0])
            valeur_colonne2 = float(mots[1])
            colonne1.append(valeur_colonne1)
            colonne2.append(valeur_colonne2)
        except ValueError:
            pass
LWS0x = np.array(colonne1)
LWS0y = np.array(colonne2)        


with open("DataBruno/LW-V-Iter0.txt", "r") as fichier:
    lignes = fichier.readlines()
colonne1 = []
colonne2 = []
for ligne in lignes:
    mots = ligne.split()
    if len(mots) >= 2:
        try:
            valeur_colonne1 = float(mots[0])
            valeur_colonne2 = float(mots[1])
            colonne1.append(valeur_colonne1)
            colonne2.append(valeur_colonne2)
        except ValueError:
            pass
LWV0x = np.array(colonne1)
LWV0y = np.array(colonne2)

with open("DataBruno/LW-V-Iter300.txt", "r") as fichier:
    lignes = fichier.readlines()
colonne1 = []
colonne2 = []
for ligne in lignes:
    mots = ligne.split()
    if len(mots) >= 2:
        try:
            valeur_colonne1 = float(mots[0])
            valeur_colonne2 = float(mots[1])
            colonne1.append(valeur_colonne1)
            colonne2.append(valeur_colonne2)
        except ValueError:
            pass
LWV300x = np.array(colonne1)
LWV300y = np.array(colonne2)

plt.figure()
plt.plot(X,U0[0,:])
plt.plot(X,LWV0y.T,'--')
plt.xlabel("$x$ (m)")
plt.ylabel("$v$ (m/s)")
plt.title("LW : $t=0 s$")
plt.savefig('Figures/Comp_LW-V_0.eps', format='eps')

plt.figure()
plt.plot(X,U0[0,:]-LWV0y.T)
plt.xlabel("$x$ (m)")
plt.ylabel("$v-v_{Bruno}$ (m/s)")
plt.title("LW : $t=0 s$")
plt.savefig('Figures/Diff_LW-V_0.eps', format='eps')

plt.figure()
plt.plot(X,U0[1,:])
plt.plot(X,LWS0y.T,'--')
plt.xlabel("$x$ (m)")
plt.ylabel("$\sigma$ (Pa)")
plt.title("LW : $t=0 s$")
plt.savefig('Figures/Comp_LW-S_0.eps', format='eps')

plt.figure()
plt.plot(X,U0[1,:]-LWS0y.T)
plt.xlabel("$x$ (m)")
plt.ylabel("$\sigma-\sigma_{Bruno}$ (Pa)")
plt.title("LW : $t=0 s$")
plt.savefig('Figures/Diff_LW-S_0.eps', format='eps')

plt.figure()
plt.plot(X,Res[0,:])
plt.plot(X,LWV300y.T,'--')
plt.xlabel("$x$ (m)")
plt.ylabel("$v$ (m/s)")
plt.title("LW : $t=$"+tmax)
plt.savefig('Figures/Comp_LW-V_300.eps', format='eps')

plt.figure()
plt.plot(X,Res[0,:]-LWV300y.T)
plt.xlabel("$x$ (m)")
plt.ylabel("$v-v_{Bruno}$ (m/s)")
plt.title("LW : $t=$"+tmax)
plt.savefig('Figures/Diff_LW-V_300.eps', format='eps')


K=4
Res=scheme1D.FD_cauchy(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")

plt.figure()
plt.plot(X,U0[0,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$v$ (m/s)")
plt.title("ADER4 : $t=$0 s")
plt.savefig('Figures/ADER4-V_0.eps', format='eps')

plt.figure()
plt.plot(X,U0[1,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$\sigma$ (Pa)")
plt.title("ADER4 : $t=$0 s")
plt.savefig('Figures/ADER4-S_0.eps', format='eps')


plt.figure()
plt.plot(X,Res[0,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$v$ (m/s)")
plt.title("ADER4 : $t=$"+tmax)
plt.savefig('Figures/ADER4-V_300.eps', format='eps')

plt.figure()
plt.plot(X,Res[1,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$\sigma$ (Pa)")
plt.title("ADER4 : $t=$"+tmax)
plt.savefig('Figures/ADER4-S_300.eps', format='eps')

with open("DataBruno/ADER4-V-Iter300.txt", "r") as fichier:
    lignes = fichier.readlines()
colonne1 = []
colonne2 = []
for ligne in lignes:
    mots = ligne.split()
    if len(mots) >= 2:
        try:
            valeur_colonne1 = float(mots[0])
            valeur_colonne2 = float(mots[1])
            colonne1.append(valeur_colonne1)
            colonne2.append(valeur_colonne2)
        except ValueError:
            pass
A4V300x = np.array(colonne1)
A4V300y = np.array(colonne2)

plt.figure()
plt.plot(X,Res[0,:])
plt.plot(X,A4V300y.T,'--')
plt.xlabel("$x$ (m)")
plt.ylabel("$v$ (m/s)")
plt.title("ADER4 : $t=$"+tmax)
plt.savefig('Figures/Comp_ADER4-V_300.eps', format='eps')

plt.figure()
plt.plot(X,Res[0,:]-A4V300y.T)
plt.xlabel("$x$ (m)")
plt.ylabel("$v-v_{Bruno}$ (m/s)")
plt.title("ADER4 : $t=$"+tmax)
plt.savefig('Figures/Diff_ADER4-V_300.eps', format='eps')
