#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 13:54:55 2025

@author: michael
"""


import os
#os.chdir("/home/darche/Documents/ModulatedMetamaterials/WAPRO_TD/1-Interface_Imerged_PS_modulated2/")
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/Homog_Ordre0Ordre2/")

import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

import sys
sys.path.append("../WAPOD/")

import inputReading
import scheme1D_TDold as scheme1D
import scheme1D_TD_Homo as schemeH
import figures as ffg
import source
    
for imod in range(2,3):
    configuration=inputReading.readfile('Demarrer.txt')
    mat=inputReading.readfile('Milieu.txt')
    frontiere=inputReading.readfile('Frontiere.txt')     
    sourcef=inputReading.readfile('Source.txt')   
    modulation=inputReading.readfile('Modulation'+str(imod)+'.txt')
    modulation=inputReading.readfile('Modulation.txt')
    ################ Parameters of the simulation   
    X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(configuration)
    fs,CFL,x_0=inputReading.source(sourcef)
    Nmat,rho,cm=inputReading.material(mat)
    alpha=inputReading.frontiere(frontiere)
    rhox,Celx=inputReading.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
    
    
        
    Tfin=0.1
    K=4
    pFilm=1
    dt=min(scheme1D.timeStep(CFL,dx,cm))
    Ntinit=2
    Uzn,DeltaT=schemeH.temporalSchemePS0(Ntinit,dt,configuration,sourcef,modulation,frontiere,mat)
    Nt2=int(np.ceil(Tfin/dt))
    
    delt=min(dt,DeltaT)
    Ntfin=int(np.ceil(Tfin/delt))
    tp=np.linspace(0,Tfin,Ntfin+2)
    testcfl=Celx[0]*delt/dx
    sce=np.zeros(Ntfin+2)
    for i in range(Ntfin+2):
        sce[i]=source.choice_timefct(sourcef,tp[i])
        
    
    h=alpha[1]-alpha[0]
    alpha0=alpha[0]
    
    
    fm=float(modulation["Freq"])
    cel=float(modulation["ModC"])
    if cel !=0:
        Nip=int(cel/h/fm)
    else:
        Nip=1
    
    K0=float(frontiere["Kn_0"])
    M0=float(frontiere["Mn_0"])
    DK=float(modulation["DeltaF"])
    DM=float(modulation["DeltaM"])
    
    
    K1=K0/(1+DK)
    K2=K0/(1-DK)
    
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
    
    if Nip%2==1:
        rhoef1=rhoet+M0/h+DM*M0/h/Nip
        Eef1=1/(1/Eet+1/K0/h+DK/K0/h/Nip)
        rhoef2=rhoet+M0/h-DM*M0/h/Nip
        Eef2=1/(1/Eet+1/K0/h-DK/K0/h/Nip)
    else:
        rhoef1=rhoet+M0/h
        Eef1=1/(1/Eet+1/K0/h)
        rhoef2=rhoet+M0/h
        Eef2=1/(1/Eet+1/K0/h)
        
    cmef1=np.sqrt(Eef1/rhoef1)
    cmef2=np.sqrt(Eef2/rhoef2)
    cmefmin=min(cmef1,cmef2)
    
    ratio=float(modulation["Ratio"])
    Eef=1/(ratio/Eef1+(1-ratio)/Eef2)
    rhoef=ratio*rhoef1+(1-ratio)*rhoef2
    
    
    c_r = cmef1*ratio+cmef2*(1-ratio)
    c_0 = np.sqrt(Eef/rhoef)
    
    if fm!=0:
        tau=1/fm
    else:
        tau=0
    
    eta0=fs/c_0*h*Nip
    eta1=(fs+fm)/c_0*h*Nip
    
    rhoef0=rhoet+M0/h
    EpsR=DM*M0/h/Nip/rhoef0
    Eef0=Eet*K0*h/(K0*h+Eet)
    EpsE=DK*Eef0/K0/h/Nip
    
    
    U00,DeltaT=schemeH.temporalSchemePS0(Ntfin+1,delt,configuration,sourcef,modulation,frontiere,mat)
    
    nSveUz0='Uz0_KM_fc'+str(fs)+'fm'+str(fm)+'.txt'
    nSveU0='U0_KM_fc'+str(fs)+'fm'+str(fm)+'.txt'
    nSvetp0='tp0_KM_fc'+str(fs)+'fm'+str(fm)+'.txt'
    nSveVfilm0='Vfilm0_KM_fc'+str(fs)+'fm'+str(fm)+'.txt'
    
    #U,DeltaT,P1,P2=schemeH.temporalSchemePS(Nt,configuration,sourcef,modulation,frontiere,mat)
    
    K=4
    pFilm=1
    
    testcfl=Celx[0]*delt/dx
    testcfl2=c_0*delt/dx
    
    
    U0=scheme1D.initPointSourceProblem(configuration,frontiere,mat,sourcef)
    Res,Vfilm,Sfilm=scheme1D.FD_sourceVII(U0,Ntfin+1,Nx,K,rhox,Celx,delt,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
    #Resp1,Vfilmp1,Sfilmp1=scheme1D.FD_sourceVII(U0,Nt2,Nx,K,rhox,Celx,dt,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)
    Ufilm=ffg.VtoU(Vfilm, delt)
    
    
    
    np.savetxt(nSveU0, U00, fmt='%18e', delimiter='\t')
    np.savetxt(nSveVfilm0, Vfilm, fmt='%18e', delimiter='\t')
    np.savetxt(nSvetp0, tp, fmt='%18e', delimiter='\t')
    
    #Ufilmp1=ffg.VtoU(Vfilmp1, dt)
    V0=ffg.UtoV(U00,delt)
    
    
    
    Tfig=0.1
    Ntfig=int(np.ceil(Tfig/delt))
    
    
    nfig=4
    
    
    plt.figure()
    plt.scatter(X[::1],Ufilm[Ntfig,::1],marker=r'$\bigcirc$',label=r'$U_h$',color='blue',s=40)
    #plt.plot(X[::1],Ufilm[Ntfig,::1],linestyle='None',marker='o',color='white',markersize=3)
    plt.plot(X,U00[Ntfig-1,:],'grey',linewidth=3,label=r'$U_0$')
    plt.xlabel(r'$X$ (m)',fontsize=16)
    plt.ylabel(r'$U$ (m/s)',fontsize=16)
    plt.title("$T=$"+str(round(tp[Ntfig], 3))+" s",fontsize=16)
    plt.legend(fontsize=14)#loc="lower right",)
    plt.tight_layout()
    namefigU='U0_KM_fc'+str(fs)+'fm'+str(fm)+'.pdf'
    plt.savefig(namefigU,format='pdf')
    
    plt.figure()
    plt.scatter(X[::nfig],Vfilm[Ntfig-1,::nfig],marker=r'$\bigcirc$',label=r'$V_h$',color='blue',s=40)
    #plt.plot(X[::1],Vfilm[Ntfin-1,::1],marker='o',color='white',linestyle=None,markersize=3)
    plt.plot(X,V0[Ntfig-1,:],'grey',linewidth=3,label=r'$V_0$')
    plt.xlabel(r'$X$ (m)',fontsize=16)
    plt.ylabel(r'$V$ (m/s)',fontsize=16)
    plt.title("$T=$"+str(round(tp[Ntfig], 3))+" s",fontsize=16)
    plt.legend(fontsize=14)#loc="lower right",)
    plt.tight_layout()
    namefigV='V0_KM_fc'+str(fs)+'fm'+str(fm)+'.pdf'
    plt.savefig(namefigV,format='pdf')
    
    
    psmax=np.where(Ufilm[Ntfig-1,:Nx//2]==max(Ufilm[Ntfig-1,:Nx//2]))
    zoomlim=[X[psmax[0][0]-200],X[psmax[0][0]+600]]
    
    
    plt.figure()
    for i in range(alpha.size):
        plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
    plt.scatter(X[::1],Ufilm[Ntfig,::1],marker=r'$\bigcirc$',label=r'$U_h$',color='blue',s=40)
    #plt.plot(X[::1],Ufilm[Ntfig,::1],linestyle='None',marker='o',color='white',markersize=3)
    plt.plot(X,U00[Ntfig-1,:],'grey',linewidth=3,label=r'$U_0$')
    plt.xlabel(r'$X$ (m)',fontsize=16)
    plt.ylabel(r'$U$ (m/s)',fontsize=16)
    plt.xlim(zoomlim)
    plt.title("$T=$"+str(round(tp[Ntfig], 3))+" s",fontsize=16)
    plt.legend(fontsize=14)#loc="lower right",)
    plt.tight_layout()
    namefigUz='Uz0_KM_fc'+str(fs)+'fm'+str(fm)+'.pdf'
    plt.savefig(namefigUz,format='pdf')
    
    plt.figure()
    for i in range(alpha.size):
        plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
    plt.scatter(X[::nfig],Vfilm[Ntfig-1,::nfig],marker=r'$\bigcirc$',label=r'$V_h$',color='blue',s=40)
    #plt.plot(X[::1],Vfilm[Ntfin-1,::1],marker='o',color='white',linestyle=None,markersize=3)
    plt.plot(X,V0[Ntfig-1,:],'grey',linewidth=3,label=r'$V_0$')
    plt.xlabel(r'$X$ (m)',fontsize=16)
    plt.ylabel(r'$V$ (m/s)',fontsize=16)
    plt.xlim(zoomlim)
    plt.title("$T=$"+str(round(tp[Ntfig], 3))+" s",fontsize=16)
    plt.legend(fontsize=14)#loc="lower right",)
    plt.tight_layout()
    namefigVz='Vz0_KM_fc'+str(fs)+'fm'+str(fm)+'.pdf'
    plt.savefig(namefigVz,format='pdf')
    # ffg.dispTK(Ufilm,tp,dx)
    
    
    plt.figure()
    plt.plot(tp[:-1],Ufilm[:,int(Nx//2)+2],linewidth=3)
    plt.plot(tp,U00[:,int(Nx//2)+2],'--',linewidth=3)
    plt.xlabel(r'$T$ (s)',fontsize=16)
    plt.ylabel(r'$U$ (m)',fontsize=16)
    plt.tight_layout()
    namefigT='Utemp0_KM_fc'+str(fs)+'fm'+str(fm)+'.pdf'
    plt.savefig(namefigT,format='pdf')
    
    
kmax=4*fs/min(cmef1,cmef2)
ffg.dispFK_tperiodic(Ufilm,delt,dx,100,kmax)
ffg.dispFK_tperiodic(U00,delt,dx,100,kmax)

