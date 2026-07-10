#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 18:01:38 2025

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
import funcPWE
   
vecteta0=[]
vecteta1=[]
 
for imod in range(1,2):
    configuration=inputReading.readfile('Demarrer.txt')
    mat=inputReading.readfile('Milieu.txt')
    frontiere=inputReading.readfile('Frontiere.txt')     
    sourcef=inputReading.readfile('Source'+str(imod)+'.txt')   
    modulation=inputReading.readfile('Modulation5.txt')
    ################ Parameters of the simulation   
    X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(configuration)
    fs,CFL,x_0=inputReading.source(sourcef)
    Nmat,rho,cm=inputReading.material(mat)
    alpha=inputReading.frontiere(frontiere)
    rhox,Celx=inputReading.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
    
    
        
    Tfin=0.5
    K=4
    pFilm=1
    dt=min(scheme1D.timeStep(CFL,dx,cm))
    Ntinit=2
    Uzn,Un,DeltaT,P1n,P2n,dP1n,dP2n=schemeH.temporalSchemePS(Ntinit,dt,configuration,sourcef,modulation,frontiere,mat,ordre0="yes")
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
    
    
    Uz,U,DeltaT,P1,P2,dP1,dP2=schemeH.temporalSchemePS(Ntfin+1,delt,configuration,sourcef,modulation,frontiere,mat,ordre0="yes")
    
    nSveUz='Uz_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.txt'
    nSveU='U_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.txt'
    nSveP1='P1_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.txt'
    nSveP2='P2_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.txt'
    nSvedP1='dP1_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.txt'
    nSvedP2='dP2_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.txt'
    nSvetp='tp_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.txt'
    nSveVfilm='Vfilm_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.txt'
    
    #U,DeltaT,P1,P2=schemeH.temporalSchemePS(Nt,configuration,sourcef,modulation,frontiere,mat)
    
    K=4
    pFilm=1
    
    testcfl=Celx[0]*delt/dx
    
    U0=scheme1D.initPointSourceProblem(configuration,frontiere,mat,sourcef)
    Res,Vfilm,Sfilm=scheme1D.FD_sourceVIItn(U0,Ntfin+1,Nx,K,rhox,Celx,delt,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)#:(U0,Nt,Nx,K,rhox,Celx,dt,dx,"no")
    #Resp1,Vfilmp1,Sfilmp1=scheme1D.FD_sourceVII(U0,Nt2,Nx,K,rhox,Celx,dt,dx,x_0,sce,configuration,mat,frontiere,modulation,pFilm)
    Ufilm=ffg.VtoU(Vfilm, delt)
    
    
    
    np.savetxt(nSveU, U, fmt='%18e', delimiter='\t')
    np.savetxt(nSveUz, Uz, fmt='%18e', delimiter='\t')
    np.savetxt(nSveVfilm, Vfilm, fmt='%18e', delimiter='\t')
    np.savetxt(nSvetp, tp, fmt='%18e', delimiter='\t')
    np.savetxt(nSveP1, P1, fmt='%18e', delimiter='\t')
    np.savetxt(nSveP2, P2, fmt='%18e', delimiter='\t')
    np.savetxt(nSvedP1, dP1, fmt='%18e', delimiter='\t')
    np.savetxt(nSvedP2, dP2, fmt='%18e', delimiter='\t')
    
    #Ufilmp1=ffg.VtoU(Vfilmp1, dt)
    V=ffg.UtoV(U,delt)
    Vz=ffg.UtoV(Uz,delt)
    Ufilm=ffg.VtoU(Vfilm, delt)
    
    
    Tfig=0.2
    Ntfig=int(np.ceil(Tfig/delt))
    U2C=schemeH.postCorrector(Ntfig-1,dx,Nx,U,P1,P2)
    V2C=schemeH.postCorrectorV(Ntfig-1,dx,Nx,U,V,P1,P2,dP1,dP2)
    
    
    nfig=4
    
    psmax=np.where(Ufilm[Ntfig-1,:Nx//2]==max(Ufilm[Ntfig-1,:Nx//2]))
    xzoommin=X[psmax[0][0]-200]
    xzoommax=X[psmax[0][0]+600]
    zoomlim=[xzoommin,xzoommax]
    
    plt.figure()
    plt.axvspan(xzoommin, xzoommax, facecolor='g', alpha=0.3)
    plt.scatter(X[::1],Ufilm[Ntfig,::1],marker=r'$\bigcirc$',label=r'$U_h$',color='blue',s=40)
    #plt.plot(X[::1],Ufilm[Ntfig,::1],linestyle='None',marker='o',color='white',markersize=3)
    plt.plot(X,Uz[Ntfig-1,:],'grey',linewidth=3,label=r'$U_0$')
    plt.plot(X,U2C,color='red',linestyle='None',marker='o',markersize=3,label=r'$U^{(2)}$')
    plt.xlabel(r'$X$ (m)',fontsize=16)
    plt.ylabel(r'$U$ (m)',fontsize=16)
    plt.title("$T=$"+str(round(tp[Ntfig], 3))+" s",fontsize=16)
    plt.legend(fontsize=14)#loc="lower right",)
    plt.tight_layout()
    namefigU='U_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.pdf'
    plt.savefig(namefigU,format='pdf')
    
    plt.figure()
    plt.axvspan(xzoommin, xzoommax, facecolor='g', alpha=0.3)
    plt.scatter(X[::nfig],Vfilm[Ntfig-1,::nfig],marker=r'$\bigcirc$',label=r'$V_h$',color='blue',s=40)
    #plt.plot(X[::1],Vfilm[Ntfin-1,::1],marker='o',color='white',linestyle=None,markersize=3)
    plt.plot(X,Vz[Ntfig-1,:],'grey',linewidth=3,label=r'$V_0$')
    plt.plot(X,V2C,color='red',linestyle='None',marker='o',markersize=3,label=r'$V^{(2)}$')
    plt.xlabel(r'$X$ (m)',fontsize=16)
    plt.ylabel(r'$V$ (m/s)',fontsize=16)
    plt.title("$T=$"+str(round(tp[Ntfig], 3))+" s",fontsize=16)
    plt.legend(fontsize=14)#loc="lower right",)
    plt.tight_layout()
    namefigV='V_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.pdf'
    plt.savefig(namefigV,format='pdf')
    
    
    
    
    
    plt.figure()
    for i in range(alpha.size):
        plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
    plt.scatter(X[::1],Ufilm[Ntfig,::1],marker=r'$\bigcirc$',label=r'$U_h$',color='blue',s=40)
    #plt.plot(X[::1],Ufilm[Ntfig,::1],linestyle='None',marker='o',color='white',markersize=3)
    plt.plot(X,Uz[Ntfig-1,:],'grey',linewidth=3,label=r'$U_0$')
    plt.plot(X,U2C,color='red',linestyle='None',marker='o',markersize=3,label=r'$U^{(2)}$')
    plt.xlabel(r'$X$ (m)',fontsize=16)
    plt.ylabel(r'$U$ (m)',fontsize=16)
    plt.xlim(zoomlim)
    plt.title("$T=$"+str(round(tp[Ntfig], 3))+" s",fontsize=16)
    plt.legend(fontsize=14)#loc="lower right",)
    plt.tight_layout()
    namefigUz='Uz_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.pdf'
    plt.savefig(namefigUz,format='pdf')
    
    plt.figure()
    for i in range(alpha.size):
        plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
    plt.scatter(X[::nfig],Vfilm[Ntfig-1,::nfig],marker=r'$\bigcirc$',label=r'$V_h$',color='blue',s=40)
    #plt.plot(X[::1],Vfilm[Ntfin-1,::1],marker='o',color='white',linestyle=None,markersize=3)
    plt.plot(X,Vz[Ntfig-1,:],'grey',linewidth=3,label=r'$V_0$')
    plt.plot(X,V2C,color='red',linestyle='None',marker='o',markersize=3,label=r'$V^{(2)}$')
    plt.xlabel(r'$X$ (m)',fontsize=16)
    plt.ylabel(r'$V$ (m/s)',fontsize=16)
    plt.xlim(zoomlim)
    plt.title("$T=$"+str(round(tp[Ntfig], 3))+" s",fontsize=16)
    plt.legend(fontsize=14)#loc="lower right",)
    plt.tight_layout()
    namefigVz='Vz_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.pdf'
    plt.savefig(namefigVz,format='pdf')
    # ffg.dispTK(Ufilm,tp,dx)
    
    
    plt.figure()
    plt.plot(tp[:-1],Ufilm[:,int(Nx//2)+2],linewidth=3,label=r'$U$')
    plt.plot(tp[:],Uz[:,int(Nx//2)+2],':',linewidth=3,label=r'$U_0$')
    # plt.plot(tp,U[:,int(Nx//2)+2],'--',linewidth=3,label=r'$U_2$')
    plt.xlabel(r'$T$ (s)',fontsize=16)
    plt.ylabel(r'$U$ (m)',fontsize=16)
    plt.legend(fontsize=14)
    plt.tight_layout()
    namefigT='Utemp_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.pdf'
    plt.savefig(namefigT,format='pdf')
    
    # en
    vecteta0.append(eta0*2*np.pi)
    vecteta1.append(eta1*2*np.pi)
    
    # plt.figure()
    # plt.plot(tp[:-1],Efilm[:,int(Nx//2)+2],linewidth=3,label=r'$U$')
    # plt.plot(tp[:],Uz[:,int(Nx//2)+2],':',linewidth=3,label=r'$U_0$')
    # # plt.plot(tp,U[:,int(Nx//2)+2],'--',linewidth=3,label=r'$U_2$')
    # plt.xlabel(r'$T$ (s)',fontsize=16)
    # plt.ylabel(r'$U$ (m)',fontsize=16)
    # plt.legend(fontsize=14)
    # plt.tight_layout()
    # namefigT='Etemp_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.pdf'
    # plt.savefig(namefigT,format='pdf')
    
    Tfig=0.2
    Ntfig=int(np.ceil(Tfig/delt))
    
    fctModE="np.sin"
    fctModR="np.sin"
    NE=6
    if tau!=0:
        omegavect,VP=funcPWE.fillMatrix(tau,Eef0,rhoef0,EpsE,EpsR,fctModE,fctModR,NE)
        kmax=max(np.real(VP[NE-1,:])/2/np.pi)
        ffg.dispFK_tperiodic(Ufilm[:Ntfig,:],delt,dx,fm,kmax)
        for i in range(NE):
            plt.plot((VP[i,:])/2/np.pi,omegavect/2/np.pi,'--r')
        namefigDispMic='DispMic_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.pdf'
        plt.savefig(namefigDispMic,format='pdf')
        ffg.dispFK_tperiodic(Uz[:Ntfig,:],delt,dx,fm,kmax)
        for i in range(NE):
            plt.plot((VP[i,:])/2/np.pi,omegavect/2/np.pi,'--r')
        namefigDispU0='DispU0_K-1e9_fc'+str(fs)+'fm'+str(fm)+'.pdf'
        plt.savefig(namefigDispU0,format='pdf')
    