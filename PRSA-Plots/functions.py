#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 16:53:54 2026

@author: michael
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg as ssl

##################
## Input reading
##################

def readfile(filename):
    with open(filename,'r') as file:
        namedict={}
        if filename[0:6]=="Milieu":
            line_nb_mat = file.readline().strip()
            part_nb_milieu = line_nb_mat.strip().split(':')
            namenbmat=part_nb_milieu[0]
            nb_mat=part_nb_milieu[1]
            namedict[namenbmat] = nb_mat
            k=0
            for line in file:
                parts = line.strip().split(':')
                nummed="milieu "+str(k)
                if parts[0].strip()==nummed:
                    if len(parts) == 2:
                        variable_name = parts[0].strip()+"_"+str(k)
                        value = parts[1].strip()
                        namedict[variable_name] = value
                    k=k+1
                else:
                    if len(parts) == 2:
                        variable_name = parts[0].strip()+"_"+str(k-1)
                        value = parts[1].strip()
                        namedict[variable_name] = value
        elif filename[0:9]=="Frontiere":
            line_nb_inter = file.readline().strip()
            part_nb_inter = line_nb_inter.strip().split(':')
            namenbinter=part_nb_inter[0]
            nb_inter=part_nb_inter[1]
            namedict[namenbinter] = nb_inter
            k=0
            for line in file:
                parts = line.strip().split(':')
                numinter="frontiere "+str(k+1)
                if parts[0].strip()==numinter:
                    if len(parts) == 2:
                        variable_name = parts[0].strip()+"_"+str(k)
                        value = parts[1].strip()
                        namedict[variable_name] = value
                    k=k+1
                else:
                    if len(parts) == 2:
                        variable_name = parts[0].strip()+"_"+str(k-1)
                        value = parts[1].strip()
                        namedict[variable_name] = value
        else :
            for line in file:
                parts = line.strip().split(':')
                if len(parts) == 2:
                    variable_name = parts[0].strip()
                    value = parts[1].strip()
                    namedict[variable_name] = value
    return namedict
        
# Geometry
def geometric(configuration):
    xmin=float(configuration["Xinf"])
    xmax=float(configuration["Xsup"])
    Nx=int(configuration["Nx"])+1
    ESIM=int(configuration["Esim"])
    Npt=int(configuration["Npt"])
    dx=(xmax-xmin)/(Nx-1)
    X=np.linspace(xmin,xmax,Nx)
    return X,xmin,xmax,Nx,dx,ESIM,Npt

def material(mat):
    Nmat=int(mat["nombre de milieux "])
    rho=np.zeros(Nmat)
    cm=np.zeros(Nmat)
    # E=np.zeros(Nmat)
    for i in range(Nmat):
        nameRho="Rho"+'_'+str(i)
        rho[i]=float(mat[nameRho])
        nameCm="Cel"+'_'+str(i)
        cm[i]=float(mat[nameCm])
        # E[i]=cm[i]**2*rho[i]
        #A=np.array([[0,-1/rho],[-E,0]])
    return Nmat,rho,cm#,E,A

def frontiere(front):
    Nf=int(front["nombre de frontieres "])
    alphaF=np.zeros(Nf+1)
    for i in range(Nf):
        nameFront="Alpha_"+str(i)
        alphaF[i]=float(front[nameFront])
    return alphaF


def affectMaterials(alphaF,rho,cm,xmin,xmax,Nx,dx):
    minf=xmin
    alphaF[-1]=xmax
    msup=alphaF[0]
    rhox=np.zeros(Nx)
    Celx=np.zeros(Nx)
    med=0
    for mx in range(Nx):
        if mx*dx+xmin<minf:
            print("erreur geometrique")
        elif mx*dx+xmin<= msup:
            rhox[mx]=rho[med]
            Celx[mx]=cm[med]
        else:
            med=med+1
            rhox[mx]=rho[med]
            Celx[mx]=cm[med]
            minf=msup
            msup=alphaF[med]
    return rhox,Celx

def affectMaterialsH(rho,cm,xmin,xmax,Nx,dx):
    rhox=np.zeros(Nx)
    Celx=np.zeros(Nx)
    for mx in range(Nx):
        rhox[mx]=rho
        Celx[mx]=cm
    return rhox,Celx

def comE(rho,c):
    E=c**2*rho
    return E

def comA(rho,c):
    E=comE(rho,c)
    A=np.array([[0,-1/rho],[-E,0]])
    return A     

def comR(rho,c):
    E=comE(rho,c)
    R=np.array([[0,-E],[-1/rho,0]])
    return R

def source(sourcef):
    fs=float(sourcef["Frequence"])
    CFL=float(sourcef["CFL"])
    if sourcef["Forme"]=="point source":
        x0=float(sourcef["Xsource"])
        return fs,CFL,x0
    elif sourcef["Forme"]=="Cauchy":
        tshift=float(sourcef["Tshift"]) 
        return fs,CFL,tshift



def fonctionModMP(t,frontiere,i): 
    contacti="Contact_"+str(i)
    if frontiere[contacti]=="masse-ressort module independant":
        Delta=float(frontiere["DeltaM_"+str(i)])
        if frontiere["Synchro_"+str(i)]=="Yes":
            Omega=2*np.pi*float(frontiere["Freq_"+str(i)])
            if frontiere["Dephasage0_"+str(i)]=="Yes":
                phi0=float(frontiere["phiM0_"+str(i)])
            elif frontiere["Dephasage0_"+str(i)]== "x2Pi":
                phi0=2*np.pi*float(frontiere["phiM0_"+str(i)])
            else:
                phi0=0    
        else:
            Omega=2*np.pi*float(frontiere["FreqM_"+str(i)])
            if frontiere["Dephasage0_"+str(i)]=="Yes":
                phi0=float(frontiere["phiM0_"+str(i)])
            elif frontiere["Dephasage0_"+str(i)]== "x2Pi":
                phi0=2*np.pi*float(frontiere["phiM0_"+str(i)])
            else:
                phi0=0
        if frontiere["TypeMod_"+str(i)]=="Sinus":
            fct=Delta*np.sin(Omega*t+phi0)
    elif frontiere[contacti]=="parfait":
        fct=0
    else:
        fct=0
    return fct

def fonctionModKP(t,frontiere,i): 
    contacti="Contact_"+str(i)
    if frontiere[contacti]=="masse-ressort module independant":
        Delta=float(frontiere["DeltaF_"+str(i)])
        if frontiere["Synchro_"+str(i)]=="Yes":
            Omega=2*np.pi*float(frontiere["Freq_"+str(i)])
            if frontiere["Dephasage0_"+str(i)]=="Yes":
                phi0=float(frontiere["phiF0_"+str(i)])
            elif frontiere["Dephasage0_"+str(i)]== "x2Pi":
                phi0=2*np.pi*float(frontiere["phiF0_"+str(i)])
            else:
                phi0=0    
        else:
            Omega=2*np.pi*float(frontiere["FreqF_"+str(i)])
            if frontiere["Dephasage0_"+str(i)]=="Yes":
                phi0=float(frontiere["phiF0_"+str(i)])
            elif frontiere["Dephasage0_"+str(i)]== "x2Pi":
                phi0=2*np.pi*float(frontiere["phiF0_"+str(i)])
            else:
                phi0=0
        if frontiere["TypeMod_"+str(i)]=="Sinus":
            fct=Delta*np.sin(Omega*t+phi0)
    elif frontiere[contacti]=="parfait":
        fct=0
    else:
        fct=0
    return fct

def propEff(rho,c,C,M,h):
    rhoeff=rho+M/h
    unsE=1/(rho*c**2)+C/h
    Eeff=1/unsE
    return rhoeff,Eeff

##################
## Data processing
##################
# Frequency-wavenumber FFT
def specfkfreqneg(U,dx,dt,st=4,sx=4):
    Nt=np.size(U[::st,0])
    Nx=np.size(U[0,::sx])
    # Zero-padding
    zpt=2**14
    zps=2**14
    Uzp=np.zeros([zpt,zps])
    Uzp[:Nt,int(zps/2-Nx/2):int(zps/2+Nx/2)]=U[::st,::sx]
    #Fourier in time domain 
    FFTt = np.fft.fft(Uzp[::], axis=0)
    #Fourier in space
    FFTx = np.fft.fft(np.fft.fftshift(FFTt, axes=0), axis=1)
    FFK = np.absolute(FFTx)
    # Get the frequency and K vectors
    frqv = np.fft.fftfreq(zpt, st*dt)
    wavv = np.fft.fftfreq(zps, sx*dx)
    # Centering
    wavv = np.fft.fftshift(wavv)
    frqv = np.fft.fftshift(frqv)
    return frqv,wavv,FFK
    
# Integration of velocity fields
def VtoU(datas,dt,u0=0):
    Nt=np.size(datas[:,0])
    Nx=np.size(datas[0,:])
    Ufilm=np.zeros([Nt,Nx])
    for i in range(0,Nt-1):
        Ufilm[i+1,:]=Ufilm[i,:]+dt/2*(datas[i+1,:]+datas[i,:])
    return Ufilm
    
### Derivative of displacement fields
# Time
def UtoV(datas,dt):
    Nt=np.size(datas[:,0])
    Nx=np.size(datas[0,:])
    Vfilm=np.zeros([Nt,Nx])
    for i in range(1,Nt-1):
        Vfilm[i,:]=(datas[i+1,:]-datas[i-1,:])/2/dt
    return Vfilm    
# Space 
def UtoS(datas,dx):
    Nt=np.size(datas[:,0])
    Nx=np.size(datas[0,:])
    Sfilm=np.zeros([Nt,Nx])
    for i in range(1,Nx-1):
        Sfilm[:,i]=(datas[:,i+1]-datas[:,i-1])/2/dx
    return Sfilm
    
### Energies
# Microstructured media
def energyMicro(Vfilm,Sfilm,CSVM,Nt,dx,rhox,Celx,nfront,modulC,modulM,frontiere):
    energyV=np.zeros(Nt)
    energyI=np.zeros(Nt)
    for ti in range(Nt):
        energyV[ti]=dx/2*np.sum(rhox*Vfilm[ti]*Vfilm[ti])+dx/2*np.sum(1/rhox/Celx/Celx*Sfilm[ti]*Sfilm[ti])
        for interi in range(int(nfront)):
            if frontiere["Contact_"+str(interi)]!="parfait":
                energyI[ti]=energyI[ti]=energyI[ti]+1/2*modulM[ti]*CSVM[ti,2+4*interi]*CSVM[ti,2+4*interi]+1/2*modulC[ti]*CSVM[ti,3+4*interi]*CSVM[ti,3+4*interi] 
    energymicro=energyV+energyI
    return energymicro

# Homogenised media
def energyHomogenisedU(V,S,Nt,dx,dt,rhoeff,Eeff):
    energyH=np.zeros(Nt)
    for ti in range(Nt):
        energyH[ti]=dx/2*np.sum(rhoeff[ti]*V[ti]*V[ti])+dx/2*np.sum(Eeff[ti]*S[ti]*S[ti])
    return energyH

def energyHomogenisedVS(V,S,Nt,dx,dt,rhoeff,Eeff):
    energyH=np.zeros(Nt)
    for ti in range(Nt):
        energyH[ti]=dx/2*np.sum(rhoeff[ti]*V[ti]*V[ti])+dx/2*np.sum(1/Eeff[ti]*S[ti]*S[ti])
    return energyH

##################
## Correctors for order-2 homogenised fields
##################
# On displacement field
def postCorrector(Nt,dx,Nx,Utot,P1tot,P2tot):
    U=Utot[Nt,:]
    P1=P1tot[Nt,:]
    P2=P2tot[Nt,:]
    dxU=np.zeros(Nx)
    dxxU=np.zeros(Nx)
    for i in range(1,Nx-1):
        dxU[i]=(U[i+1]-U[i-1])/2/dx
    for i in range(2,Nx-2):
        dxxU[i]=(dxU[i+1]-dxU[i-1])/2/dx
    U2=U+P1*dxU+P2*dxxU
    return U2
    
# On velocity field   
def postCorrectorV(Nt,dx,Nx,Utot,Vtot,P1tot,P2tot,dP1tot,dP2tot):
    U=Utot[Nt,:]
    V=Vtot[Nt,:]
    P1=P1tot[Nt,:]
    P2=P2tot[Nt,:]
    dP1=dP1tot[Nt,:]
    dP2=dP2tot[Nt,:]
    dxU=np.zeros(Nx)
    dxxU=np.zeros(Nx)
    dxV=np.zeros(Nx)
    dxxV=np.zeros(Nx)
    for i in range(1,Nx-1):
        dxU[i]=(U[i+1]-U[i-1])/2/dx
        dxV[i]=(V[i+1]-V[i-1])/2/dx
    for i in range(2,Nx-2):
        dxxU[i]=(dxU[i+1]-dxU[i-1])/2/dx
        dxxV[i]=(dxV[i+1]-dxV[i-1])/2/dx
    V2=V+P1*dxV+P2*dxxV+dP1*dxU+dP2*dxxU
    return V2

##################
## Plane Wave Expansion
##################
# Heavyside
def fctheavyside(t,tau,phi):
    Tphi=tau*phi
    res=t.copy()
    for i in range(t.shape[0]):
        if t[i]<Tphi:
            res[i]=+1
        else:
            res[i]=-1
    return res    
### Some time-modulation of the material properties 
# Mass density
def rhot(t,tau,rho,Delta,fctM):
    if fctM=="np.sin":
        return rho*(1+Delta*np.sin((2*np.pi*t/tau)))
    if fctM=="np.cos":
        return rho*(1+Delta*np.cos((2*np.pi*t/tau)))
    elif fctM=="np.sign":
        rhot0=np.sin(2*np.pi*t/tau)
        return rho*(1+Delta*np.sign(rhot0))
    elif fctM=="heavy":
        return rho*(1+Delta*fctheavyside(t,tau))
    else:
        return 0

# Young modulus 
def Et(t,tau,E,delta,fctM):
    if fctM=="np.sin":
        return 1/(1/E*(1+delta*np.sin((2*np.pi*t/tau))))
    if fctM=="np.cos":
        return 1/(1/E*(1+delta*np.cos((2*np.pi*t/tau))))
    elif fctM=="np.sign":
        e0=np.sin((2*np.pi*t/tau))
        return 1/(1/E*(1+delta*np.sign(e0)))
    elif fctM=="heavy":
        return 1/(1/E*(1+delta*fctheavyside(t,tau)))
    else:
        return 0

# Fourier decomposition
def SF(fct,n,tau,val,delta,fctMod,M):
    t=np.linspace(0,tau,M)
    dt=t[1]-t[0]
    an=dt/tau*np.sum((fct(t[1:],tau,val,delta,fctMod)+fct(t[:-1],tau,val,delta,fctMod))/2*np.exp(-2j*n*np.pi/tau*(t[1:]+t[:-1])/2))   
    return an
    
# Dispersion curves for time-modulated
def fillMatrix(tau,valE,valR,deltaE,deltaR,fctModE,fctModR,NE,Nom=50,N=16,Ni=500):
    MatrixA=np.zeros([2*N+1,2*N+1],dtype=complex)
    MatrixB=np.zeros([2*N+1,2*N+1],dtype=complex)
    VP=np.zeros([NE,Nom+1],dtype=complex)
    om=np.linspace(-np.pi/tau,np.pi/tau,Nom+1)
    for omeg in range(Nom+1):
        for n in range(-N,N+1):
            nn=n+N
            for p in range(-N,N+1):
                pp=p+N
                anE=SF(Et,p-n,tau,valE,deltaE,fctModE,Ni)
                anR=SF(rhot,p-n,tau,valR,deltaR,fctModR,Ni)
                MatrixA[pp,nn]=(2*np.pi*n/tau-om[omeg])*(2*np.pi*p/tau-om[omeg])*anR
                MatrixB[pp,nn]=anE
        EigV=ssl.eigsh(MatrixA,k=NE,M=MatrixB,which='SM',return_eigenvectors=False)
        VP[:,omeg]=np.sqrt(EigV)
    VP=np.nan_to_num(VP) 
    VP.sort(axis=0)
    return om,VP

##################
## Make plots
##################
# Plot configuration
def format_fig(x,y,title=None,legend=None,sizef=16):
    plt.xlabel(x,fontsize=sizef)
    plt.ylabel(y,fontsize=sizef)
    if title!=None:
        plt.title(title,fontsize=sizef)
    if title!=None:
        plt.legend(fontsize=sizef)
    plt.tight_layout()
    return

# Space-time representation
def sismo(x,t,data1,data2,numLines,pFilm=1,Nt=None):
    if Nt==None:
        Nt=np.size(data1[0,:])
    maxDatas=np.max([np.max(data1[:Nt:int(Nt/numLines),:]),-np.min(data1[:Nt:int(Nt/numLines),:])])
    datas1=data1[:Nt,:]
    datas2=data2[:Nt,:]
    plt.figure()
    for j in range(numLines):
        plt.plot(x,(datas2[j*int(Nt/numLines),:]+2*j*1.02*maxDatas)*pFilm*int(Nt/numLines),'#1f77b4',linewidth=2)
        plt.plot(x,(datas1[j*int(Nt/numLines),:]+2*j*1.02*maxDatas)*pFilm*int(Nt/numLines),':r',linewidth=3)
    plt.yticks([i * 2 * 1.02 * maxDatas * pFilm * int(Nt/numLines) for i in range(0,numLines,2)],np.round(t[:Nt:2*pFilm*int(Nt/numLines)], 3))
    format_fig(r"$X$ (m)",r"$T$ (s)")
    
# Dispersion diagram
def dispFK_tperiodic2(datas,dt,dx,fm,kmax):
    frqv,wavv,FFK=specfkfreqneg(datas,dx,dt)
    numk=np.size(wavv)
    numf=np.size(frqv)
    df=frqv[2]-frqv[1]
    Nf=int(fm/df)+1
    dk=wavv[2]-wavv[1]
    Nk=int(kmax/dk)+1
    plt.figure()
    plt.pcolor(wavv[numk//2:numk//2+Nk:],frqv[numf//2-Nf//2:numf//2+Nf//2:],np.log(FFK[numf//2-Nf//2:numf//2+Nf//2:,:Nk:]),cmap="Blues",rasterized=True,edgecolors='face')
    plt.colorbar()
    plt.xlim([0,kmax])
    format_fig(r"$k$ (1/m)",r"$f$ (Hz)")
    return
