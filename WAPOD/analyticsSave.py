#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 13:57:13 2023

@author: michael
"""

#import os
#os.chdir("/home/michael/Documents/Post-doc IMI/Simulations/1D/1-SMI_CauchyError/")

import numpy as np
import inputReading
#import matplotlib.pyplot as plt
import source
import modulation as modulateP

#
def invPara(K,M):
    if K==0:
        omeginvc=0
        invK=0
    else:
        omeginvc=np.sqrt(M/4/K)
        invK=1/K
    return invK,omeginvc

def fR01(Z,rho,cm,omega,M,K):    
    invK,omeginvc=invPara(K,M)
    denI=(Z[0]+Z[1])*(1-(omega*omeginvc)**2)+1j*omega*(Z[0]*Z[1]*invK+M)
    numR=(-Z[0]+Z[1])*(1-(omega*omeginvc)**2)+1j*omega*(-Z[0]*Z[1]*invK+M)
    FR01=numR/denI
    return FR01

def fT01(Z,rho,cm,omega,M,K):    
    invK,omeginvc=invPara(K,M)
    denI=(Z[0]+Z[1])*(1-(omega*omeginvc)**2)+1j*omega*(Z[0]*Z[1]*invK+M)
    numT=2*rho[0]*cm[1]*(1+(omega*omeginvc)**2)
    FT01=numT/denI
    return FT01

def analytics1InterfaceCauchy(Tmax,X,sourcef,frontiere,NFourier,mat,Kf):
    Force=float(sourcef["Force"])
    f0=float(sourcef["Frequence"])
    Fmax=Kf*f0
    DFourier=Fmax/NFourier
    Alpha=float(frontiere["Alpha_0"])
    Nmat,rho,cm=inputReading.material(mat)
    Z=rho*cm
    K=float(frontiere["Kn_0"])
    M=float(frontiere["Mn_0"])
    TempV=0
    TempS=0
    sce=source.jkps_t
    for J in range(NFourier):
        F = J * DFourier
        Omega = 2 * np.pi * F
        TF=2*source.TFsce(sce,J,NFourier,f0,Kf,Tmax)
        K0= Omega / cm[0]
        K1= Omega / cm[1]
        R01=fR01(Z,rho,cm,Omega,M,K)
        T01=fT01(Z,rho,cm,Omega,M,K)
        if sourcef["Onde"]=="V":
            cste=1
        else:
            cste=1/rho[0]/cm[0]
        if X < Alpha:
            # onde incidente
            sftI=Omega * Tmax - K0 * X
            expShiftI= np.exp(1j*sftI)
            VI=1/cm[0]* TF*expShiftI
            SI=-rho[0]* TF*expShiftI
            # onde reflechie
            sftR=Omega * Tmax + K0 * (X-2*Alpha)
            expShiftR= np.exp(1j*sftR)
            VR=-1/cm[0]* TF*expShiftR*R01
            SR=-rho[0]* TF*expShiftR*R01
            # integration numerique
            if ((J == 0) or (J == NFourier)):
                TempV= TempV + 0.5 * np.real(VI+VR)
                TempS= TempS + 0.5 * np.real(SI+SR)
            else:
                TempV= TempV + 1 * np.real(VI+VR)
                TempS= TempS + 1 * np.real(SI+SR)
        else:
            # onde transmise
            sftT=Omega * Tmax - K0*Alpha-K1*(X-Alpha)
            expShiftT= np.exp(1j*sftT)
            VT=1/cm[1]* TF*expShiftT*T01
            ST=-rho[1]* TF*expShiftT*T01
            # integration numerique
            if ((J == 1) or (J == NFourier)):
                TempV= TempV + 0.5 * np.real(VT)
                TempS= TempS + 0.5 * np.real(ST)
            else:
                TempV= TempV + 1 * np.real(VT)
                TempS= TempS + 1 * np.real(ST)
    Kte = DFourier * Force
    SolV = TempV * Kte * cste
    SolS = TempS * Kte * cste
    return SolV,SolS

def analytics1InterfacePSV(Tmax,X,sourcef,frontiere,NFourier,mat,Kf):
    Force=float(sourcef["Force"])
    f0=float(sourcef["Frequence"])
    x0=float(sourcef["Xsource"])
    Fmax=Kf*f0
    DFourier=Fmax/NFourier
    Alpha=float(frontiere["Alpha_0"])
    Nmat,rho,cm=inputReading.material(mat)
    Z=rho*cm
    K=float(frontiere["Kn1"])
    M=float(frontiere["Mn1"])
    TempV=0
    TempS=0
    sce=source.jkps_t
    for J in range(NFourier):
        F = J * DFourier
        Omega = 2 * np.pi * F
        TF=2*source.TFsce(sce,J,NFourier,f0,Kf,Tmax)
        K0= Omega / cm[0]
        K1= Omega / cm[1]
        R01=fR01(Z,rho,cm,Omega,M,K)
        T01=fT01(Z,rho,cm,Omega,M,K)
        if sourcef["Onde"]=="V":
            cste=1
        else:
            cste=-1/rho[0]/cm[0]
        if X < Alpha:
            # onde incidente
            sftI=Omega * Tmax - K0 * np.abs(X-x0)
            expShiftI= np.exp(1j*sftI)
            VI=1/cm[0]/2* TF*expShiftI
            SI=-rho[0]/2* TF*expShiftI*np.sign(X-x0)
            # onde reflechie
            sftR=Omega * Tmax + K0 * (X+x0-2*Alpha)
            expShiftR= np.exp(1j*sftR)
            VR=-1/cm[0]/2* TF*expShiftR*R01
            SR=-rho[0]/2* TF*expShiftR*R01
            # integration numerique
            if ((J == 0) or (J == NFourier)):
                TempV= TempV + 0.5 * np.real(VI+VR)
                TempS= TempS + 0.5 * np.real(SI+SR)
            else:
                TempV= TempV + 1 * np.real(VI+VR)
                TempS= TempS + 1 * np.real(SI+SR)
        else:
            # onde transmise
            sftT=Omega * Tmax - K0 *(Alpha-x0)-K1*(X-Alpha)
            expShiftT= np.exp(1j*sftT)
            VT=1/cm[1]/2* TF*expShiftT*T01
            ST=-rho[1]/2* TF*expShiftT*T01
            # integration numerique
            if ((J == 1) or (J == NFourier)):
                TempV= TempV + 0.5 * np.real(VT)
                TempS= TempS + 0.5 * np.real(ST)
            else:
                TempV= TempV + 1 * np.real(VT)
                TempS= TempS + 1 * np.real(ST)
    Kte = DFourier * Force
    SolV = TempV * Kte * cste
    SolS = TempS * Kte * cste
    return SolV,SolS

def invRiemannR(U, S, Rho, Cel):
    iRR= 0.5 * (U - S / (Rho * Cel))
    return iRR
def invRiemannL(U, S, Rho, Cel):
    iRL= 0.5 * (U + S / (Rho * Cel))
    return iRL

def functCauchyK(mat,modulation,frontiere,sourcef,KteV=1,KteS=-1):
    Nmat,rho,cm=inputReading.material(mat)
    alpha=inputReading.frontiere(frontiere)
    namorK = "beta_"+str(0)
    if frontiere.get(namorK)==None:
        etaK=0.
    else:
        etaK= float(frontiere[namorK])
    rc0=rho[0]*cm[0]
    rc1=rho[1]*cm[1]
    Tinit=float(sourcef["Tshift"])
    if float(frontiere["Kn_0"])==0:
        Kn=np.infty
    else:
        Kn=float(frontiere["Kn_0"])
    def fct(y,t):
        Xs=alpha[0]-cm[0]*(t-Tinit)
        sce=source.choice_timefct(sourcef, Tinit-Xs/cm[0])
        U=KteV/cm[0]*1*sce
        S=KteS*rho[0]*1*sce
        K=Kn/(1+modulateP.fonctionModK(t-Tinit,alpha[0],modulation,frontiere,0))
        Cprime=1/Kn*modulateP.fonctionModKDer(t-Tinit,alpha[0],modulation,frontiere,0)
        fac=1+rc1/rc0+rc1*(Cprime+etaK)
        return K/rc1*(2*invRiemannR(U, S, rho[0], cm[0])-fac*y)
    return fct

def functSourceK(mat,modulation,frontiere,sourcef,KteV=1,KteS=-1):
    Nmat,rho,cm=inputReading.material(mat)
    alpha=inputReading.frontiere(frontiere)
    namorK = "beta_"+str(0)
    if frontiere.get(namorK)==None:
        etaK=0.
    else:
        etaK= float(frontiere[namorK])
    rc0=rho[0]*cm[0]
    rc1=rho[1]*cm[1]
    x0=float(sourcef["Xsource"])
    Force=float(sourcef["Force"])
    if float(frontiere["Kn_0"])==0:
        Kn=np.infty
    else:
        Kn=float(frontiere["Kn_0"])
    def fct(y,t):
        sce=source.choice_timefct(sourcef, t-(alpha[0]-x0)/cm[0])
        U=Force/cm[0]*sce
        #S=KteS*rho[0]*1*sce
        K=Kn/(1+modulateP.fonctionModK(t,alpha[0],modulation,frontiere,0))
        Cprime=1/Kn*modulateP.fonctionModKDer(t,alpha[0],modulation,frontiere,0)
        fac=1+rc1/rc0+rc1*(Cprime+etaK)
        return K/rc1*(U-fac*y)
    return fct


def functCauchyM(mat,modulation,frontiere,sourcef,KteV=1,KteS=-1):
    Nmat,rho,cm=inputReading.material(mat)
    alpha=inputReading.frontiere(frontiere)
    rc0=rho[0]*cm[0]
    rc1=rho[1]*cm[1]
    namorM = "zeta_"+str(0)
    if frontiere.get(namorM)==None:
        etaM=0.
    else:
        etaM= float(frontiere[namorM])
    Tinit=float(sourcef["Tshift"])
    Mn=float(frontiere["Mn_0"])
    def fct(z,t):
        Xs=alpha[0]-cm[0]*(t-Tinit)
        sce=source.choice_timefct(sourcef, Tinit-Xs/cm[0])
        U=KteV/cm[0]*1*sce
        S=KteS*rho[0]*1*sce
        M=Mn*(1+modulateP.fonctionModM(t-Tinit,alpha[0],modulation,frontiere,0))
        Mprime=Mn*modulateP.fonctionModMDer(t-Tinit,alpha[0],modulation,frontiere,0)
        fac=1+rc0/rc1+(Mprime+etaM)/rc1
        return -rc1/M*(2*rc0*invRiemannR(U, S, rho[0], cm[0])+fac*z)
    return fct

def functSourceM(mat,modulation,frontiere,sourcef,KteV=1,KteS=-1):
    Nmat,rho,cm=inputReading.material(mat)
    alpha=inputReading.frontiere(frontiere)
    namorM = "zeta_"+str(0)
    if frontiere.get(namorM)==None:
        etaM=0.
    else:
        etaM= float(frontiere[namorM])
    rc0=rho[0]*cm[0]
    rc1=rho[1]*cm[1]
    x0=float(sourcef["Xsource"])
    Force=float(sourcef["Force"])
    Mn=float(frontiere["Mn_0"])
    def fct(z,t):
        sce=source.choice_timefct(sourcef, t-(alpha[0]-x0)/cm[0])
        U=Force*rho[0]*sce
        #S=KteS*rho[0]*1*sce
        M=Mn*(1+modulateP.fonctionModM(t,alpha[0],modulation,frontiere,0))
        Mprime=Mn*modulateP.fonctionModMDer(t,alpha[0],modulation,frontiere,0)
        fac=1+rc0/rc1+(Mprime+etaM)/rc1
        return -rc1/M*(U+fac*z)
    return fct
# def functSource(y,t,modulation,frontiere,i):
#     r0c0=rho[0]*cm[0]
#     r1c1=rho[1]*cm[1]
#     K=fctK(t,modulation)
#     invKprim=invKp(t,modulation)
#     fct=K/r1c1*(ht(t)+(1+r1c1/r0c0+r1c1*invKprim(t)))*y


    
def preProcessSMM(modulation, mat,frontiere,sourcef,dt,Nt,NOde):
    Kn=float(frontiere["Kn_0"])
    Mn=float(frontiere["Mn_0"])
    if sourcef['Forme'] == 'Cauchy':
        Tinit=float(sourcef["Tshift"])
    else:
        Tinit=0
    tRK=np.linspace(Tinit,Tinit+Nt*dt,NOde+1)
    TabUp=np.zeros(NOde+1)
    # for i in range(Nt):
    #     sce.append((sourcef,t[i]))
    if Kn != 0 and Mn==0:
        Vstart=0
        if sourcef['Forme'] == 'Cauchy':
            fctC=functCauchyK(mat,modulation,frontiere,sourcef)
            TabY=intRK4(fctC, Vstart,tRK) 
        if sourcef['Forme'] == 'point source':
            fctC=functSourceK(mat,modulation,frontiere,sourcef)
            TabY=intRK4(fctC, Vstart,tRK) 
        # else:
        #     TabY=intRK4(functSource, Vstart, tRK);
        for I in range(NOde+1):
            TabUp[I] = TabY[I]
    elif Mn!=0 and Kn==0:
        Vstart=0
        if sourcef['Forme'] == 'Cauchy':
            fctC=functCauchyM(mat,modulation,frontiere,sourcef)
            TabY=intRK4(fctC, Vstart,tRK) 
        if sourcef['Forme'] == 'point source':
            fctC=functSourceM(mat,modulation,frontiere,sourcef)
            TabY=intRK4(fctC, Vstart,tRK) 
        # else:
        #     TabY=intRK4(functSource, Vstart, tRK);
        for I in range(NOde+1):
            TabUp[I] = TabY[I]
    return TabUp
    
         
def intRK4(f,y0,t):
    y = np.zeros(len(t))
    y[0] = y0    
    for i in range(0,len(t)-1):
        h = t[i+1]-t[i]
        F1 = h*f(y[i],t[i])
        F2 = h*f((y[i]+F1/2),(t[i]+h/2))
        F3 = h*f((y[i]+F2/2),(t[i]+h/2))
        F4 = h*f((y[i]+F3),(t[i]+h))
        y[i+1] = y[i] + 1/6*(F1 + 2*F2 + 2*F3 + F4)
    return y     

def calculUP(T,Tinit,Dode,TabUp):
    J = int(np.trunc((T - Tinit) / Dode))
    Eps = (T - Tinit  - J * Dode) / Dode
    Result = (1 - Eps) * TabUp[J] + Eps * TabUp[J + 1]
    return Result

def analytics_Cauchy_1IModK(Tmx,X,sourcef,modulation,frontiere,mat,NOde,dt,Nt):
    Force=float(sourcef["Force"])
    f0=float(sourcef["Frequence"])
    KteV=1
    KteS=-1
    Tinit=float(sourcef["Tshift"])
    Tmax=Tmx+Tinit
    Alpha=float(frontiere["Alpha_0"])
    Nmat,rho,cm=inputReading.material(mat)
    TempV=0
    TempS=0
    DOde= Nt *dt / NOde
    sce=source.temporalFct[sourcef["Temporel"]]
    SolV =np.zeros(len(X))
    SolS=np.zeros(len(X))
    trace=preProcessSMM(modulation,mat,frontiere,sourcef,dt,Nt,NOde)
    for i in range(len(X)):
        if X[i] <= Alpha:
            F2=sce(Tinit - (X[i] - cm[0] * (Tmax - Tinit)) / cm[0],f0)
            TempV = KteV * 1 / cm[0] * F2 * Force
            TempS = KteS * rho[0] * F2 * Force
            J0R =  invRiemannR(TempV, TempS, rho[0], cm[0])
            Ta = Tmax - (Alpha - X[i]) / cm[0]
            if Ta >= Tinit:
                Up = calculUP(Ta,Tinit,DOde,trace)
                Sp = - rho[1] * cm[1] * Up
                Sm = Sp
                F2=sce(Tinit - (Alpha - cm[0] * (Ta - Tinit)) / cm[0],f0)
                TempV = KteV * 1 / cm[0] * F2 * Force
                TempS = KteS * rho[0] * F2 * Force
                Um = Sm / (rho[0] * cm[0]) + 2 * invRiemannR(TempV, TempS, rho[0], cm[0])
                J0L = invRiemannL(Um, Sm, rho[0], cm[0])
            else:
                F2=sce(Tinit - (X[i] + cm[0] * (Tmax - Tinit)) / cm[0],f0)
                TempV = KteV * 1 / cm[0] * F2 * Force
                TempS = KteS * rho[0] * F2 * Force
                J0L = invRiemannL(TempV, TempS, rho[0], cm[0])
            SolV[i] = J0R + J0L
            SolS[i] = rho[0] *cm[0] * (J0L - J0R)
        else:
            Tb = Tmax - (X[i] - Alpha) / cm[1]
            if Tb >= Tinit:
                Up = calculUP(Tb,Tinit,DOde,trace)
            else:
                Up = 0.
            SolV[i] = Up
            SolS[i] = - rho[1] * cm[1] * Up
    return SolV,SolS

def analytics_Cauchy_1IModM(Tmx,X,sourcef,modulation,frontiere,mat,NOde,dt,Nt):
    Force=float(sourcef["Force"])
    f0=float(sourcef["Frequence"])
    KteV=1
    KteS=-1
    Tinit=float(sourcef["Tshift"])
    Tmax=Tmx+Tinit
    Alpha=float(frontiere["Alpha_0"])
    Nmat,rho,cm=inputReading.material(mat)
    TempV=0
    TempS=0
    DOde= Nt *dt / NOde
    sce=source.temporalFct[sourcef["Temporel"]]
    SolV =np.zeros(len(X))
    SolS=np.zeros(len(X))
    trace=preProcessSMM(modulation,mat,frontiere,sourcef,dt,Nt,NOde)
    for i in range(len(X)):
        if X[i] <= Alpha:
            F2=sce(Tinit - (X[i] - cm[0] * (Tmax - Tinit)) / cm[0],f0)
            TempV = KteV * 1 / cm[0] * F2 * Force
            TempS = KteS * rho[0] * F2 * Force
            J0R =  invRiemannR(TempV, TempS, rho[0], cm[0])
            Ta = Tmax - (Alpha - X[i]) / cm[0]
            if Ta >= Tinit:
                Sp = calculUP(Ta,Tinit,DOde,trace)
                Up = - 1/rho[1]/ cm[1] * Sp
                Um = Up
                F2=sce(Tinit - (Alpha - cm[0] * (Ta - Tinit)) / cm[0],f0)
                TempV = KteV * 1 / cm[0] * F2 * Force
                TempS = KteS * rho[0] * F2 * Force
                Sm = Um * (rho[0] * cm[0]) - 2 * rho[0] * cm[0] * invRiemannR(TempV, TempS, rho[0], cm[0])
                J0L = invRiemannL(Um, Sm, rho[0], cm[0])
            else:
                F2=sce(Tinit - (X[i] + cm[0] * (Tmax - Tinit)) / cm[0],f0)
                TempV = KteV * 1 / cm[0] * F2 * Force
                TempS = KteS * rho[0] * F2 * Force
                J0L = invRiemannL(TempV, TempS, rho[0], cm[0])
            SolV[i] =  (J0R + J0L)
            SolS[i] = rho[0]*cm[0] *(J0L - J0R)
        else:
            Tb = Tmax - (X[i] - Alpha) / cm[1]
            if Tb >= Tinit:
                Up = calculUP(Tb,Tinit,DOde,trace)
            else:
                Up = 0.
            SolV[i] = -1/rho[1] / cm[1]*Up
            SolS[i] =  Up
    return SolV,SolS
    # return trace
    
    
def analytics_SourceV_1IModK(Tmx,X,sourcef,modulation,frontiere,mat,NOde,dt,Nt):
    Tinit=0
    Force=float(sourcef["Force"])
    f0=float(sourcef["Frequence"])
    x0=float(sourcef["Xsource"])
    Tmax=Tmx
    Alpha=float(frontiere["Alpha_0"])
    Nmat,rho,cm=inputReading.material(mat)
    DOde= Nt *dt / NOde
    sce=source.temporalFct[sourcef["Temporel"]]
    SolV =np.zeros(len(X))
    SolS=np.zeros(len(X))
    trace=preProcessSMM(modulation,mat,frontiere,sourcef,dt,Nt,NOde)
    for i in range(len(X)):
        if X[i] <= x0:
            Ta = Tmax - (Alpha - X[i]) / cm[0]
            if Ta >= 0:
                Up = calculUP(Ta,Tinit,DOde,trace)
                DA1 = - rho[1] * cm[1] /rho[0] /cm[0]* Up
                F2=sce(Tmax - (Alpha -X[i]+Alpha-x0)/ cm[0],f0)
                DA2=Force/2/cm[0]*F2
                DA=DA1+DA2
            else:
                DA=0
            F2=sce(Tmax - np.abs(X[i]-x0)/ cm[0],f0)
            Vs=Force/2/cm[0]*F2
            SolV[i]=DA+Vs
            SolS[i]=rho[0]*cm[0]*(DA+Vs)
        elif X[i]<Alpha:
            Ta = Tmax - (Alpha - X[i]) / cm[0]
            if Ta>=0:
                Up = calculUP(Ta,Tinit,DOde,trace)
                DA1 = - rho[1] * cm[1] /rho[0] / cm[0]* Up
                F2=sce(Tmax - (Alpha -X[i]+Alpha-x0)/ cm[0],f0)
                DA2=Force/2/cm[0]*F2
                DA=DA1+DA2
            else:
                DA=0
            F2=sce(Tmax - (X[i]-x0)/ cm[0],f0)
            Vs=Force/2/cm[0]*F2
            SolV[i]=DA+Vs
            SolS[i]=rho[0]*cm[0]*(DA-Vs)
        else:
            Tb=Tmax-(X[i]-Alpha)/cm[1]
            if Tb>=0:
                Up = calculUP(Tb,Tinit,DOde,trace)
            else:
                Up=0.
            SolV[i]=Up
            SolS[i]=-rho[1]*cm[1]*Up
    return SolV,SolS

def analytics_SourceV_1IModM(Tmx,X,sourcef,modulation,frontiere,mat,NOde,dt,Nt):
    Tinit=0
    Force=float(sourcef["Force"])
    f0=float(sourcef["Frequence"])
    x0=float(sourcef["Xsource"])
    Tmax=Tmx
    Alpha=float(frontiere["Alpha_0"])
    Nmat,rho,cm=inputReading.material(mat)
    DOde= Nt *dt / NOde
    sce=source.temporalFct[sourcef["Temporel"]]
    SolV =np.zeros(len(X))
    SolS=np.zeros(len(X))
    trace=preProcessSMM(modulation,mat,frontiere,sourcef,dt,Nt,NOde)
    for i in range(len(X)):
        if X[i] <= x0:
            Ta = Tmax - (Alpha - X[i]) / cm[0]
            if Ta >= 0:
                Sp = calculUP(Ta,Tinit,DOde,trace)
                DA1 = - 1/rho[1] / cm[1] *Sp
                F2=sce(Tmax - (Alpha -X[i]+Alpha-x0)/ cm[0],f0)
                DA2=-Force/2/cm[0]*F2
                DA=DA1+DA2
            else:
                DA=0
            F2=sce(Tmax - np.abs(X[i]-x0)/ cm[0],f0)
            Vs=Force/2/cm[0]*F2
            SolV[i]=DA+Vs
            SolS[i]=rho[0]*cm[0]*(DA+Vs)
        elif X[i]<Alpha:
            Ta = Tmax - (Alpha - X[i]) / cm[0]
            if Ta>=0:
                Sp = calculUP(Ta,Tinit,DOde,trace)
                DA1 = - 1/rho[1] / cm[1] * Sp
                F2=sce(Tmax - (Alpha -X[i]+Alpha-x0)/ cm[0],f0)
                DA2=-Force/2/cm[0]*F2
                DA=DA1+DA2
            else:
                DA=0
            F2=sce(Tmax - (X[i]-x0)/ cm[0],f0)
            Vs=Force/2/cm[0]*F2
            SolV[i]=DA+Vs
            SolS[i]=rho[0]*cm[0]*(DA-Vs)
        else:
            Tb=Tmax-(X[i]-Alpha)/cm[1]
            if Tb>=0:
                Sp = calculUP(Tb,Tinit,DOde,trace)
            else:
                Sp=0.
            SolV[i]=-1/rho[1]/cm[1]*Sp
            SolS[i]=Sp
    return SolV,SolS