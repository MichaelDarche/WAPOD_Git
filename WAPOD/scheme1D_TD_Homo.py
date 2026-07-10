#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 11:32:11 2025

@author: michael
"""


import numpy as np
import sympy as sp
import source

def m1Kmod(modulation,frontiere):
    t=sp.symbols("t")
    contacti="Contact_"+str(1)
    h=float(frontiere["Alpha_1"])-float(frontiere["Alpha_0"])
    alpha=float(frontiere["Alpha_0"])
    if frontiere[contacti]=="masse-ressort module":
        ampF=float(modulation["DeltaF"])
        ampM=float(modulation["DeltaM"])
        if modulation["Synchro"]=="Yes":
            FreqM=float(modulation["Freq"])
            FreqF=float(modulation["Freq"])
            FreqB=float(modulation["Freq"])
            FreqZ=float(modulation["Freq"])
            celM=float(modulation["ModC"])
            celF=float(modulation["ModC"])
            celB=float(modulation["ModC"])
            celZ=float(modulation["ModC"])
            nfrB = "beta_"+str(0)
            if frontiere.get(nfrB)==None:
                ampB=0.
            else:
                ampB= float(modulation["DeltaB"])
            nfrZ = "zeta_"+str(0)
            if frontiere.get(nfrZ)==None:
                ampZ=0.
            else:
                ampZ= float(modulation["DeltaZ"])
        else:
            FreqM=float(modulation["FreqM"])
            FreqF=float(modulation["FreqF"])
            nfrB = "beta_"+str(0)
            if frontiere.get(nfrB)==None:
                FreqB=0.
                ampB=0.
            else:
                FreqB= float(modulation["FreqB"])
                ampB= float(modulation["DeltaB"])
            nfrZ = "zeta_"+str(0)
            if frontiere.get(nfrZ)==None:
                FreqZ=0.
                ampZ=0.
            else:
                FreqZ= float(modulation["FreqZ"])
                ampZ= float(modulation["DeltaZ"])
            celM=float(modulation["ModM"])
            celF=float(modulation["ModF"])
            celB=float(modulation["ModB"])
            celZ=float(modulation["ModZ"])
        if celM==0:
            phiMm=0
        else:
            phiMm=FreqM*2*np.pi*alpha/celM
        if celF==0:
            phiFm=0
        else:
            phiFm=FreqF*2*np.pi*alpha/celF
        if celB==0:
            phiBm=0
        else:
            phiBm=FreqB*2*np.pi*alpha/celB
        if celZ==0:
            phiZm=0
        else:
            phiZm=FreqZ*2*np.pi*alpha/celZ
        if modulation["TypeMod"]=="Sinus":
            fctM=sp.sin(2*np.pi*FreqM*t+phiMm)
            fctF=sp.sin(2*np.pi*FreqF*t+phiFm)
            fctB=sp.sin(2*np.pi*FreqB*t+phiBm)
            fctZ=sp.sin(2*np.pi*FreqZ*t+phiZm)
        elif modulation["TypeMod"]=="Square":
            fctM=sp.sign(sp.sin((2*np.pi*FreqM*t+phiMm)))
            fctF=sp.sign(sp.sin((2*np.pi*FreqF*t+phiFm)))
            fctB=sp.sign(sp.sin((2*np.pi*FreqB*t+phiBm)))
            fctZ=sp.sign(sp.sin((2*np.pi*FreqZ*t+phiZm)))
        elif modulation["TypeMod"]=="Square NS":
            ratio=float(modulation["Ratio"])
            fctM=-sp.sign(FreqM*t-sp.floor(FreqM*t)-ratio)
            fctF=-sp.sign(FreqF*t-sp.floor(FreqF*t)-ratio)
            fctB=-sp.sign(FreqB*t-sp.floor(FreqB*t)-ratio)
            fctZ=-sp.sign(FreqZ*t-sp.floor(FreqZ*t)-ratio)
        nameKni="Kn_"+str(1)
        nameMni="Mn_"+str(1)
        K=float(frontiere[nameKni])
        if K==0:
            Fmod=0
        else:
            Fmod=(1+ampF*fctF)/K
        M=float(frontiere[nameMni])
        namorK = "beta_"+str(1)
        if frontiere.get(namorK)==None:
            etaK=0.
        else:
            etaK= float(frontiere[namorK])
        namorM = "beta_"+str(1)
        if frontiere.get(namorK)==None:
            etaM=0.
        else:
            etaM= float(frontiere[namorM])
        Mmod=(1+ampM*fctM)*M
        Bmod=(1+ampB*fctB)*etaK
        Zmod=(1+ampZ*fctZ)*etaM
    if frontiere[contacti]=="masse-ressort module dephase":
        ampF=float(modulation["DeltaF"])
        ampM=float(modulation["DeltaM"])
        ampF=float(modulation["DeltaF"])
        ampM=float(modulation["DeltaM"])
        if modulation["Synchro"]=="Yes":
            Freq=float(modulation["Freq"])
            cel=float(modulation["ModC"])
            if cel==0:
                phim=0
            else:
                phim=Freq/cel
            Nip=int(cel/h/Freq)
        else:
            print("Homogenization not ready...")
        fct=[]
        if modulation["TypeMod"]=="Sinus":
            for i in range(Nip):
                fct.append(sp.sin(2*np.pi*(Freq*t+phim*(alpha+i*h))))
        elif modulation["TypeMod"]=="Square":
            for i in range(Nip):
                fct.append(sp.sign(sp.sin(2*np.pi*(Freq*t+phim*(alpha+i*h)))))
        elif modulation["TypeMod"]=="Square NS":
            ratio=float(modulation["Ratio"])
            for i in range(Nip):
                sigph=Freq*t+phim*(alpha+i*h)%1
                fct.append(-sp.sign(sigph-ratio))
        K=[]
        M=[]
        QK=[]
        QM=[]
        Mmod=0
        Fmod=0
        Bmod=0
        Zmod=0
        for i in range(Nip):
            nameKni="Kn_"+str(i)
            nameMni="Mn_"+str(i)
            namorK = "beta_"+str(i)
            if frontiere.get(namorK)==None:
                etaK=0.
            else:
                etaK= float(frontiere[namorK])
            namorM = "beta_"+str(i)
            if frontiere.get(namorK)==None:
                etaM=0.
            else:
                etaM= float(frontiere[namorM])
            K.append(float(frontiere[nameKni]))
            M.append(float(frontiere[nameMni]))
            QK.append(etaK)
            QM.append(etaM)
            if K[i]==0:
                Fmod=Fmod
            else:
                Fmod=Fmod+1/K[i]*(1+ampF*fct[i])
            Mmod=Mmod+M[i]*(1+ampM*fct[i])
            Bmod=Bmod+QK[i]*(1+ampB*fct[i])
            Zmod=Zmod+QM[i]*(1+ampZ*fct[i])            
    return Fmod

def m1KMmod(modulation,frontiere):
    t=sp.symbols("t")
    contacti="Contact_"+str(1)
    h=float(frontiere["Alpha_1"])-float(frontiere["Alpha_0"])
    alpha=float(frontiere["Alpha_0"])
    if frontiere[contacti]=="masse-ressort module":
        ampF=float(modulation["DeltaF"])
        ampM=float(modulation["DeltaM"])
        if modulation["Synchro"]=="Yes":
            FreqM=float(modulation["Freq"])
            FreqF=float(modulation["Freq"])
            FreqB=float(modulation["Freq"])
            FreqZ=float(modulation["Freq"])
            celM=float(modulation["ModC"])
            celF=float(modulation["ModC"])
            celB=float(modulation["ModC"])
            celZ=float(modulation["ModC"])
            nfrB = "beta_"+str(0)
            if frontiere.get(nfrB)==None:
                ampB=0.
            else:
                ampB= float(modulation["DeltaB"])
            nfrZ = "zeta_"+str(0)
            if frontiere.get(nfrZ)==None:
                ampZ=0.
            else:
                ampZ= float(modulation["DeltaZ"])
        else:
            FreqM=float(modulation["FreqM"])
            FreqF=float(modulation["FreqF"])
            nfrB = "beta_"+str(0)
            if frontiere.get(nfrB)==None:
                FreqB=0.
                ampB=0.
            else:
                FreqB= float(modulation["FreqB"])
                ampB= float(modulation["DeltaB"])
            nfrZ = "zeta_"+str(0)
            if frontiere.get(nfrZ)==None:
                FreqZ=0.
                ampZ=0.
            else:
                FreqZ= float(modulation["FreqZ"])
                ampZ= float(modulation["DeltaZ"])
            celM=float(modulation["ModM"])
            celF=float(modulation["ModF"])
            celB=float(modulation["ModB"])
            celZ=float(modulation["ModZ"])
        if celM==0:
            phiMm=0
        else:
            phiMm=FreqM*2*np.pi*alpha/celM
        if celF==0:
            phiFm=0
        else:
            phiFm=FreqF*2*np.pi*alpha/celF
        if celB==0:
            phiBm=0
        else:
            phiBm=FreqB*2*np.pi*alpha/celB
        if celZ==0:
            phiZm=0
        else:
            phiZm=FreqZ*2*np.pi*alpha/celZ
        if modulation["TypeMod"]=="Sinus":
            fctM=sp.sin(2*np.pi*FreqM*t+phiMm)
            fctF=sp.sin(2*np.pi*FreqF*t+phiFm)
            fctB=sp.sin(2*np.pi*FreqB*t+phiBm)
            fctZ=sp.sin(2*np.pi*FreqZ*t+phiZm)
        elif modulation["TypeMod"]=="Square":
            fctM=sp.sign(sp.sin((2*np.pi*FreqM*t+phiMm)))
            fctF=sp.sign(sp.sin((2*np.pi*FreqF*t+phiFm)))
            fctB=sp.sign(sp.sin((2*np.pi*FreqB*t+phiBm)))
            fctZ=sp.sign(sp.sin((2*np.pi*FreqZ*t+phiZm)))
        elif modulation["TypeMod"]=="Square NS":
            ratio=float(modulation["Ratio"])
            fctM=-sp.sign(FreqM*t-sp.floor(FreqM*t)-ratio)
            fctF=-sp.sign(FreqF*t-sp.floor(FreqF*t)-ratio)
            fctB=-sp.sign(FreqB*t-sp.floor(FreqB*t)-ratio)
            fctZ=-sp.sign(FreqZ*t-sp.floor(FreqZ*t)-ratio)
        nameKni="Kn_"+str(1)
        nameMni="Mn_"+str(1)
        K=float(frontiere[nameKni])
        if K==0:
            Fmod=0
        else:
            Fmod=(1+ampF*fctF)/K
        M=float(frontiere[nameMni])
        namorK = "beta_"+str(1)
        if frontiere.get(namorK)==None:
            etaK=0.
        else:
            etaK= float(frontiere[namorK])
        namorM = "beta_"+str(1)
        if frontiere.get(namorK)==None:
            etaM=0.
        else:
            etaM= float(frontiere[namorM])
        Mmod=(1+ampM*fctM)*M
        Bmod=(1+ampB*fctB)*etaK
        Zmod=(1+ampZ*fctZ)*etaM
    if frontiere[contacti]=="masse-ressort module dephase":
        ampF=float(modulation["DeltaF"])
        ampM=float(modulation["DeltaM"])
        ampF=float(modulation["DeltaF"])
        ampM=float(modulation["DeltaM"])
        if modulation["Synchro"]=="Yes":
            Freq=float(modulation["Freq"])
            cel=float(modulation["ModC"])
            if cel==0:
                phim=0
            else:
                phim=Freq/cel
            Nip=int(cel/h/Freq)
        else:
            print("Homogenization not ready...")
        fct=[]
        if modulation["TypeMod"]=="Sinus":
            for i in range(Nip):
                fct.append(sp.sin(2*np.pi*(Freq*t+phim*(alpha+i*h))))
        elif modulation["TypeMod"]=="Square":
            for i in range(Nip):
                fct.append(sp.sign(sp.sin(2*np.pi*(Freq*t+phim*(alpha+i*h)))))
        elif modulation["TypeMod"]=="Square NS":
            ratio=float(modulation["Ratio"])
            for i in range(Nip):
                sigph=Freq*t+phim*(alpha+i*h)%1
                fct.append(-sp.sign(sigph-ratio))
        K=[]
        M=[]
        QK=[]
        QM=[]
        Mmod=0
        Fmod=0
        Bmod=0
        Zmod=0
        for i in range(Nip):
            nameKni="Kn_"+str(i)
            nameMni="Mn_"+str(i)
            namorK = "beta_"+str(i)
            if frontiere.get(namorK)==None:
                etaK=0.
            else:
                etaK= float(frontiere[namorK])
            namorM = "beta_"+str(i)
            if frontiere.get(namorK)==None:
                etaM=0.
            else:
                etaM= float(frontiere[namorM])
            K.append(float(frontiere[nameKni]))
            M.append(float(frontiere[nameMni]))
            QK.append(etaK)
            QM.append(etaM)
            if K[i]==0:
                Fmod=Fmod
            else:
                Fmod=Fmod+1/K[i]*(1+ampF*fct[i])
            Mmod=Mmod+M[i]*(1+ampM*fct[i])
            Bmod=Bmod+QK[i]*(1+ampB*fct[i])
            Zmod=Zmod+QM[i]*(1+ampZ*fct[i])        
    if frontiere[contacti]=="parfait":
        Fmod=0
        Mmod=0
        Bmod=0
        Zmod=0
    return Fmod,Mmod

def m1KMmodPeriod(frontiere,nbpCell):
    t=sp.symbols("t")
    h=float(frontiere["Alpha_"+str(nbpCell)])-float(frontiere["Alpha_0"])
    alpha=float(frontiere["Alpha_0"])
    ampFi=[]
    ampMi=[]
    ampBi=[]
    ampZi=[]
    Freqi=[]
    fctF=[]
    fctM=[]
    fctB=[]
    fctZ=[]
    PhiF0=[]
    PhiM0=[]
    PhiZ0=[]
    PhiB0=[]
    K=[]
    M=[]
    QM=[]
    QK=[]
    Fmod=0
    Bmod=0
    Mmod=0
    Zmod=0
    for i in range(nbpCell-1):
        ampFi.append(float(frontiere["DeltaF_"+str(i)]))
        ampMi.append(float(frontiere["DeltaM_"+str(i)]))
        nDeltaZ="DeltaZ_"+str(i)
        nDeltaB="DeltaB_"+str(i)
        if frontiere.get(nDeltaZ)==None:
            ampZi.append(0.)
            PhiZ0.append(0.)
        else:
            ampZi.append(sp.Float(frontiere[nDeltaZ]))
            PhiZ0.append(sp.Float(frontiere["phiZ0_"+str(i)]))
        if frontiere.get(nDeltaB)==None:
            ampBi.append(0.)
            PhiB0.append(0.)
        else:
            ampBi.append(sp.Float(frontiere[nDeltaB]))
            PhiB0.append(sp.Float(frontiere["phiB0_"+str(i)]))
        Freqi.append(float(frontiere["Freq_"+str(i)]))
        PhiF0.append(float(frontiere["phiF0_"+str(i)]))
        PhiM0.append(float(frontiere["phiM0_"+str(i)]))
        if frontiere["TypeMod_"+str(i)]=="Sinus":
            fctF.append(sp.sin(2*np.pi*(Freqi[i]*t+PhiF0[i])))
            fctM.append(sp.sin(2*np.pi*(Freqi[i]*t+PhiM0[i])))
            fctB.append(sp.sin(2*np.pi*(Freqi[i]*t+PhiB0[i])))
            fctZ.append(sp.sin(2*np.pi*(Freqi[i]*t+PhiZ0[i])))
        elif frontiere["TypeMod_"+str(i)]=="Cosinus":
            fctF.append(sp.cos(2*np.pi*(Freqi[i]*t+PhiF0[i])))
            fctM.append(sp.cos(2*np.pi*(Freqi[i]*t+PhiM0[i])))
            fctB.append(sp.cos(2*np.pi*(Freqi[i]*t+PhiB0[i])))
            fctZ.append(sp.cos(2*np.pi*(Freqi[i]*t+PhiZ0[i])))
        elif frontiere["TypeMod_"+str(i)]=="Square":
            fctF.append(sp.sign(sp.sin(2*np.pi*(Freqi[i]*t+PhiF0[i]))))
            fctM.append(sp.sign(sp.sin(2*np.pi*(Freqi[i]*t+PhiM0[i]))))
            fctB.append(sp.sign(sp.sin(2*np.pi*(Freqi[i]*t+PhiB0[i]))))
            fctZ.append(sp.sign(sp.sin(2*np.pi*(Freqi[i]*t+PhiZ0[i]))))
        elif frontiere["TypeMod_"+str(i)]=="Square NS":
            ratio=float(frontiere["Ratio_"+str(i)])
            sigph=Freqi[i]*t%1
            fctF.append(-sp.sign(sigph-ratio))
            fctM.append(-sp.sign(sigph-ratio))
            fctB.append(-sp.sign(sigph-ratio))
            fctZ.append(-sp.sign(sigph-ratio))
        nameKni="Kn_"+str(i)
        namorK = "beta_"+str(i)
        if frontiere.get(namorK)==None:
            etaK=0.
        else:
            etaK= float(frontiere[namorK])
        K.append(float(frontiere[nameKni]))
        QK.append(etaK)
        nameMni="Mn_"+str(i)
        namorM = "zeta_"+str(i)
        if frontiere.get(namorM)==None:
            etaM=0.
        else:
            etaM= float(frontiere[namorM])
        M.append(float(frontiere[nameMni]))
        QM.append(etaM)
        if K[i]==0:
            Fmod=Fmod
        else:
            Fmod=Fmod+1/K[i]*(1+ampFi[i]*fctF[i])
        Mmod=Mmod+M[i]*(1+ampMi[i]*fctM[i])
        Bmod=Bmod+QK[i]*(1+ampBi[i]*fctB[i]) 
        Zmod=Bmod+QM[i]*(1+ampZi[i]*fctZ[i])         
    return Fmod,Mmod,Bmod,Zmod

def m1KMmodP(modulation,frontiere,nbpCell=1):
    t=sp.symbols("t")
    ampFi=[]
    ampMi=[]
    ampBi=[]
    ampZi=[]
    Freqi=[]
    fct=[]
    K=[]
    M=[]
    QM=[]
    QK=[]
    Fmod=0
    Bmod=0
    Mmod=0
    Zmod=0
    for i in range(nbpCell):
        ampFi.append(float(modulation["DeltaF"]))
        ampMi.append(float(modulation["DeltaM"]))
        ampBi.append(float(modulation["DeltaF"]))
        ampZi.append(float(modulation["DeltaM"]))
        Freqi.append(float(modulation["Freq"]))
        if modulation["TypeMod"]=="Sinus":
            fct.append(sp.sin(2*np.pi*(Freqi[i]*t)))
        elif modulation["TypeMod"]=="Cosinus":
            fct.append(sp.cos(2*np.pi*(Freqi[i]*t)))
        elif modulation["TypeMod"]=="Square":
            fct.append(sp.sign(sp.sin(2*np.pi*(Freqi[i]*t))))
        elif modulation["TypeMod"]=="Square NS":
            ratio=float(frontiere["Ratio"])
            sigph=Freqi[i]*t%1
            fct.append(-sp.sign(sigph-ratio))
        nameKni="Kn_"+str(i)
        namorK = "beta_"+str(i)
        if frontiere.get(namorK)==None:
            etaK=0.
        else:
            etaK= float(frontiere[namorK])
        K.append(float(frontiere[nameKni]))
        QK.append(etaK)
        nameMni="Mn_"+str(i)
        namorM = "zeta_"+str(i)
        if frontiere.get(namorM)==None:
            etaM=0.
        else:
            etaM= float(frontiere[namorM])
        M.append(float(frontiere[nameMni]))
        QM.append(etaM)
        if K[i]==0:
            Fmod=Fmod
        else:
            Fmod=Fmod+1/K[i]*(1+ampFi[i]*fct[i])
        Mmod=Mmod+M[i]*(1+ampMi[i]*fct[i])
        Bmod=Bmod+QK[i]*(1+ampBi[i]*fct[i]) 
        Zmod=Bmod+QM[i]*(1+ampZi[i]*fct[i])         
    return Fmod,Mmod,Bmod,Zmod

def mat0(modulation,frontiere,milieu):
    h=float(frontiere["Alpha_1"])-float(frontiere["Alpha_0"])
    rho=float(milieu["Rho_0"])
    c=float(milieu["Cel_0"])
    E=rho*c*c
    Fmod=m1Kmod(modulation,frontiere)
    E0=1/(Fmod/h+1/E)
    return rho,E0

    
def mat00(modulation,frontiere,milieu):
    h=float(frontiere["Alpha_1"])-float(frontiere["Alpha_0"])
    rho=float(milieu["Rho_0"])
    c=float(milieu["Cel_0"])
    E=rho*c*c
    Fmod,Mmod,Bmod,Zmod=m1KMmodP(modulation,frontiere)
    E0=1/(Fmod/h+1/E)
    rho0=rho+Mmod/h
    return rho0,E0

def mat00Period(frontiere,nbpCell,milieu):
    h=float(frontiere["Alpha_"+str(nbpCell)])-float(frontiere["Alpha_0"])
    rho=float(milieu["Rho_0"])
    c=float(milieu["Cel_0"])
    E=rho*c*c
    Fmod,Mmod,Bmod,Zmod=m1KMmodPeriod(frontiere,nbpCell)
    E0=1/(Fmod/h+1/E)
    rho0=rho+Mmod/h
    return rho0,E0

def CPP1(modulation,frontiere,milieu,ordre0="no"):
    h=float(frontiere["Alpha_1"])-float(frontiere["Alpha_0"])
    rho=float(milieu["Rho_0"])
    c=float(milieu["Cel_0"])
    E=rho*c*c
    Fmod = m1Kmod(modulation,frontiere)
    if ordre0=="yes":
        CPP1=0
    else:
        CPP1=-1+1/(1+E/h*Fmod)
    return CPP1

def B1B3B5(modulation,frontiere,milieu,ordre0="no"):
    t=sp.symbols("t")
    CDP=CPP1(modulation,frontiere,milieu,ordre0)
    B1=1/12*(CDP)**2
    if modulation["TypeMod"]=="Square" or modulation["TypeMod"]=="Square NS":
        CDPder=0
    else:
        CDPder=sp.diff(CDP,t)
    B3=1/6*CDP*CDPder
    B5=1/12*CDP*sp.diff(CDPder,t)
    return B1,B3,B5


def coefsbc(modulation,frontiere,milieu,ordre0="no"):
    h=float(frontiere["Alpha_1"])-float(frontiere["Alpha_0"])
    B1,B3,B5=B1B3B5(modulation,frontiere,milieu,ordre0)
    rho,E0=mat0(modulation,frontiere,milieu)
    coefc=sp.sqrt(h**2*B5+E0/rho)
    b1=h**2*B1
    b3=h**2*B3
    return coefc,b1,b3

def minmaxbc(modulation,frontiere,mat,ordre0="no",nt=4000):
    t=sp.symbols('t')
    fk=float(modulation["Freq"])
    coefc,b1,b3=coefsbc(modulation,frontiere,mat,ordre0)
    t=sp.symbols('t')
    funcc=sp.lambdify(t,coefc,modules="numpy")
    funcb1=sp.lambdify(t,b1, modules="numpy")
    maxc=funcc(0)
    minb1=funcb1(0)
    if fk!=0:
        t=np.linspace(0,1/fk,num=nt)
        for i in range(nt):
            maxc=max(maxc,funcc(t[i]))
            minb1=min(minb1,funcb1(t[i]))
    print(minb1,maxc)
    return minb1,maxc

def deltat(dx,modulation,frontiere,mat,nt=1000,CFL=0.95,dtopti="no",ordre0="no"):
    minb1,maxc=minmaxbc(modulation,frontiere,mat,ordre0)
    if dtopti=="yes":
        deltaT=np.sqrt(dx**2+4*minb1)/maxc
    else:
        deltaT=CFL*dx/maxc
    return deltaT

def deltatp(dx,frontiere,nbpCell,mat,nt=1000,CFL=0.95):
    t=sp.symbols('t')
    rho,E=mat00Period(frontiere,nbpCell,mat)
    c=sp.sqrt(E/rho)
    funcc=sp.lambdify(t,c,modules="numpy")
    maxc=0
    fk=float(frontiere["Freq_0"])
    for i in range(nbpCell-1):
        fk=min(fk,float(frontiere["Freq_"+str(i)]))
    if fk!=0:
        t=np.linspace(0,1/fk,num=nt)
        for i in range(nt):
            maxc=max(maxc,funcc(t[i]))
    deltaT=CFL*dx/maxc
    return deltaT
    

def coefsABC(modulation,frontiere,milieu,DeltaX,DeltaT,ordre0="no"):
    coefc,b1,b3=coefsbc(modulation,frontiere,milieu,ordre0)
    A=coefc*DeltaT/DeltaX
    B=b1/DeltaX**2
    C=b3*DeltaT/DeltaX**2
    return A,B,C

def funcPQ(modulation,frontiere,milieu,DeltaX,DeltaT,ordre0="no"):
    A,B,C=coefsABC(modulation,frontiere,milieu,DeltaX,DeltaT,ordre0)
    P=1+2*B+C
    Q=-(B+C/2)    
    return P,Q



def funcRS(modulation,frontiere,milieu,DeltaX,DeltaT,ordre0="no"):
    A,B,C=coefsABC(modulation,frontiere,milieu,DeltaX,DeltaT,ordre0)
    Rn=2*(1-A**2+2*B)
    Sn=A**2-2*B
    return Rn,Sn


def funcPQ0MK(modulation,frontiere,milieu,DeltaX,DeltaT):
    rho0,E0=mat00(modulation,frontiere,milieu)
    P=rho0
    Q=0  
    return P,Q

def funcRS0MKE(modulation,frontiere,milieu,DeltaX,DeltaT):
    rho0,E0=mat00(modulation,frontiere,milieu)
    A=E0*(DeltaT/DeltaX)**2
    Rn=-2*A
    Sn=A
    return Rn,Sn

def funcRS0MKrhom(modulation,frontiere,milieu,DeltaX,DeltaT):
    rho0,E0=mat00(modulation,frontiere,milieu)
    Rn=rho0
    Sn=0
    return Rn,Sn

def funcRS0MKrhop(modulation,frontiere,milieu,DeltaX,DeltaT):
    rho0,E0=mat00(modulation,frontiere,milieu)
    Rn=rho0
    Sn=0
    return Rn,Sn

def funcRSm10MK(modulation,frontiere,milieu,DeltaX,DeltaT):
    rho0,E0=mat00(modulation,frontiere,milieu)
    Rm1=-rho0
    Sm1=0
    return Rm1,Sm1

def funcRSm1(modulation,frontiere,milieu,DeltaX,DeltaT,ordre0="no"):
    A,B,C=coefsABC(modulation,frontiere,milieu,DeltaX,DeltaT,ordre0)
    Rm1=-(1+2*B-C)
    Sm1=B-C/2
    return Rm1,Sm1

def PfuncPQ0MK(frontiere,npbCell,milieu,DeltaX,DeltaT):
    rho0,E0=mat00Period(frontiere,npbCell,milieu)
    P=rho0
    Q=0  
    return P,Q

def PfuncRS0MKE(frontiere,npbCell,milieu,DeltaX,DeltaT):
    rho0,E0=mat00Period(frontiere,npbCell,milieu)
    A=E0*(DeltaT/DeltaX)**2
    Rn=-2*A
    Sn=A
    return Rn,Sn

def PfuncRS0MKrhom(frontiere,npbCell,milieu,DeltaX,DeltaT):
    rho0,E0=mat00Period(frontiere,npbCell,milieu)
    Rn=rho0
    Sn=0
    return Rn,Sn

def PfuncRS0MKrhop(frontiere,npbCell,milieu,DeltaX,DeltaT):
    rho0,E0=mat00Period(frontiere,npbCell,milieu)
    Rn=rho0
    Sn=0
    return Rn,Sn

def PfuncRSm10MK(frontiere,npbCell,milieu,DeltaX,DeltaT):
    rho0,E0=mat00Period(frontiere,npbCell,milieu)
    Rm1=-rho0
    Sm1=0
    return Rm1,Sm1
    

def matrixRS(R,S,Nx):
    matM=np.zeros([Nx+1,Nx+1])
    for i in range(Nx):
        matM[i,i]=R
        if i==0:
            matM[i+1,i]=S
        elif i==Nx:
            matM[i-1,i]=S
        else:
            matM[i+1,i]=S
            matM[i-1,i]=S
    return matM


  
def matrixRSm1(R,S,Nx):
    matN=np.zeros([Nx+1,Nx+1])
    for i in range(Nx):
        matN[i,i]=R
        if i==0:
            matN[i+1,i]=S
        elif i==Nx:
            matN[i-1,i]=S
        else:
            matN[i+1,i]=S
            matN[i-1,i]=S
    return matN    
 


def dirac_source(x, x_0, dx):
    if abs(x-x_0) < dx/2:
        return 1.0/dx
    elif abs(x-x_0) == dx/2:
        return 1.0/dx/2
    else:
        return 0.0  
    
def matF(sourcef,Nt,DeltaT,Nx,DeltaX):
    x_0=float(sourcef["Xsource"])
    vectDirac=np.zeros(Nx+1)
    for i in range(Nx+1):
        vectDirac[i]=dirac_source(i*DeltaX,x_0,DeltaX)
    matFS=np.zeros([Nt+1,Nx+1])
    for j in range(Nt+1):
        f=source.choice_timefct(sourcef,j*DeltaT)
        matFS[j,:]=vectDirac*f*DeltaT**2
    return matFS 



def vectD(Un,Unm1,F,N,Nm1):
    D=np.dot(N,Un)+np.dot(Nm1,Unm1)+F
    return D
    



def matrixImplHomog(P,Q,Nx):
    matA=np.zeros([Nx+1,Nx+1])
    for i in range(Nx+1):
        matA[i,i]=P
        if i==0:
            matA[i,i+1]=Q
        elif i==Nx:
            matA[i,i-1]=Q
        else:
            matA[i,i+1]=Q
            matA[i,i-1]=Q
    return matA
        

    
def algorithmThomas(A,D):
    ## Solving A x = D with A tridiagonal
    Nx=A.shape[0]
    diag=np.zeros(Nx)
    diaginf=np.zeros(Nx)
    diagsup=np.zeros(Nx)
    for i in range(Nx):
        diag[i]=A[i,i]
        if i<Nx-1:
            diagsup[i]=A[i,i+1]
        if i>=1:
            diaginf[i]=A[i,i-1]
    p=np.zeros(Nx)
    q=np.zeros(Nx)
    p[0]=diagsup[0]/diag[0]
    q[0]=D[0]/diag[0]
    for i in range(1,Nx):
        den=+diag[i]-diaginf[i]*p[i-1]
        p[i]=diagsup[i]/den
        q[i]=(D[i]-diaginf[i]*q[i-1])/den
    X=np.zeros(Nx)
    X[Nx-1]=q[Nx-1]
    for i in range(Nx-2,-1,-1):
        X[i]=q[i]-p[i]*X[i+1]
    return X
        
def P1P2(Nx,DeltaX,time,modulation,frontiere,milieu):
    t=sp.symbols("t")
    CDP=CPP1(modulation,frontiere,milieu)
    dCDP=sp.diff(CDP,t)
    CPP=sp.lambdify(t,CDP)
    dCPP=sp.lambdify(t,dCDP)
    alpha0=float(frontiere["Alpha_0"])
    h=float(frontiere["Alpha_1"])-alpha0
    hP1=np.zeros(Nx+1)
    hP2=np.zeros(Nx+1)
    hdP1=np.zeros(Nx+1)
    hdP2=np.zeros(Nx+1)
    valCPP=CPP(time)
    valdCPP=dCPP(time)
    for i in range(Nx+1):
        Y=(i*DeltaX+alpha0)%h
        y=Y/h
        hP1[i]=valCPP*(y-1/2)*h
        hP2[i]=-1/12*valCPP*(6*y**2-6*y+1)*h**2
        hdP1[i]=valdCPP*(y-1/2)*h
        hdP2[i]=-1/12*valdCPP*(6*y**2-6*y+1)*h**2
    return hP1,hP2,hdP1,hdP2



def temporalSchemePS(Nt,DeltaT,config,source,modulation,frontiere,milieu,ordre0="no"):
    Nx=int(config["Nx"])
    L=float(config["Xsup"])-float(config["Xinf"])
    DeltaX=L/Nx
    DeltaTT=deltat(DeltaX,modulation,frontiere,milieu)
    sources=matF(source,Nt,DeltaT,Nx,DeltaX)
    U=np.zeros([Nt+1,Nx+1])
    P1=np.zeros([Nt+1,Nx+1])
    P2=np.zeros([Nt+1,Nx+1])
    dP1=np.zeros([Nt+1,Nx+1])
    dP2=np.zeros([Nt+1,Nx+1])
    Unm1=np.zeros(Nx+1)
    Un=np.zeros(Nx+1)
    U[0,:]=Un
    t=np.linspace(0,Nt*DeltaT,Nt)
    P,Q=funcPQ(modulation,frontiere,milieu,DeltaX,DeltaT)
    R,S=funcRS(modulation,frontiere,milieu,DeltaX,DeltaT)
    Rm1,Sm1=funcRSm1(modulation,frontiere,milieu,DeltaX,DeltaT)
    t=sp.symbols('t')
    funcP=sp.lambdify(t,P,modules="numpy")
    funcQ=sp.lambdify(t,Q,modules="numpy")
    funcR=sp.lambdify(t,R,modules="numpy")
    funcS=sp.lambdify(t,S,modules="numpy")
    funcRm1=sp.lambdify(t,Rm1,modules="numpy")
    funcSm1=sp.lambdify(t,Sm1,modules="numpy")
    if ordre0=="yes":
        U0=np.zeros([Nt+1,Nx+1])
        Un0=np.zeros(Nx+1)
        Unm10=np.zeros(Nx+1)
        U0[0,:]=Un0
        P0,Q0=funcPQ(modulation,frontiere,milieu,DeltaX,DeltaT,ordre0)
        R0,S0=funcRS(modulation,frontiere,milieu,DeltaX,DeltaT,ordre0)
        Rm10,Sm10=funcRSm1(modulation,frontiere,milieu,DeltaX,DeltaT,ordre0)
        t=sp.symbols('t')
        funcP0=sp.lambdify(t,P0,modules="numpy")
        funcQ0=sp.lambdify(t,Q0,modules="numpy")
        funcR0=sp.lambdify(t,R0,modules="numpy")
        funcS0=sp.lambdify(t,S0,modules="numpy")
        funcRm10=sp.lambdify(t,Rm10,modules="numpy")
        funcSm10=sp.lambdify(t,Sm10,modules="numpy")
    for it in range(0,Nt):
        print(it)
        time=(it)*DeltaT
        A=matrixImplHomog(funcP(time),funcQ(time),Nx)
        N=matrixRS(funcR(time),funcS(time),Nx)
        Nm1=matrixRSm1(funcRm1(time),funcSm1(time),Nx)
        D=vectD(Un,Unm1,sources[it,:],N,Nm1)
        Unp1=algorithmThomas(A, D)
        U[it+1,:]=Unp1
        Unm1=Un
        Un=Unp1
        if ordre0=="yes":
            A0=matrixImplHomog(funcP0(time),funcQ0(time),Nx)
            N0=matrixRS(funcR0(time),funcS0(time),Nx)
            Nm10=matrixRSm1(funcRm10(time),funcSm10(time),Nx)
            D0=vectD(Un0,Unm10,sources[it,:],N0,Nm10)
            Unp10=algorithmThomas(A0, D0)
            U0[it+1,:]=Unp10
            Unm10=Un0
            Un0=Unp10
        P1[it,:],P2[it,:],dP1[it,:],dP2[it,:]=P1P2(Nx,DeltaX,time,modulation,frontiere,milieu)
    if ordre0=="yes":
        return U0,U,DeltaTT,P1,P2,dP1,dP2
    else:
        return U,DeltaTT,P1,P2

def temporalSchemePS0MK(Nt,DeltaT,config,source,modulation,frontiere,milieu):
    Nx=int(config["Nx"])
    L=float(config["Xsup"])-float(config["Xinf"])
    rhoinit=float(milieu["Rho_1"])
    DeltaX=L/Nx
    DeltaTT=deltat(DeltaX,modulation,frontiere,milieu)
    DeltaTT=min(DeltaT,DeltaTT)
    sources=matF(source,Nt,DeltaT,Nx,DeltaX)
    rho0,E0=mat00(modulation,frontiere,milieu)
    U=np.zeros([Nt+1,Nx+1])
    Unm1=np.zeros(Nx+1)
    Un=np.zeros(Nx+1)
    U[0,:]=Un
    t=np.linspace(0,Nt*DeltaT,Nt)
    P,Q=funcPQ0MK(modulation,frontiere,milieu,DeltaX,DeltaT)
    R,S=funcRS0MKE(modulation,frontiere,milieu,DeltaX,DeltaT)
    R2m,S2m=funcRS0MKrhom(modulation,frontiere,milieu,DeltaX,DeltaT)
    R2p,S2p=funcRS0MKrhom(modulation,frontiere,milieu,DeltaX,DeltaT)
    Rm1,Sm1=funcRSm10MK(modulation,frontiere,milieu,DeltaX,DeltaT)
    t=sp.symbols('t')
    funcP=sp.lambdify(t,P,modules="numpy")
    funcQ=sp.lambdify(t,Q,modules="numpy")
    funcR=sp.lambdify(t,R,modules="numpy")
    funcS=sp.lambdify(t,S,modules="numpy")
    funcR2m=sp.lambdify(t,R2m,modules="numpy")
    funcS2m=sp.lambdify(t,S2m,modules="numpy")
    funcR2p=sp.lambdify(t,R2p,modules="numpy")
    funcS2p=sp.lambdify(t,S2p,modules="numpy")
    funcRm1=sp.lambdify(t,Rm1,modules="numpy")
    funcSm1=sp.lambdify(t,Sm1,modules="numpy")
    for it in range(0,Nt):
        print(it)
        time=(it)*DeltaT
        A=matrixImplHomog(funcP(time+DeltaTT/2),funcQ(time+DeltaTT/2),Nx)
        N=matrixRS(funcR(time),funcS(time),Nx)
        N2p=matrixRS(funcR2p(time+DeltaTT/2),funcS2p(time+DeltaTT/2),Nx)
        N2m=matrixRS(funcR2m(time-DeltaTT/2),funcS2m(time-DeltaTT/2),Nx)
        Nm1=matrixRSm1(funcRm1(time-DeltaTT/2),funcSm1(time-DeltaTT/2),Nx)
        D=vectD(Un,Unm1,funcP(time+DeltaTT/2)/funcP(time)*rhoinit*sources[it,:],N+N2p+N2m,Nm1)
        Unp1=np.linalg.solve(A, D)#algorithmThomas(A, D)
        U[it+1,:]=Unp1
        Unm1=Un
        Un=Unp1
    return U,DeltaTT  

def temporalSchemePS0MKP(Nt,DeltaT,config,source,npbCell,frontiere,milieu):
    Nx=int(config["Nx"])
    L=float(config["Xsup"])-float(config["Xinf"])
    rhoinit=float(milieu["Rho_1"])
    DeltaX=L/Nx
    DeltaTT=deltatp(DeltaX,frontiere,npbCell,milieu)
    DeltaTT=min(DeltaT,DeltaTT)
    sources=matF(source,Nt,DeltaT,Nx,DeltaX)
    rho0,E0=mat00Period(frontiere,npbCell,milieu)
    print(rho0,E0)
    U=np.zeros([Nt+1,Nx+1])
    Unm1=np.zeros(Nx+1)
    Un=np.zeros(Nx+1)
    U[0,:]=Un
    t=np.linspace(0,Nt*DeltaT,Nt)
    P,Q=PfuncPQ0MK(frontiere,npbCell,milieu,DeltaX,DeltaT)
    R,S=PfuncRS0MKE(frontiere,npbCell,milieu,DeltaX,DeltaT)
    R2m,S2m=PfuncRS0MKrhom(frontiere,npbCell,milieu,DeltaX,DeltaT)
    R2p,S2p=PfuncRS0MKrhom(frontiere,npbCell,milieu,DeltaX,DeltaT)
    Rm1,Sm1=PfuncRSm10MK(frontiere,npbCell,milieu,DeltaX,DeltaT)
    t=sp.symbols('t')
    funcP=sp.lambdify(t,P,modules="numpy")
    funcQ=sp.lambdify(t,Q,modules="numpy")
    funcR=sp.lambdify(t,R,modules="numpy")
    funcS=sp.lambdify(t,S,modules="numpy")
    funcR2m=sp.lambdify(t,R2m,modules="numpy")
    funcS2m=sp.lambdify(t,S2m,modules="numpy")
    funcR2p=sp.lambdify(t,R2p,modules="numpy")
    funcS2p=sp.lambdify(t,S2p,modules="numpy")
    funcRm1=sp.lambdify(t,Rm1,modules="numpy")
    funcSm1=sp.lambdify(t,Sm1,modules="numpy")
    for it in range(0,Nt):
        print(it)
        time=(it)*DeltaT
        A=matrixImplHomog(funcP(time+DeltaTT/2),funcQ(time+DeltaTT/2),Nx)
        N=matrixRS(funcR(time),funcS(time),Nx)
        N2p=matrixRS(funcR2p(time+DeltaTT/2),funcS2p(time+DeltaTT/2),Nx)
        N2m=matrixRS(funcR2m(time-DeltaTT/2),funcS2m(time-DeltaTT/2),Nx)
        Nm1=matrixRSm1(funcRm1(time-DeltaTT/2),funcSm1(time-DeltaTT/2),Nx)
        D=vectD(Un,Unm1,funcP(time+DeltaTT/2)/funcP(time)*rhoinit*sources[it,:],N+N2p+N2m,Nm1)
        Unp1=np.linalg.solve(A, D)#algorithmThomas(A, D)
        U[it+1,:]=Unp1
        Unm1=Un
        Un=Unp1
    return U,DeltaTT  

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

