#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 17:09:49 2023

@author: michael
"""

import numpy as np
import sympy as sp

def hPropBicouche(cm1,rho1,cm2,rho2,l1,l2,CFL):
    L=l1+l2
    E1=rho1*cm1**2
    E2=rho2*cm2**2
    rhoh=(l1*rho1+l2*rho2)/L
    Eh=L/(l1/E1+l2/E2)
    cmh=np.sqrt(Eh/rhoh)
    CFLh=CFL*cmh/max(cm1,cm2)
    return cmh,rhoh,CFLh

def hPropCracks(cm,rho,K,M,l):
    E=rho*cm**2
    rhoh=rho+M/l
    Eh=K*l/(1+K*l/E)
    cmh=np.sqrt(Eh/rhoh)
    fCFLh=cmh/cm
    return cmh,rhoh,fCFLh

def propHomo(t,rho,cm,frontiere,modulation):
    contacti="Contact_"+str(1)
    E=rho*cm**2
    alpha=float(frontiere["Alpha_1"])
    h=float(frontiere["Alpha_2"])-float(frontiere["Alpha_1"])
    if frontiere[contacti]=="parfait":
        rhoeff=rho
        Eeff=E
        Ceff=np.sqrt(Eeff/rhoeff)
    if frontiere[contacti]=="masse-ressort":
        nameKni="Kn_"+str(1)
        nameMni="Mn_"+str(1)
        K=float(frontiere[nameKni])
        M=float(frontiere[nameMni])
        rhoeff=rho+M/h
        Eeff=(K*h)/(1+K*h/E)
        Ceff=np.sqrt(Eeff/rhoeff)
    if frontiere[contacti]=="masse-ressort module":
        ampF=float(modulation["DeltaF"])
        ampM=float(modulation["DeltaM"])
        if modulation["Synchro"]=="Yes":
            FreqM=float(modulation["Freq"])
            FreqF=float(modulation["Freq"])
            celM=float(modulation["ModC"])
            celF=float(modulation["ModC"])
        else:
            FreqM=float(modulation["FreqM"])
            FreqF=float(modulation["FreqF"])
            celM=float(modulation["ModM"])
            celF=float(modulation["ModF"])
        if celM==0:
            phiMm=0
        else:
            phiMm=FreqM*2*np.pi*alpha/celM
        if celF==0:
            phiFm=0
        else:
            phiFm=FreqF*2*np.pi*alpha/celF
        if modulation["TypeMod"]=="Sinus":
            fctM=np.sin(2*np.pi*FreqM*t+phiMm)
            fctF=np.sin(2*np.pi*FreqF*t+phiFm)
        elif modulation["TypeMod"]=="Square":
            fctM=np.sign(np.sin(2*np.pi*FreqM*t+phiMm))
            fctF=np.sign(np.sin(2*np.pi*FreqF*t+phiFm))
        elif modulation["TypeMod"]=="Square NS":
            ratio=float(modulation["Ratio"])
            fctM=-np.sign(FreqM*t-int(FreqM*t)-ratio)
            fctF=-np.sign(FreqF*t-int(FreqF*t)-ratio)
        nameKni="Kn_"+str(1)
        nameMni="Mn_"+str(1)
        K=float(frontiere[nameKni])
        M=float(frontiere[nameMni])
        if K==0:
            Fmod=0
        else:
            Fmod=(1+ampF*fctF)/K
        Mmod=(1+ampM*fctM)*M
        rhoeff=rho+Mmod/h
        Eeff=(h)/(Fmod+h/E)
        Ceff=np.sqrt(Eeff/rhoeff)
    elif frontiere[contacti]=="masse-ressort module 2":
        ampF=-float(modulation["DeltaF"])
        ampM=-float(modulation["DeltaM"])
        if modulation["Synchro"]=="Yes":
            FreqM=float(modulation["Freq"])
            FreqF=float(modulation["Freq"])
            celM=float(modulation["ModC"])
            celF=float(modulation["ModC"])
        else:
            FreqM=float(modulation["FreqM"])
            FreqF=float(modulation["FreqF"])
            celM=float(modulation["ModM"])
            celF=float(modulation["ModF"])
        if celM==0:
            phiMm=0
        else:
            phiMm=FreqM*2*np.pi*alpha/celM
        if celF==0:
            phiFm=0
        else:
            phiFm=FreqF*2*np.pi*alpha/celF
        if modulation["TypeMod"]=="Sinus":
            fctM=np.sin(2*np.pi*FreqM*t+phiMm)
            fctF=np.sin(2*np.pi*FreqF*t+phiFm)
        elif modulation["TypeMod"]=="Square":
            fctM=np.sign(np.sin(2*np.pi*FreqM*t+phiMm))
            fctF=np.sign(np.sin(2*np.pi*FreqF*t+phiFm))
        elif modulation["TypeMod"]=="Square NS":
            ratio=float(modulation["Ratio"])
            fctM=-np.sign(FreqM*t-int(FreqM*t)-ratio)
            fctF=-np.sign(FreqF*t-int(FreqF*t)-ratio)
        nameKni="Kn_"+str(1)
        nameMni="Mn_"+str(1)
        K=float(frontiere[nameKni])
        M=float(frontiere[nameMni])
        if K==0:
            Fmod=0
        else:
            Fmod=(1+ampF*fctF)/K
        Mmod=(1+ampM*fctM)*M
        rhoeff=rho+Mmod/h
        Eeff=(h)/(Fmod+h/E)
        Ceff=np.sqrt(Eeff/rhoeff)
    if frontiere[contacti]=="masse-ressort module dephase":
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
                fct.append(np.sin(2*np.pi*(Freq*t+phim*(alpha+i*h))))
        elif modulation["TypeMod"]=="Square":
            for i in range(Nip):
                fct.append(np.sign(np.sin(2*np.pi*(Freq*t+phim*(alpha+i*h)))))
        elif modulation["TypeMod"]=="Square NS":
            ratio=float(modulation["Ratio"])
            sigph=Freq*t-phim*(alpha+i*h)%1
            for i in range(Nip):
                fct.append(np.where(sigph<ratio,1,-1))
        K=[]
        M=[]
        Mmod=0
        Fmod=0
        for i in range(Nip):
            nameKni="Kn_"+str(i)
            nameMni="Mn_"+str(i)
            K.append(float(frontiere[nameKni]))
            M.append(float(frontiere[nameMni]))
            if K[i]==0:
                Fmod=Fmod
            else:
                Fmod=Fmod+1/K[i]*(1+ampF*fct[i])
            Mmod=Mmod+M[i]*(1+ampM*fct[i])
        rhoeff=rho+Mmod/Nip*h
        Eeff=(Nip*h)/(Fmod+Nip*h/E)
        Ceff=np.sqrt(Eeff/rhoeff)
    return rhoeff,Ceff

def propHomoFormP(Xmin,milieu,frontiere,nbpCell):
    t=sp.symbols("t")
    Mmod=0
    Fmod=0
    Bmod=0
    Zmod=0
    rhomoy=0
    Emoym1=0
    periodFront="Alpha_"+str(nbpCell-1)
    h=float(frontiere[periodFront])
    for i in range(0,nbpCell):
        contacti="Contact_"+str(i)
        alpha= float(frontiere["Alpha_"+str(i)])
        if i==0:
            li=alpha-Xmin
        else:
            li=alpha-float(frontiere["Alpha_"+str(i-1)])
        rhoi=float(milieu["Rho_"+str(i)])
        celi=float(milieu["Cel_"+str(i)])
        Ei=rhoi*celi**2
        if frontiere[contacti]=="masse-ressort module independant":
            ampF=float(frontiere["DeltaF_"+str(i)])
            ampM=float(frontiere["DeltaM_"+str(i)])
            if frontiere["Synchro_"+str(i)]=="Yes":
                FreqM=float(frontiere["Freq_"+str(i)])
                FreqF=float(frontiere["Freq_"+str(i)])
                FreqB=float(frontiere["Freq_"+str(i)])
                FreqZ=float(frontiere["Freq_"+str(i)])
                phiF0=float(frontiere["phiF0_"+str(i)])
                phiM0=float(frontiere["phiM0_"+str(i)])
                nfrB = "beta_"+str(i)
                if frontiere.get(nfrB)==None:
                    ampB=0.
                    phiB0=0.
                else:
                    ampB= float(frontiere["DeltaB_"+str(i)])
                    phiB0=float(frontiere["phiB0_"+str(i)])
                nfrZ = "zeta_"+str(i)
                if frontiere.get(nfrZ)==None:
                    ampZ=0.
                    phiZ0=0.
                else:
                    ampZ= float(frontiere["DeltaZ_"+str(i)])
                    phiZ0=float(frontiere["phiZ0_"+str(i)])
            else:
                FreqM=float(frontiere["FreqM_"+str(i)])
                FreqF=float(frontiere["FreqF_"+str(i)])
                nfrB = "beta_"+str(i)
                if frontiere.get(nfrB)==None:
                    FreqB=0.
                    ampB=0.
                    phiB0=0.
                else:
                    FreqB= float(frontiere["FreqB_"+str(i)])
                    ampB= float(frontiere["DeltaB_"+str(i)])
                    phiB0=float(frontiere["phiB0_"+str(i)])
                nfrZ = "zeta_"+str(i)
                if frontiere.get(nfrZ)==None:
                    FreqZ=0.
                    ampZ=0.
                    phiZ0=0.
                else:
                    FreqZ= float(frontiere["FreqZ_"+str(i)])
                    ampZ= float(frontiere["DeltaZ_"+str(i)])
                    phiZ0=float(frontiere["phiZ0_"+str(i)])
                phiF0=float(frontiere["phiF0_"+str(i)])
                phiM0=float(frontiere["phiM0_"+str(i)])
            phiMm=2*np.pi*phiM0
            phiFm=2*np.pi*phiF0
            phiBm=2*np.pi*phiB0
            phiZm=2*np.pi*phiZ0
            if frontiere["TypeMod_"+str(i)]=="Sinus":
                fctM=sp.sin(2*np.pi*FreqM*t+phiMm)
                fctF=sp.sin(2*np.pi*FreqF*t+phiFm)
                fctB=sp.sin(2*np.pi*FreqB*t+phiBm)
                fctZ=sp.sin(2*np.pi*FreqZ*t+phiZm)
            elif frontiere["TypeMod_"+str(i)]=="Cosinus":
                fctM=sp.cos(2*np.pi*FreqM*t+phiMm)
                fctF=sp.cos(2*np.pi*FreqF*t+phiFm)
                fctB=sp.cos(2*np.pi*FreqB*t+phiBm)
                fctZ=sp.cos(2*np.pi*FreqZ*t+phiZm)
            elif frontiere["TypeMod_"+str(i)]=="Square":
                fctM=sp.sign(sp.sin((2*np.pi*FreqM*t+phiMm)))
                fctF=sp.sign(sp.sin((2*np.pi*FreqF*t+phiFm)))
                fctB=sp.sign(sp.sin((2*np.pi*FreqB*t+phiBm)))
                fctZ=sp.sign(sp.sin((2*np.pi*FreqZ*t+phiZm)))
            elif frontiere["TypeMod_"+str(i)]=="Square NS":
                ratio=float(frontiere["Ratio"+str(i)])
                fctM=-sp.sign(FreqM*t-sp.floor(FreqM*t)-ratio)
                fctF=-sp.sign(FreqF*t-sp.floor(FreqF*t)-ratio)
                fctB=-sp.sign(FreqB*t-sp.floor(FreqB*t)-ratio)
                fctZ=-sp.sign(FreqZ*t-sp.floor(FreqZ*t)-ratio)
            nameKni="Kn_"+str(i)
            nameMni="Mn_"+str(i)
            K=float(frontiere[nameKni])
            if K==0:
                Fmodi=0
            else:
                Fmodi=(1+ampF*fctF)/K
            M=float(frontiere[nameMni])
            namorK = "beta_"+str(i)
            if frontiere.get(namorK)==None:
                etaK=0.
            else:
                etaK= float(frontiere[namorK])
            namorM = "zeta_"+str(i)
            if frontiere.get(namorM)==None:
                etaM=0.
            else:
                etaM= float(frontiere[namorM])
            Mmodi=(1+ampM*fctM)*M
            Bmodi=(1+ampB*fctB)*etaK
            Zmodi=(1+ampZ*fctZ)*etaM
            Mmod=Mmod+Mmodi
            Fmod=Fmod+Fmodi
            Bmod=Bmod+Bmodi
            Zmod=Zmod+Zmodi
        rhomoy=rhomoy+rhoi*li/h
        Emoym1=Emoym1+li/h/Ei
    rhoeff=rhomoy+Mmod/h
    if Fmod==0:
        Eeff=1/(Emoym1)
    else:
        unSurEeff=Emoym1+Fmod/h
        Eeff=1/(unSurEeff)
    Ceff=sp.sqrt(Eeff/rhoeff)
    QCeff=Bmod/h
    QMeff=Zmod/h
    return rhoeff,Eeff,Ceff,QCeff,QMeff
    
def propHomoForm(rho,cm,frontiere,modulation):
    t=sp.symbols("t")
    contacti="Contact_"+str(1)
    E=rho*cm**2
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
        rhoeff=rho+Mmod/h
        Eeff=(h)/(Fmod+h/E)
        Ceff=sp.sqrt(Eeff/rhoeff)
        QCeff=Bmod/h
        QMeff=Zmod/h
    if frontiere[contacti]=="masse-ressort module independant":
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
        rhoeff=rho+Mmod/h
        Eeff=(h)/(Fmod+h/E)
        Ceff=sp.sqrt(Eeff/rhoeff)
        QCeff=Bmod/h
        QMeff=Zmod/h
    if frontiere[contacti]=="masse-ressort module 2":
        ampF=-float(modulation["DeltaF"])
        ampM=-float(modulation["DeltaM"])
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
        rhoeff=rho+Mmod/h
        Eeff=(h)/(Fmod+h/E)
        Ceff=sp.sqrt(Eeff/rhoeff)
        QCeff=Bmod/h
        QMeff=Zmod/h
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
        rhoeff=rho+Mmod/Nip/h
        Eeff=1/(Fmod/Nip/h+1/E)
        Ceff=sp.sqrt(Eeff/rhoeff)
        QCeff=Bmod/alpha
        QMeff=Zmod/alpha
    return rhoeff,Eeff,Ceff,QCeff,QMeff


def updateCFL(Nt,dt,t,rho,cm,frontiere,modulation):
    fCFL=1
    if frontiere["Contact_"+str(0)]=="masse-ressort":
        rhoeff,Ceff=propHomo(0,rho[0],cm[0],frontiere,modulation)
        fCFL=Ceff/cm[0]
    for ts in range(Nt):
        if frontiere["Contact_"+str(0)]=="masse-ressort module":
            rhoeff,Ceff=propHomo(ts*dt,rho[0],cm[0],frontiere,modulation)
            fCFL=max(Ceff/cm[0], fCFL)
        if frontiere["Contact_"+str(0)]=="masse-ressort module dephase":
            rhoeff,Ceff=propHomo(ts*dt,rho[0],cm[0],frontiere,modulation)
            fCFL=max(Ceff/cm[0], fCFL)
    return fCFL

def updateCFLP(Nt,dt,t,Xmin,nbpCell,frontiere,milieu):
    fCFL=1
    rhoeff,Eeff,Ceff,QCeff,QMeff=propHomoFormP(Xmin,milieu,frontiere,nbpCell)
    cmedmax=0
    for i in range(nbpCell):
        cmed=float(milieu["Cel_"+str(i)])
        cmedmax=max(cmedmax,cmed)
    t = sp.symbols('t')
    cel=sp.lambdify(t,Ceff)
    maxcel=0
    for ts in range(Nt):
        maxcel=max(maxcel,cel(ts*dt))
    fCFL=max(maxcel/cmedmax, fCFL)
    return fCFL


def timeProp(Nt,dt,t,rho,cm,frontiere,modulation):
    rhoeffT=np.zeros(Nt)
    CeffT=np.zeros(Nt)
    if frontiere["Contact_"+str(0)]=="masse-ressort":
        rhoeffT,CeffT=propHomo(0,rho[0],cm[0],frontiere,modulation)
    for ts in range(Nt):
        if frontiere["Contact_"+str(0)]=="masse-ressort module":
            rhoeffT[ts],CeffT[ts]=propHomo(ts*dt,rho[0],cm[0],frontiere,modulation)
        if frontiere["Contact_"+str(0)]=="masse-ressort module dephase":
            rhoeffT[ts],CeffT[ts]=propHomo(ts*dt,rho[0],cm[0],frontiere,modulation)
    return rhoeffT,CeffT






