#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 11:10:30 2023

@author: michael
"""

import numpy as np



def fonctionModK(t,x,modulation,frontiere,i): 
    contacti="Contact_"+str(i)
    Delta=float(modulation["DeltaF"])
    if frontiere[contacti]=="masse-ressort module":
        k=0
        if modulation["Synchro"]=="Yes":
            Omega=2*np.pi*float(modulation["Freq"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phi0"])
            else:
                phi0=0    
        else:
            Omega=2*np.pi*float(modulation["FreqF"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phiF0"])
            else:
                phi0=0
    if frontiere[contacti]=="masse-ressort module independant":
        k=0
        if frontiere["Synchro"+str(i)]=="Yes":
            Omega=2*np.pi*float(frontiere["Freq"+str(i)])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phiF0"])
            else:
                phi0=0    
        else:
            Omega=2*np.pi*float(modulation["FreqF"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phiF0"])
            else:
                phi0=0
    if frontiere[contacti]=="masse-ressort module key":
        k=0
        if modulation["Synchro"]=="Yes":
            Omega=2*np.pi*float(modulation["Freq"])
            cel=2800
            DelX=float(frontiere["Alpha_1"])-float(frontiere["Alpha_0"])
            tau=DelX/cel
            phi0=-Omega*tau
            Delta=-Delta
    elif frontiere[contacti]=="masse-ressort module dephase":
        if modulation["Synchro"]=="Yes":
            Omega=2*np.pi*float(modulation["Freq"])
            if modulation["Dephasage"]=="xct":
                celF=float(modulation["ModC"])
                if celF==0:
                    Omega=0
                    k=0
                else:
                    k=Omega/celF
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phi0"])
            else:
                phi0=0    
        else:
            Omega=-2*np.pi*float(modulation["FreqF"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phiF0"])
            else:
                phi0=0  
            if modulation["Dephasage"]=="xct":
                celF=float(modulation["ModF"])
                if celF==0:
                    Omega=0
                    k=0
                else:
                    k=Omega/celF   
            else:
                k=0
    if modulation["TypeMod"]=="Sinus":
        fct=Delta*np.sin(Omega*t-k*x+phi0)
    if modulation["TypeMod"]=="Cosinus":
        fct=Delta*np.cos(Omega*t-k*x+phi0)
    if modulation["TypeMod"]=="Cosinus QP":
        fct=Delta*(0.4*np.cos(Omega*t-k*x)+0.13*np.cos(Omega/np.pi*t-k*x)+0.25*np.cos(Omega*np.sqrt(2)*t-k*x)+0.22*np.cos(Omega*np.sqrt(3)*t-k*x))
    if modulation["TypeMod"]=="QP-f":
        fct=Delta/2*(np.sin(Omega*t-k*x+phi0)+np.sin(4*Omega/np.pi*t-k*x))
    if modulation["TypeMod"]=="QP-r":
        fct=Delta/2*(np.sin(Omega*t-k*x+phi0)+np.sin(Omega*np.sqrt(2)*t-k*x))
    elif modulation["TypeMod"]=="Square":
        fct=Delta*np.sign(np.sin(Omega*t-k*x+phi0))
    elif modulation["TypeMod"]=="SquareSym":
        fct=Delta*np.sign(np.cos(Omega*t-k*x+phi0))
    elif modulation["TypeMod"]=="Square NS":
        ratio=float(modulation["Ratio"])
        fct=-Delta*np.sign(Omega/2/np.pi*t-np.floor(Omega/2/np.pi*t)-ratio)
    return fct

def fonctionModM(t,x,modulation,frontiere,i): 
    contacti="Contact_"+str(i)
    Delta=float(modulation["DeltaM"])
    if frontiere[contacti]=="masse-ressort module":
        k=0
        if modulation["Synchro"]=="Yes":
            Omega=2*np.pi*float(modulation["Freq"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phiM0"])
            else:
                phi0=0    
        else:
            Omega=2*np.pi*float(modulation["FreqM"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phiM0"])
            else:
                phi0=0
    elif frontiere[contacti]=="masse-ressort module dephase":
        if modulation["Synchro"]=="Yes":
            Omega=2*np.pi*float(modulation["Freq"])
            if modulation["Dephasage"]=="xct":
                celM=float(modulation["ModC"])
                if celM==0:
                    Omega=0
                    k=0
                else:
                    k=Omega/celM
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phi0"])
            else:
                phi0=0    
        else:
            Omega=-2*np.pi*float(modulation["FreqM"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phiM0"])
            else:
                phi0=0  
            if modulation["Dephasage"]=="xct":
                celM=float(modulation["ModM"])
                if celM==0:
                    Omega=0
                    k=0
                else:
                    k=Omega/celM   
            else:
                k=0
    if modulation["TypeMod"]=="Sinus":
        fct=Delta*np.sin(Omega*t-k*x+phi0)
    if modulation["TypeMod"]=="Cosinus":
        fct=Delta*np.cos(Omega*t-k*x+phi0)
    if modulation["TypeMod"]=="Cosinus QP":
        fct=Delta*(0.4*np.cos(Omega*t-k*x)+0.13*np.cos(Omega/np.pi*t-k*x)+0.25*np.cos(Omega*np.sqrt(2)*t-k*x)+0.22*np.cos(Omega*np.sqrt(3)*t-k*x))
    if modulation["TypeMod"]=="QP-f":
        fct=Delta/2*(np.sin(Omega*t-k*x+phi0)+np.sin(4*Omega/np.pi*t-k*x))
    if modulation["TypeMod"]=="QP-r":
        fct=Delta/2*(np.sin(Omega*t-k*x+phi0)+np.sin(Omega*np.sqrt(2)*t-k*x))
    elif modulation["TypeMod"]=="Square":
        fct=Delta*np.sign(np.sin(Omega*t-k*x+phi0))
    elif modulation["TypeMod"]=="SquareSym":
        fct=Delta*np.sign(np.cos(Omega*t-k*x+phi0))
    elif modulation["TypeMod"]=="Square NS":
        ratio=float(modulation["Ratio"])
        fct=-Delta*np.sign(Omega/2/np.pi*t-np.floor(Omega/2/np.pi*t)-ratio)
    return fct

def fonctionModKP(t,x,frontiere,i): 
    contacti="Contact_"+str(i)
    if frontiere[contacti]=="masse-ressort module independant":
        k=0
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
            fct=Delta*np.sin(Omega*t-k*x+phi0)
        if frontiere["TypeMod_"+str(i)]=="Cosinus":
            fct=Delta*np.cos(Omega*t-k*x+phi0)
        if frontiere["TypeMod_"+str(i)]=="Cosinus QP":
            fct=Delta*(0.4*np.cos(Omega*t-k*x)+0.13*np.cos(Omega/np.pi*t-k*x)+0.25*np.cos(Omega*np.sqrt(2)*t-k*x)+0.22*np.cos(Omega*np.sqrt(3)*t-k*x))
        if frontiere["TypeMod_"+str(i)]=="QP-f":
            fct=Delta/2*(np.sin(Omega*t-k*x+phi0)+np.sin(4*Omega/np.pi*t-k*x))
        if frontiere["TypeMod_"+str(i)]=="QP-r":
            fct=Delta/2*(np.sin(Omega*t-k*x+phi0)+np.sin(Omega*np.sqrt(2)*t-k*x))
        elif frontiere["TypeMod_"+str(i)]=="Square":
            fct=Delta*np.sign(np.sin(Omega*t-k*x+phi0))
        elif frontiere["TypeMod_"+str(i)]=="SquareSym":
            fct=Delta*np.sign(np.cos(Omega*t-k*x+phi0))
        elif frontiere["TypeMod_"+str(i)]=="Square NS":
            ratio=float(frontiere["Ratio_"+str(i)])
            fct=-Delta*np.sign(Omega/2/np.pi*t-np.floor(Omega/2/np.pi*t)-ratio)
    elif frontiere[contacti]=="parfait":
        fct=0
    return fct

def fonctionModMP(t,x,frontiere,i): 
    contacti="Contact_"+str(i)
    if frontiere[contacti]=="masse-ressort module independant":
        k=0
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
            fct=Delta*np.sin(Omega*t-k*x+phi0)
        if frontiere["TypeMod_"+str(i)]=="Cosinus":
            fct=Delta*np.cos(Omega*t-k*x+phi0)
        if frontiere["TypeMod_"+str(i)]=="Cosinus QP":
            fct=Delta*(0.4*np.cos(Omega*t-k*x)+0.13*np.cos(Omega/np.pi*t-k*x)+0.25*np.cos(Omega*np.sqrt(2)*t-k*x)+0.22*np.cos(Omega*np.sqrt(3)*t-k*x))
        if frontiere["TypeMod_"+str(i)]=="QP-f":
            fct=Delta/2*(np.sin(Omega*t-k*x+phi0)+np.sin(4*Omega/np.pi*t-k*x))
        if frontiere["TypeMod_"+str(i)]=="QP-r":
            fct=Delta/2*(np.sin(Omega*t-k*x+phi0)+np.sin(Omega*np.sqrt(2)*t-k*x))
        elif frontiere["TypeMod_"+str(i)]=="Square":
            fct=Delta*np.sign(np.sin(Omega*t-k*x+phi0))
        elif frontiere["TypeMod_"+str(i)]=="SquareSym":
            fct=Delta*np.sign(np.cos(Omega*t-k*x+phi0))
        elif frontiere["TypeMod_"+str(i)]=="Square NS":
            ratio=float(frontiere["Ratio_"+str(i)])
            fct=-Delta*np.sign(Omega/2/np.pi*t-np.floor(Omega/2/np.pi*t)-ratio)
    elif frontiere[contacti]=="parfait":
        fct=0
    return fct

def fonctionModZ(t,x,modulation,frontiere,i): 
    contacti="Contact_"+str(i)
    Delta=float(modulation["DeltaZ"])
    if frontiere[contacti]=="masse-ressort module":
        k=0
        if modulation["Synchro"]=="Yes":
            Omega=2*np.pi*float(modulation["Freq"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phi0"])
            else:
                phi0=0    
        else:
            Omega=2*np.pi*float(modulation["FreqZ"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phiZ0"])
            else:
                phi0=0
    elif frontiere[contacti]=="masse-ressort module dephase":
        if modulation["Synchro"]=="Yes":
            Omega=2*np.pi*float(modulation["Freq"])
            if modulation["Dephasage"]=="xct":
                celZ=float(modulation["ModC"])
                if celZ==0:
                    Omega=0
                    k=0
                else:
                    k=Omega/celZ
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phi0"])
            else:
                phi0=0    
        else:
            Omega=-2*np.pi*float(modulation["FreqZ"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phiZ0"])
            else:
                phi0=0  
            if modulation["Dephasage"]=="xct":
                celZ=float(modulation["ModZ"])
                if celZ==0:
                    Omega=0
                    k=0
                else:
                    k=Omega/celZ   
            else:
                k=0
    if modulation["TypeMod"]=="Sinus":
        fct=Delta*np.sin(Omega*t-k*x+phi0)
    elif modulation["TypeMod"]=="Cosinus":
        fct=Delta*np.cos(Omega*t-k*x+phi0)
    if modulation["TypeMod"]=="Cosinus QP":
        fct=Delta*(0.4*np.cos(Omega*t-k*x)+0.13*np.cos(Omega/np.pi*t-k*x)+0.25*np.cos(Omega*np.sqrt(2)*t-k*x)+0.22*np.cos(Omega*np.sqrt(3)*t-k*x))
    elif modulation["TypeMod"]=="QP-f":
        fct=Delta/2*(np.sin(Omega*t-k*x+phi0)+np.sin(4*Omega/np.pi*t-k*x))
    elif modulation["TypeMod"]=="QP-r":
        fct=Delta/2*(np.sin(Omega*t-k*x+phi0)+np.sin(Omega*np.sqrt(2)*t-k*x))
    elif modulation["TypeMod"]=="Square":
        fct=Delta*np.sign(np.sin(Omega*t-k*x+phi0))
    elif modulation["TypeMod"]=="SquareSym":
        fct=Delta*np.sign(np.cos(Omega*t-k*x+phi0))
    elif modulation["TypeMod"]=="Square NS":
        ratio=float(modulation["Ratio"])
        fct=-Delta*np.sign(Omega/2/np.pi*t-np.floor(Omega/2/np.pi*t)-ratio)
    return fct

def fonctionModB(t,x,modulation,frontiere,i): 
    contacti="Contact_"+str(i)
    Delta=float(modulation["DeltaB"])
    if frontiere[contacti]=="masse-ressort module":
        k=0
        if modulation["Synchro"]=="Yes":
            Omega=2*np.pi*float(modulation["Freq"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phi0"])
            else:
                phi0=0    
        else:
            Omega=2*np.pi*float(modulation["FreqB"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phiB0"])
            else:
                phi0=0
    elif frontiere[contacti]=="masse-ressort module dephase":
        if modulation["Synchro"]=="Yes":
            Omega=2*np.pi*float(modulation["Freq"])
            if modulation["Dephasage"]=="xct":
                celB=float(modulation["ModC"])
                if celB==0:
                    Omega=0
                    k=0
                else:
                    k=Omega/celB
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phi0"])
            else:
                phi0=0    
        else:
            Omega=-2*np.pi*float(modulation["FreqB"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phiB0"])
            else:
                phi0=0  
            if modulation["Dephasage"]=="xct":
                celB=float(modulation["ModB"])
                if celB==0:
                    Omega=0
                    k=0
                else:
                    k=Omega/celB   
            else:
                k=0
    if modulation["TypeMod"]=="Sinus":
        fct=Delta*np.sin(Omega*t-k*x+phi0)
    if modulation["TypeMod"]=="Cosinus":
        fct=Delta*np.cos(Omega*t-k*x+phi0)
    if modulation["TypeMod"]=="Cosinus QP":
        fct=Delta*(0.4*np.cos(Omega*t-k*x)+0.13*np.cos(Omega/np.pi*t-k*x)+0.25*np.cos(Omega*np.sqrt(2)*t-k*x)+0.22*np.cos(Omega*np.sqrt(3)*t-k*x))
    if modulation["TypeMod"]=="QP-f":
        fct=Delta/2*(np.sin(Omega*t-k*x+phi0)+np.sin(4*Omega/np.pi*t-k*x))
    if modulation["TypeMod"]=="QP-r":
        fct=Delta/2*(np.sin(Omega*t-k*x+phi0)+np.sin(Omega*np.sqrt(2)*t-k*x))
    elif modulation["TypeMod"]=="Square":
        fct=Delta*np.sign(np.sin(Omega*t-k*x+phi0))
    elif modulation["TypeMod"]=="SquareSym":
        fct=Delta*np.sign(np.cos(Omega*t-k*x+phi0))
    elif modulation["TypeMod"]=="Square NS":
        ratio=float(modulation["Ratio"])
        fct=-Delta*np.sign(Omega/2/np.pi*t-np.floor(Omega/2/np.pi*t)-ratio)
    return fct

def fonctionModKDer(t,x,modulation,frontiere,i): 
    contacti="Contact_"+str(i)
    Delta=float(modulation["DeltaF"])
    if frontiere[contacti]=="masse-ressort module":
        k=0
        if modulation["Synchro"]=="Yes":
            Omega=2*np.pi*float(modulation["Freq"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phi0"])
            else:
                phi0=0    
        else:
            Omega=2*np.pi*float(modulation["FreqF"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phiF0"])
            else:
                phi0=0
    elif frontiere[contacti]=="masse-ressort module dephase":
        if modulation["Synchro"]=="Yes":
            Omega=2*np.pi*float(modulation["Freq"])
            if modulation["Dephasage"]=="xct":
                celF=float(modulation["ModC"])
                if celF==0:
                    Omega=0
                    k=0
                else:
                    k=Omega/celF
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phi0"])
            else:
                phi0=0    
        else:
            Omega=2*np.pi*float(modulation["FreqF"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phiF0"])
            else:
                phi0=0  
            if modulation["Dephasage"]=="xct":
                celF=float(modulation["ModF"])
                if celF==0:
                    Omega=0
                    k=0
                else:
                    k=Omega/celF        
    if modulation["TypeMod"]=="Sinus":
        fct=Omega*Delta*np.cos(Omega*t-k*x+phi0)
    if modulation["TypeMod"]=="Cosinus":
        fct=-Omega*Delta*np.sin(Omega*t-k*x+phi0)
    if modulation["TypeMod"]=="Cosinus QP":
        fct=-Omega*Delta*(0.4*np.sin(Omega*t-k*x)+0.13/np.pi*np.sin(Omega/np.pi*t-k*x)+0.25*np.sqrt(2)*np.sin(Omega*np.sqrt(2)*t-k*x)+0.22*np.sqrt(3)*np.sin(Omega*np.sqrt(3)*t-k*x))
    if modulation["TypeMod"]=="QP-f":
        fct=Delta/2*(np.cos(Omega*t-k*x+phi0)*Omega+np.sin(4*Omega/np.pi*t-k*x)*4*Omega/np.pi)
    if modulation["TypeMod"]=="QP-r":
        fct=Delta/2*(Omega*np.cos(Omega*t-k*x+phi0)+Omega*np.sqrt(2)*np.cos(Omega*np.sqrt(2)*t-k*x))
    elif modulation["TypeMod"]=="Square":
        fct=0*Delta*np.sign(np.sin(Omega*t-k*x+phi0))
    elif modulation["TypeMod"]=="SquareSym":
        fct=0*Delta*np.sign(np.cos(Omega*t-k*x+phi0))
    elif modulation["TypeMod"]=="Square NS":
        ratio=float(modulation["Ratio"])
        fct=-0*Delta*np.sign(Omega/2/np.pi*t-np.floor(Omega/2/np.pi*t)-ratio)   
    return fct

def fonctionModMDer(t,x,modulation,frontiere,i): 
    contacti="Contact_"+str(i)
    Delta=float(modulation["DeltaM"])
    if frontiere[contacti]=="masse-ressort module":
        k=0
        if modulation["Synchro"]=="Yes":
            Omega=2*np.pi*float(modulation["Freq"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phi0"])
            else:
                phi0=0    
        else:
            Omega=2*np.pi*float(modulation["FreqM"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phiF0"])
            else:
                phi0=0
    elif frontiere[contacti]=="masse-ressort module dephase":
        if modulation["Synchro"]=="Yes":
            Omega=2*np.pi*float(modulation["Freq"])
            if modulation["Dephasage"]=="xct":
                celM=float(modulation["ModC"])
                if celM==0:
                    Omega=0
                    k=0
                else:
                    k=Omega/celM
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phi0"])
            else:
                phi0=0    
        else:
            Omega=-2*np.pi*float(modulation["FreqM"])
            if modulation["Dephasage0"]=="Yes":
                phi0=float(modulation["phiM0"])
            else:
                phi0=0  
            if modulation["Dephasage"]=="xct":
                celM=float(modulation["ModM"])
                if celM==0:
                    Omega=0
                    k=0
                else:
                    k=Omega/celM        
    if modulation["TypeMod"]=="Sinus":
        fct=Omega*Delta*np.cos(Omega*t-k*x+phi0)
    if modulation["TypeMod"]=="Cosinus":
        fct=-Omega*Delta*np.sin(Omega*t-k*x+phi0)
    if modulation["TypeMod"]=="Cosinus QP":
        fct=-Omega*Delta*(0.4*np.sin(Omega*t-k*x)+0.13/np.pi*np.sin(Omega/np.pi*t-k*x)+0.25*np.sqrt(2)*np.sin(Omega*np.sqrt(2)*t-k*x)+0.22*np.sqrt(3)*np.sin(Omega*np.sqrt(3)*t-k*x))
    if modulation["TypeMod"]=="QP-f":
        fct=Delta/2*(np.cos(Omega*t-k*x+phi0)*Omega+np.sin(4*Omega/np.pi*t-k*x)*4*Omega/np.pi)
    if modulation["TypeMod"]=="QP-r":
        fct=Delta/2*(Omega*np.cos(Omega*t-k*x+phi0)+Omega*np.sqrt(2)*np.cos(Omega*np.sqrt(2)*t-k*x))
    elif modulation["TypeMod"]=="Square":
        fct=0*Delta*np.sign(np.sin(Omega*t-k*x+phi0))
    elif modulation["TypeMod"]=="SquareSym":
        fct=0*Delta*np.sign(np.cos(Omega*t-k*x+phi0))
    elif modulation["TypeMod"]=="Square NS":
        ratio=float(modulation["Ratio"])
        fct=-0*Delta*np.sign(Omega/2/np.pi*t-np.floor(Omega/2/np.pi*t)-ratio)   
    return fct
