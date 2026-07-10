#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 17:28:26 2023

@author: michael
"""
import inputReading
import numpy as np
import sympy as sp
import math as math

# Perfect contact
def matC(rho,c,order):
    A=inputReading.comA(rho,c)
    matrice=np.zeros([2*(order+1),2*(order+1)])
    matrice[0:2,0:2]=np.identity(2)
    for m in range(1,order+1):
        matrice[2*m:2*m+2,2*m:2*m+2]=np.linalg.matrix_power(-A,m)
    return matrice 

# SM contact
# SpringMass Matrix
def SMMatrix(spring,mass):
    # B=np.zeros([2,2])
    B=np.array([[0,1/spring],[mass,0]])
    return B

def CMMatrix(flex,mass):
    # B=np.zeros([2,2])
    B=sp.Matrix([[0,flex],[mass,0]])
    return B
#    

def dampMatrix(betaC,zetaM):
    # B=np.zeros([2,2])
    D=sp.Matrix([[0,betaC],[zetaM,0]])
    return D

def matC_SM(rho,c,order,B,amorK=0,amorM=0):
    D=dampMatrix(amorK,amorM)
    A=inputReading.comA(rho,c)
    F=np.dot(B,-A)
    matrice=np.zeros([2*(order+1),2*(order+1)])
    matrice[0:2,0:2]=np.identity(2)+D/2
    for m in range(1,order+1):
        matrice[2*m:2*m+2,2*m:2*m+2]=np.linalg.matrix_power(-A,m)-np.dot(D,np.linalg.matrix_power(-A,m))/2
        matrice[2*m-2:2*m,2*m:2*m+2]=np.dot(F,np.linalg.matrix_power(-A,m-1))/2
    return matrice
    

def fctModM(modulation):    
    if modulation["TypeMod"]=="Sinus":
        t, x ,Omega,phi= sp.symbols('t x OmegaM phiM')
        fct=sp.sin(Omega*t+phi)
    if modulation["TypeMod"]=="Cosinus":
        t, x ,Omega,phi= sp.symbols('t x OmegaM phiM')
        fct=sp.cos(Omega*t+phi)    
    if modulation["TypeMod"]=="Cosinus QP":
        t, x ,Omega,phi= sp.symbols('t x OmegaM phiM')
        fct=(0.4*sp.cos(Omega*t+phi)+0.13*sp.cos(Omega/sp.pi*t)+0.25*sp.cos(Omega*sp.sqrt(2)*t)+0.22*sp.cos(Omega*sp.sqrt(3)*t))
    if modulation["TypeMod"]=="QP-f":
        t, x ,Omega,phi= sp.symbols('t x OmegaM phiM')
        fct=(sp.sin(Omega*t+phi)+sp.sin(Omega/sp.pi*t))/2
    if modulation["TypeMod"]=="QP-r":
        t, x ,Omega,phi= sp.symbols('t x OmegaM phiM')
        fct=(sp.sin(Omega*t+phi)+sp.sin(Omega*sp.sqrt(2)*t))/2
    if modulation["TypeMod"]=="Square":
        t, x ,Omega,phi= sp.symbols('t x OmegaM phiM')
        fct=sp.sign(sp.sin(Omega*t+phi))
    if modulation["TypeMod"]=="SquareSym":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=sp.sign(sp.cos(Omega*t+phi))
    if modulation["TypeMod"]=="Square NS":
        t, x ,Omega,phi= sp.symbols('t x OmegaM phiM')
        fct=-sp.sign(Omega/2/sp.pi*t-sp.floor(Omega/2/sp.pi*t)-phi)
    return fct,t,x,Omega,phi

def fctModF(modulation):    
    if modulation["TypeMod"]=="Sinus":
        t, x ,Omega,phi= sp.symbols('t x OmegaF phiF')
        fct=sp.sin(Omega*t+phi)
    if modulation["TypeMod"]=="Cosinus":
        t, x ,Omega,phi= sp.symbols('t x OmegaF phiF')
        fct=sp.cos(Omega*t+phi)
    if modulation["TypeMod"]=="Cosinus QP":
        t, x ,Omega,phi= sp.symbols('t x OmegaM phiM')
        fct=(0.4*sp.cos(Omega*t+phi)+0.13*sp.cos(Omega/sp.pi*t)+0.25*sp.cos(Omega*sp.sqrt(2)*t)+0.22*sp.cos(Omega*sp.sqrt(3)*t))
    if modulation["TypeMod"]=="QP-f":
        t, x ,Omega,phi= sp.symbols('t x OmegaF phiF')
        fct=(sp.sin(Omega*t+phi)+sp.sin(Omega/sp.pi*t))/2
    if modulation["TypeMod"]=="QP-r":
        t, x ,Omega,phi= sp.symbols('t x OmegaF phiF')
        fct=(sp.sin(Omega*t+phi)+sp.sin(Omega*sp.sqrt(2)*t))/2
    if modulation["TypeMod"]=="Square":
        t, x ,Omega,phi= sp.symbols('t x OmegaF phiF')
        fct=sp.sign(sp.sin(Omega*t+phi))
    if modulation["TypeMod"]=="SquareSym":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=sp.sign(sp.cos(Omega*t+phi))
    if modulation["TypeMod"]=="Square NS":
        t, x ,Omega,phi= sp.symbols('t x OmegaF phiF')
        fct=-sp.sign(Omega/2/sp.pi*t-sp.floor(Omega/2/sp.pi*t)-phi)
    return fct,t,x,Omega,phi

def fctModZ(modulation):    
    if modulation["TypeMod"]=="Sinus":
        t, x ,Omega,phi= sp.symbols('t x OmegaZ phiZ')
        fct=sp.sin(Omega*t+phi)
    if modulation["TypeMod"]=="Cosinus":
        t, x ,Omega,phi= sp.symbols('t x OmegaZ phiZ')
        fct=sp.cos(Omega*t+phi)    
    if modulation["TypeMod"]=="Cosinus QP":
        t, x ,Omega,phi= sp.symbols('t x OmegaZ phiZ')
        fct=(0.4*sp.cos(Omega*t+phi)+0.13*sp.cos(Omega/sp.pi*t)+0.25*sp.cos(Omega*sp.sqrt(2)*t)+0.22*sp.cos(Omega*sp.sqrt(3)*t))
    if modulation["TypeMod"]=="QP-f":
        t, x ,Omega,phi= sp.symbols('t x OmegaZ phiZ')
        fct=(sp.sin(Omega*t+phi)+sp.sin(Omega/sp.pi*t))/2
    if modulation["TypeMod"]=="QP-r":
        t, x ,Omega,phi= sp.symbols('t x OmegaZ phiZ')
        fct=(sp.sin(Omega*t+phi)+sp.sin(Omega*sp.sqrt(2)*t))/2
    if modulation["TypeMod"]=="Square":
        t, x ,Omega,phi= sp.symbols('t x OmegaZ phiZ')
        fct=sp.sign(sp.sin(Omega*t+phi))
    if modulation["TypeMod"]=="SquareSym":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=sp.sign(sp.cos(Omega*t+phi))
    if modulation["TypeMod"]=="Square NS":
        t, x ,Omega,phi= sp.symbols('t x OmegaZ phiZ')
        fct=-sp.sign(Omega/2/sp.pi*t-sp.floor(Omega/2/sp.pi*t)-phi)
    return fct,t,x,Omega,phi

def fctModB(modulation):    
    if modulation["TypeMod"]=="Sinus":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=sp.sin(Omega*t+phi)
    if modulation["TypeMod"]=="Cosinus":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=sp.cos(Omega*t+phi)
    if modulation["TypeMod"]=="Cosinus QP":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=(0.4*sp.cos(Omega*t+phi)+0.13*sp.cos(Omega/sp.pi*t)+0.25*sp.cos(Omega*sp.sqrt(2)*t)+0.22*sp.cos(Omega*sp.sqrt(3)*t))
    if modulation["TypeMod"]=="QP-f":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=(sp.sin(Omega*t+phi)+sp.sin(Omega/sp.pi*t))/2
    if modulation["TypeMod"]=="QP-r":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=(sp.sin(Omega*t+phi)+sp.sin(Omega*sp.sqrt(2)*t))/2
    if modulation["TypeMod"]=="Square":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=sp.sign(sp.sin(Omega*t+phi))
    if modulation["TypeMod"]=="SquareSym":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=sp.sign(sp.cos(Omega*t+phi))
    if modulation["TypeMod"]=="Square NS":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=-sp.sign(Omega/2/sp.pi*t-sp.floor(Omega/2/sp.pi*t)-phi)
    return fct,t,x,Omega,phi

def fctModMidpt(frontiere,Ni):    
    if frontiere["TypeMod_"+str(Ni)]=="Sinus":
        t, x ,Omega,phi= sp.symbols('t x OmegaM phiM')
        fct=sp.sin(Omega*t+phi)
    if frontiere["TypeMod_"+str(Ni)]=="Cosinus":
        t, x ,Omega,phi= sp.symbols('t x OmegaM phiM')
        fct=sp.cos(Omega*t+phi)    
    if frontiere["TypeMod_"+str(Ni)]=="Cosinus QP":
        t, x ,Omega,phi= sp.symbols('t x OmegaM phiM')
        fct=(0.4*sp.cos(Omega*t+phi)+0.13*sp.cos(Omega/sp.pi*t)+0.25*sp.cos(Omega*sp.sqrt(2)*t)+0.22*sp.cos(Omega*sp.sqrt(3)*t))
    if frontiere["TypeMod_"+str(Ni)]=="QP-f":
        t, x ,Omega,phi= sp.symbols('t x OmegaM phiM')
        fct=(sp.sin(Omega*t+phi)+sp.sin(Omega/sp.pi*t))/2
    if frontiere["TypeMod_"+str(Ni)]=="QP-r":
        t, x ,Omega,phi= sp.symbols('t x OmegaM phiM')
        fct=(sp.sin(Omega*t+phi)+sp.sin(Omega*sp.sqrt(2)*t))/2
    if frontiere["TypeMod_"+str(Ni)]=="Square":
        t, x ,Omega,phi= sp.symbols('t x OmegaM phiM')
        fct=sp.sign(sp.sin(Omega*t+phi))
    if frontiere["TypeMod_"+str(Ni)]=="SquareSym":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=sp.sign(sp.cos(Omega*t+phi))
    if frontiere["TypeMod_"+str(Ni)]=="Square NS":
        t, x ,Omega,phi= sp.symbols('t x OmegaM phiM')
        fct=-sp.sign(Omega/2/sp.pi*t-sp.floor(Omega/2/sp.pi*t)-phi)
    return fct,t,x,Omega,phi

def fctModFidpt(frontiere,Ni):    
    if frontiere["TypeMod_"+str(Ni)]=="Sinus":
        t, x ,Omega,phi= sp.symbols('t x OmegaF phiF')
        fct=sp.sin(Omega*t+phi)
    if frontiere["TypeMod_"+str(Ni)]=="Cosinus":
        t, x ,Omega,phi= sp.symbols('t x OmegaF phiF')
        fct=sp.cos(Omega*t+phi)
    if frontiere["TypeMod_"+str(Ni)]=="Cosinus QP":
        t, x ,Omega,phi= sp.symbols('t x OmegaM phiM')
        fct=(0.4*sp.cos(Omega*t+phi)+0.13*sp.cos(Omega/sp.pi*t)+0.25*sp.cos(Omega*sp.sqrt(2)*t)+0.22*sp.cos(Omega*sp.sqrt(3)*t))
    if frontiere["TypeMod_"+str(Ni)]=="QP-f":
        t, x ,Omega,phi= sp.symbols('t x OmegaF phiF')
        fct=(sp.sin(Omega*t+phi)+sp.sin(Omega/sp.pi*t))/2
    if frontiere["TypeMod_"+str(Ni)]=="QP-r":
        t, x ,Omega,phi= sp.symbols('t x OmegaF phiF')
        fct=(sp.sin(Omega*t+phi)+sp.sin(Omega*sp.sqrt(2)*t))/2
    if frontiere["TypeMod_"+str(Ni)]=="Square":
        t, x ,Omega,phi= sp.symbols('t x OmegaF phiF')
        fct=sp.sign(sp.sin(Omega*t+phi))
    if frontiere["TypeMod_"+str(Ni)]=="SquareSym":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=sp.sign(sp.cos(Omega*t+phi))
    if frontiere["TypeMod_"+str(Ni)]=="Square NS":
        t, x ,Omega,phi= sp.symbols('t x OmegaF phiF')
        fct=-sp.sign(Omega/2/sp.pi*t-sp.floor(Omega/2/sp.pi*t)-phi)
    return fct,t,x,Omega,phi

def fctModZidpt(frontiere,Ni):    
    if frontiere["TypeMod_"+str(Ni)]=="Sinus":
        t, x ,Omega,phi= sp.symbols('t x OmegaZ phiZ')
        fct=sp.sin(Omega*t+phi)
    if frontiere["TypeMod_"+str(Ni)]=="Cosinus":
        t, x ,Omega,phi= sp.symbols('t x OmegaZ phiZ')
        fct=sp.cos(Omega*t+phi)    
    if frontiere["TypeMod_"+str(Ni)]=="Cosinus QP":
        t, x ,Omega,phi= sp.symbols('t x OmegaZ phiZ')
        fct=(0.4*sp.cos(Omega*t+phi)+0.13*sp.cos(Omega/sp.pi*t)+0.25*sp.cos(Omega*sp.sqrt(2)*t)+0.22*sp.cos(Omega*sp.sqrt(3)*t))
    if frontiere["TypeMod_"+str(Ni)]=="QP-f":
        t, x ,Omega,phi= sp.symbols('t x OmegaZ phiZ')
        fct=(sp.sin(Omega*t+phi)+sp.sin(Omega/sp.pi*t))/2
    if frontiere["TypeMod_"+str(Ni)]=="QP-r":
        t, x ,Omega,phi= sp.symbols('t x OmegaZ phiZ')
        fct=(sp.sin(Omega*t+phi)+sp.sin(Omega*sp.sqrt(2)*t))/2
    if frontiere["TypeMod_"+str(Ni)]=="Square":
        t, x ,Omega,phi= sp.symbols('t x OmegaZ phiZ')
        fct=sp.sign(sp.sin(Omega*t+phi))
    if frontiere["TypeMod_"+str(Ni)]=="SquareSym":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=sp.sign(sp.cos(Omega*t+phi))
    if frontiere["TypeMod_"+str(Ni)]=="Square NS":
        t, x ,Omega,phi= sp.symbols('t x OmegaZ phiZ')
        fct=-sp.sign(Omega/2/sp.pi*t-sp.floor(Omega/2/sp.pi*t)-phi)
    return fct,t,x,Omega,phi

def fctModBidpt(frontiere,Ni):    
    if frontiere["TypeMod_"+str(Ni)]=="Sinus":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=sp.sin(Omega*t+phi)
    if frontiere["TypeMod_"+str(Ni)]=="Cosinus":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=sp.cos(Omega*t+phi)
    if frontiere["TypeMod_"+str(Ni)]=="Cosinus QP":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=(0.4*sp.cos(Omega*t+phi)+0.13*sp.cos(Omega/sp.pi*t)+0.25*sp.cos(Omega*sp.sqrt(2)*t)+0.22*sp.cos(Omega*sp.sqrt(3)*t))
    if frontiere["TypeMod_"+str(Ni)]=="QP-f":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=(sp.sin(Omega*t+phi)+sp.sin(Omega/sp.pi*t))/2
    if frontiere["TypeMod_"+str(Ni)]=="QP-r":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=(sp.sin(Omega*t+phi)+sp.sin(Omega*sp.sqrt(2)*t))/2
    if frontiere["TypeMod_"+str(Ni)]=="Square":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=sp.sign(sp.sin(Omega*t+phi))
    if frontiere["TypeMod_"+str(Ni)]=="SquareSym":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=sp.sign(sp.cos(Omega*t+phi))
    if frontiere["TypeMod_"+str(Ni)]=="Square NS":
        t, x ,Omega,phi= sp.symbols('t x OmegaB phiB')
        fct=-sp.sign(Omega/2/sp.pi*t-sp.floor(Omega/2/sp.pi*t)-phi)
    return fct,t,x,Omega,phi


def dtF(formalF,t,n,modulation):
    dnB={}
    Fn=formalF
    for i in range(n+1):
        if i==0:
            dnB[f"d^{0}"] = Fn      
        else:
            if modulation["TypeMod"]=="Square":
                Fn=0
            elif modulation["TypeMod"]=="SquareSym":
                Fn=0
            elif modulation["TypeMod"]=="Square NS":
                Fn=0
            else:
                Fn=Fn.diff(t)
            dnB[f"d^{i}"] = Fn
    return dnB

def dtA(formalE,t,n,modulation):
    dnE={}
    En=formalE
    for i in range(n+1):
        if i==0:
            dnE[f"d^{-1}"] = 0*En
            dnE[f"d^{0}"] = En      
        else:
            if modulation["TypeMod"]=="Square":
                En=0
            elif modulation["TypeMod"]=="SquareSym":
                En=0
            elif modulation["TypeMod"]=="Square NS":
                En=0
            else:
                En=En.diff(t)
            dnE[f"d^{i}"] = En
    return dnE

def dtFid(formalF,t,n,frontiere,Ni):
    dnB={}
    Fn=formalF
    for i in range(n+1):
        if i==0:
            dnB[f"d^{0}"] = Fn      
        else:
            if frontiere["TypeMod_"+str(Ni)]=="Square":
                Fn=0
            elif frontiere["TypeMod_"+str(Ni)]=="SquareSym":
                Fn=0
            elif frontiere["TypeMod_"+str(Ni)]=="Square NS":
                Fn=0
            else:
                Fn=Fn.diff(t)
            dnB[f"d^{i}"] = Fn
    return dnB

def dtAid(formalE,t,n,frontiere,Ni):
    dnE={}
    En=formalE
    for i in range(n+1):
        if i==0:
            dnE[f"d^{-1}"] = 0*En
            dnE[f"d^{0}"] = En      
        else:
            if frontiere["TypeMod_"+str(Ni)]=="Square":
                En=0
            elif frontiere["TypeMod_"+str(Ni)]=="SquareSym":
                En=0
            elif frontiere["TypeMod_"+str(Ni)]=="Square NS":
                En=0
            else:
                En=En.diff(t)
            dnE[f"d^{i}"] = En
    return dnE


def formatC_SM_modulated_indep(rho,c,order,formalB,x,t,phiM,OmegaM,phiF,OmegaF,DF,DM,frontiere,Ni,formalE):#,alpha,ts,LambdaM,OmegaM,ampF,ampM):
    A=inputReading.comA(rho,c)
    matrice=sp.zeros(2*(order+1))
    matrice=matrice.subs(0,sp.Float(0))
    coefPascaln=np.zeros(order+1)
    coefPascaln[0]=1
    # formalE=immergedInterfaces.CMMatrix(Bmod,Zmod)
    dnB=dtFid(formalB, t, order+1,frontiere,Ni)
    dnA=dtAid(formalE,t,order+1,frontiere,Ni)
    matCP=matC(rho,c,order)
    # D=dampMatrix(amorK,amorM)
    for m in range(0,order):
        coefPascaln1=np.zeros(order+1)
        for i in range(m+2):
            coefPascaln1[i]=coefPascaln[i-1]+coefPascaln[i]
        for c in range(0,m+2):
            n=m+1-c
            dn=f"d^{n}"
            en=f"d^{n-1}"
            if c==m:
                matrice[2*m:2*m+2,2*c:2*c+2]=coefPascaln1[c]/2*np.dot(dnB[dn],np.linalg.matrix_power(-A,c))+sp.Matrix(matCP[2*m:2*m+2,2*c:2*c+2])+coefPascaln[c]/2*np.dot(dnA[en],np.linalg.matrix_power(-A,c))
            else:
                matrice[2*m:2*m+2,2*c:2*c+2]=coefPascaln1[c]/2*np.dot(dnB[dn],np.linalg.matrix_power(-A,c))+sp.Matrix(matCP[2*m:2*m+2,2*c:2*c+2])+coefPascaln[c]/2*np.dot(dnA[en],np.linalg.matrix_power(-A,c))
        coefPascaln=coefPascaln1
    m= order
    coefPascaln1=np.zeros(order+1)
    coefPascaln1[0]=1
    for i in range(1,m+1):
        coefPascaln1[i]=coefPascaln[i-1]+coefPascaln[i]
    for c in range(0,m+1):
        n=m+1-c
        dn=f"d^{n}"
        if c==m:
            matrice[2*m:2*m+2,2*c:2*c+2]=coefPascaln1[c]/2*np.dot(dnB[dn],np.linalg.matrix_power(-A,c))+sp.Matrix(matCP[2*m:2*m+2,2*c:2*c+2])+coefPascaln[c]/2*np.dot(dnA[en],np.linalg.matrix_power(-A,c))
        else:
            matrice[2*m:2*m+2,2*c:2*c+2]=coefPascaln1[c]/2*np.dot(dnB[dn],np.linalg.matrix_power(-A,c))+sp.Matrix(matCP[2*m:2*m+2,2*c:2*c+2])+coefPascaln[c]/2*np.dot(dnA[en],np.linalg.matrix_power(-A,c))
    coefPascaln=coefPascaln1
    # matOrderformal=Nmatrix(alpha,ts,ampF,ampM,OmegaM,LambdaM)
    # CMat=matCP+np.array(matOrderformal)
    return matrice

def formatC_SM_modulated(rho,c,order,formalB,x,t,phiM,OmegaM,phiF,OmegaF,DF,DM,modulation,formalE):#,alpha,ts,LambdaM,OmegaM,ampF,ampM):
    A=inputReading.comA(rho,c)
    matrice=sp.zeros(2*(order+1))
    matrice=matrice.subs(0,sp.Float(0))
    coefPascaln=np.zeros(order+1)
    coefPascaln[0]=1
    # formalE=immergedInterfaces.CMMatrix(Bmod,Zmod)
    dnB=dtF(formalB, t, order+1,modulation)
    dnA=dtA(formalE,t,order+1,modulation)
    matCP=matC(rho,c,order)
    # D=dampMatrix(amorK,amorM)
    for m in range(0,order):
        coefPascaln1=np.zeros(order+1)
        for i in range(m+2):
            coefPascaln1[i]=coefPascaln[i-1]+coefPascaln[i]
        for c in range(0,m+2):
            n=m+1-c
            dn=f"d^{n}"
            en=f"d^{n-1}"
            if c==m:
                matrice[2*m:2*m+2,2*c:2*c+2]=coefPascaln1[c]/2*np.dot(dnB[dn],np.linalg.matrix_power(-A,c))+sp.Matrix(matCP[2*m:2*m+2,2*c:2*c+2])+coefPascaln[c]/2*np.dot(dnA[en],np.linalg.matrix_power(-A,c))
            else:
                matrice[2*m:2*m+2,2*c:2*c+2]=coefPascaln1[c]/2*np.dot(dnB[dn],np.linalg.matrix_power(-A,c))+sp.Matrix(matCP[2*m:2*m+2,2*c:2*c+2])+coefPascaln[c]/2*np.dot(dnA[en],np.linalg.matrix_power(-A,c))
        coefPascaln=coefPascaln1
    m= order
    coefPascaln1=np.zeros(order+1)
    coefPascaln1[0]=1
    for i in range(1,m+1):
        coefPascaln1[i]=coefPascaln[i-1]+coefPascaln[i]
    for c in range(0,m+1):
        n=m+1-c
        dn=f"d^{n}"
        if c==m:
            matrice[2*m:2*m+2,2*c:2*c+2]=coefPascaln1[c]/2*np.dot(dnB[dn],np.linalg.matrix_power(-A,c))+sp.Matrix(matCP[2*m:2*m+2,2*c:2*c+2])+coefPascaln[c]/2*np.dot(dnA[en],np.linalg.matrix_power(-A,c))
        else:
            matrice[2*m:2*m+2,2*c:2*c+2]=coefPascaln1[c]/2*np.dot(dnB[dn],np.linalg.matrix_power(-A,c))+sp.Matrix(matCP[2*m:2*m+2,2*c:2*c+2])+coefPascaln[c]/2*np.dot(dnA[en],np.linalg.matrix_power(-A,c))
    coefPascaln=coefPascaln1
    # matOrderformal=Nmatrix(alpha,ts,ampF,ampM,OmegaM,LambdaM)
    # CMat=matCP+np.array(matOrderformal)
    return matrice


def valueModMatrix(Fmatrix,t,ts):
    matOrder=Fmatrix(ts)
    CMat=np.array(matOrder)
    return CMat
    
def normalizeMatrix(Mat0,Mat1):
    maxi=0.
    Mat0N=np.zeros(Mat0.shape)
    Mat1N=np.zeros(Mat1.shape)
    maxi0=np.abs(Mat0.max())
    maxi1=np.abs(Mat1.max())
    mini0=np.abs(Mat0.min())
    mini1=np.abs(Mat1.min())
    maxi=max(maxi0,maxi1,mini0,mini1)
    # for i in range(Mat0.shape[0]):
        # for j in range(Mat0.shape[0]):
        #     if np.abs(Mat0[i,j])>maxi:
        #         maxi=np.abs(Mat0[i,j])
        #     if np.abs(Mat1[i,j])>maxi:
        #         maxi=np.abs(Mat1[i,j])
    if maxi !=0.:
        Mat0N=Mat0/maxi
        Mat1N=Mat1/maxi
    return Mat0N, Mat1N

def invertMatrix(A):
    sizeOfMatrix=np.shape(A)
    if sizeOfMatrix[0]==sizeOfMatrix[1]:
        return np.linalg.inv(A)
    else:
        return np.linalg.pinv(A)
    
def matTaylor(dim,order,alpha,x):
    dist=x-alpha
    matrixT=np.zeros([dim,dim*(order+1)])
    for i in range(order+1):
        for j in range(dim):
            coef=1/math.factorial(i)*dist**i
            matrixT[j,dim*i+j]=coef
    return matrixT

def matCInter(mat0,mat1):
    mat0N, Mat1N=normalizeMatrix(mat0,mat1)
    mat1Inv=invertMatrix(Mat1N)
    return np.dot(mat1Inv,mat0N)

def matMR(dim, order,alpha,x,Npt,matCL,matCR):
    M=np.zeros([2*2*Npt,2*(order+1)])
    for xi in range(Npt):
        matT=matTaylor(dim, order, alpha, x[xi])
        M[2*xi:2*xi+2,:]=matT
    for xi in range(Npt,2*Npt): 
        matT=matTaylor(dim, order, alpha, x[xi])
        MCC=matCInter(matCL,matCR)
        M[2*xi:2*xi+2,:]=np.dot(matT,MCC)
    # print(M)
    return M


def matMRSce(dim, order,alpha,x,Npt,matCL,matCR):
    M=np.zeros([2*2*Npt,2*(order+1)])
    S=np.zeros([2*2*Npt,2*(order+1)])
    for xi in range(Npt):
        matT=matTaylor(dim, order, alpha, x[xi])
        M[2*xi:2*xi+2,:]=matT
        S[2*xi:2*xi+2,:]=np.zeros([dim,dim*(order+1)])
    for xi in range(Npt,2*Npt): 
        matT=matTaylor(dim, order, alpha, x[xi])
        MCC=matCInter(matCL,matCR)
        M[2*xi:2*xi+2,:]=np.dot(matT,MCC)
        S[2*xi:2*xi+2,:]=np.dot(matT,invertMatrix(matCL))
    # print(M)
    return M,S

def matMLSce(dim, order,alpha,x,Npt,matCL,matCR):
    M=np.zeros([2*2*Npt,2*(order+1)])
    S=np.zeros([2*2*Npt,2*(order+1)])
    for xi in range(Npt):
        matT=matTaylor(dim, order, alpha, x[xi])
        MCC=matCInter(matCL,matCR)
        M[2*xi:2*xi+2,:]=np.dot(matT,MCC)
        S[2*xi:2*xi+2,:]=np.dot(matT,invertMatrix(matCR))
    for xi in range(Npt,2*Npt): 
        matT=matTaylor(dim, order, alpha, x[xi])
        M[2*xi:2*xi+2,:]=matT
        S[2*xi:2*xi+2,:]=np.zeros([dim,dim*(order+1)])
    # print(M)
    return M,S

# def matML(dim, order,alpha,x,Npt,matCL,matCR):
#     M=np.zeros([2*2*Npt,2*(order+1)])
#     for xi in range(Npt):
#         matT=matTaylor(dim, order, alpha, x[xi])
#         MCC=matCInter(matCR,matCL)
#         M[2*xi:2*xi+2,:]=np.dot(matT,MCC)
#     for xi in range(Npt,2*Npt): 
#         matT=matTaylor(dim, order, alpha, x[xi])
#         M[2*xi:2*xi+2,:]=matT
#     return M

def matMm1(M):
    return invertMatrix(M)

def matR(dim,order,alpha,x,Npt,matCL,matCR):
    M=matMR(dim, order,alpha,x,Npt,matCL,matCR)
    Mm1=matMm1(M)
    # print(Mm1)
    return Mm1

# def matL(dim,order,alpha,x,Npt,matCL,matCR):
#     M=matML(dim, order,alpha,x,Npt,matCL,matCR)
#     Mm1=matMm1(M)
#     return Mm1

def matF(F,M,rho,c,f):
    t=sp.symbols('t')
    omega=2*np.pi*f
    f_in=sp.I*omega*sp.exp(-sp.I*omega*t)
    v_in=-f_in
    s_in=rho*c*f_in
    matF=[[1/2*sp.diff(F*s_in,t)+v_in],[1/2*sp.diff(M*v_in,t)+s_in]]
    return matF

def matdF(vecF,order):
    t=sp.symbols('t')
    Fk=np.zeros(2*(order+1))
    Fk[0]=vecF
    for j in range(1,2*(order+1)):
        Fk[j]=sp.diff(Fk[j-1],t)
    return Fk
        
        


def modifiedValuesRsource(U,Uint,mR,S0,Fk,j,Npt,dim,order, alpha, x):
    Uints=Uint-np.dot(S0,Fk)
    Ue=np.dot(mR,Uints)
    Ur=U.copy()
    for i in range(Npt,2*Npt):
        mtay=matTaylor(dim, order, alpha, x[i])
        Umod=np.dot(mtay, Ue)
        Ur[0,j-Npt+1+i]=Umod[0,0]
        Ur[1,j-Npt+1+i]=Umod[1,0]
    return Ur

def modifiedValuesLsource(U,Uint,mL,S1,Fk,j,Npt,dim,order, alpha, x):
    Uints=Uint+np.dot(S1,Fk)
    Ue=np.dot(mL,Uints)
    Ur=U.copy()
    for i in range(Npt):
        mtay=matTaylor(dim, order, alpha, x[i])
        Umod=np.dot(mtay,Ue)
        Ur[0,j-Npt+1+i]=Umod[0,0]
        Ur[1,j-Npt+1+i]=Umod[1,0]
    return Ur

def modifiedValuesR(U,Uint,mR,j,Npt,dim,order, alpha, x):
    Ue=np.dot(mR,Uint)
    Ur=U.copy()
    for i in range(Npt,2*Npt):
        mtay=matTaylor(dim, order, alpha, x[i])
        Umod=np.dot(mtay, Ue)
        Ur[0,j-Npt+1+i]=Umod[0,0]
        Ur[1,j-Npt+1+i]=Umod[1,0]
    return Ur

def modifiedValuesL(U,Uint,mL,j,Npt,dim,order, alpha, x):
    Ue=np.dot(mL,Uint)
    Ur=U.copy()
    for i in range(Npt):
        mtay=matTaylor(dim, order, alpha, x[i])
        Umod=np.dot(mtay,Ue)
        Ur[0,j-Npt+1+i]=Umod[0,0]
        Ur[1,j-Npt+1+i]=Umod[1,0]
    return Ur

def modifiedValuesRbis(U,Uint,mR,j,Npt,dim,order, alpha, x):
    Ue=np.dot(mR,Uint)
    Ur=U.copy()
    for i in range(Npt,2*Npt):
        mtay=matTaylor(dim, order, alpha, x[i])
        Umod=np.dot(mtay, Ue)
        Ur[0,j-Npt+i]=Umod[0,0]
        Ur[1,j-Npt+i]=Umod[1,0]
    return Ur

def modifiedValuesLbis(U,Uint,mL,j,Npt,dim,order, alpha, x):
    Ue=np.dot(mL,Uint)
    Ur=U.copy()
    for i in range(Npt):
        mtay=matTaylor(dim, order, alpha, x[i])
        Umod=np.dot(mtay,Ue)
        Ur[0,j-Npt+i]=Umod[0,0]
        Ur[1,j-Npt+i]=Umod[1,0]
    return Ur
# def modifiedValuesR(U,xstep,alpha):
#     mT=matP(U,xstep,dim,order,alpha,x,Npt,rhoL,cL,rhoR,cR,dx)
#     Ur[:,j+Npt-2]=
#     Ur[:,j+Npt-1]=
#     Ur[:,j+Npt]=
#     U[:,j-Npt+1:j+Npt]=Ur
#     return U
    
    

# ESIM=3
# order=2*ESIM-1
# Npt=3
# Xp=np.linspace(-Npt+1, Npt,2*Npt)
# rhoL=rho[0]
# rhoR=rho[1]
# cL=cm[0]
# cR=cm[1]

# matCL=matC(rhoL,cL,order)
# matCR=matC(rhoR,cR,order)
# testcm0=matCInter(matCL,matCR)

# Testlor=matTaylor(2,order,Npt/2,Xp[3])

# test=matMR(2, order,0.5,Xp,Npt,matCL,matCR)
# #Mm1=invertMatrix(test)

# mT=matR(2,order,0.5,Xp,Npt,matCL,matCR)
        