#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 09:21:01 2023

@author: michael
"""
import numpy as np


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

def materialWillis(mat):
    Nmat=int(mat["nombre de milieux "])
    rho=np.zeros(Nmat)
    cm=np.zeros(Nmat)
    S1=np.zeros(Nmat)
    # E=np.zeros(Nmat)
    for i in range(Nmat):
        nameRho="Rho"+'_'+str(i)
        rho[i]=float(mat[nameRho])
        nameCm="Cel"+'_'+str(i)
        cm[i]=float(mat[nameCm])
        nameS1="S1"+'_'+str(i)
        S1[i]=float(mat[nameS1])/2
        # E[i]=cm[i]**2*rho[i]
        #A=np.array([[0,-1/rho],[-E,0]])
    return Nmat,rho,cm,S1#,E,A

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

def affectMaterialsWillis(alphaF,rho,cm,S1,xmin,xmax,Nx,dx):
    minf=xmin
    alphaF[-1]=xmax
    msup=alphaF[0]
    rhox=np.zeros(Nx)
    Celx=np.zeros(Nx)
    S1x=np.zeros(Nx)
    med=0
    for mx in range(Nx):
        if mx*dx+xmin<minf:
            print("erreur geometrique")
        elif mx*dx+xmin<= msup:
            rhox[mx]=rho[med]
            Celx[mx]=cm[med]
            S1x[mx]=S1[med]
        else:
            med=med+1
            rhox[mx]=rho[med]
            Celx[mx]=cm[med]
            S1x[mx]=S1[med]
            minf=msup
            msup=alphaF[med]
    return rhox,Celx,S1x

def comE(rho,c):
    E=c**2*rho
    return E

def comA(rho,c):
    E=comE(rho,c)
    A=np.array([[0,-1/rho],[-E,0]])
    return A     

def comAWillis(rho,c,S1):
    E=comE(rho,c)
    A=np.array([[-S1/rho,-1/rho],[-(E+S1**2/rho),-S1/rho]])
    return A

def c_Willis(rho,c,S1):
    A=comAWillis(rho,c,S1)
    eigenvalues, eigenvectors = np.linalg.eig(A)
    return eigenvalues
    
         

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
    