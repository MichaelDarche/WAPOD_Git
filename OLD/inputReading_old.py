#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 09:21:01 2023

@author: michael
"""
import numpy as np
# with open('Demarrer.txt', 'r') as file:
#     configuration = {}
#     for line in file:
#         parts = line.strip().split(':')
#         if len(parts) == 2:
#             variable_name = parts[0].strip()
#             value = parts[1].strip()
#             configuration[variable_name] = value
# #% Milieu.txt file
# with open('Milieu.txt', 'r') as file:
#     mat = {}
#     for line in file:
#         parts = line.strip().split(':')
#         if len(parts) == 2:
#             variable_name = parts[0].strip()
#             value = parts[1].strip()
#             mat[variable_name] = value
# #% Frontiere.txt file
# with open('Frontiere.txt', 'r') as file:
#     frontiere = {}
#     for line in file:
#         parts = line.strip().split(':')
#         if len(parts) == 2:
#             variable_name = parts[0].strip()
#             value = parts[1].strip()
#             frontiere[variable_name] = value
# #% Source.txt file
# with open('Source.txt', 'r') as file:
#     sourcef = {}
#     for line in file:
#         parts = line.strip().split(':')
#         if len(parts) == 2:
#             variable_name = parts[0].strip()
#             value = parts[1].strip()
#             sourcef[variable_name] = value
            
def readfile(filename):
    with open(filename,'r') as file:
        if filename =="Milieu.txt"
        namedict={}
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
    dx=(xmax-xmin)/(Nx-1)
    X=np.linspace(xmin,xmax,Nx)
    return X,xmin,xmax,Nx,dx

def material(mat):
    Nmat=int(mat["nombre de milieux"])
    rho=float(mat["Rho"])
    cm=float(mat["Cel"])
    E=cm**2*rho
    A=np.array([[0,-1/rho],[-E,0]])
    return Nmat,rho,cm,E,A

def source(sourcef):
    fs=float(sourcef["Frequence"])
    CFL=float(sourcef["CFL"])
    tshift=float(sourcef["Tshift"])
    return fs,CFL,tshift