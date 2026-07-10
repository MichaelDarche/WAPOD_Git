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
        namedict={}
        if filename=="Milieu.txt":
            line_nb_mat = file.readline().strip()
            part_nb_milieu = line_nb_mat.strip().split(':')
            nb_mat=part_nb_milieu[1]
            print(nb_mat)
            for line in file:
                while line.strip()!="milieu 0":
                    parts = line.strip().split(':')
                    if len(parts) == 2:
                        variable_name = parts[0].strip()
                        value = parts[1].strip()
                        namedict[variable_name] = value
                for i in range(int(nb_mat)):
                    name_mat="milieu "+str(i+1)
                    while line.strip()!=name_mat:
                        parts = line.strip().split(':')
                        if len(parts) == 2:
                            variable_name = parts[0].strip()+'_'+str(i)
                            value = parts[1].strip()
                            namedict[variable_name] = value
        elif filename=="Frontiere.txt":
            line_nb_inter = file.readline().strip()
            part_nb_inter = line_nb_mat.strip().split(':')
            nb_inter=part_nb_inter[1]
            for line in file:
                while line.strip()!="frontiere 1":
                    parts = line.strip().split(':')
                    if len(parts) == 2:
                        variable_name = parts[0].strip()
                        value = parts[1].strip()
                        namedict[variable_name] = value
                    for i in range(int(nb_mat)):
                        name_mat="milieu "+str(i+1)
                        while line.strip()!=name_mat:
                            parts = line.strip().split(':')
                            if len(parts) == 2:
                                variable_name = parts[0].strip()+'_'+str(i)
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
    dx=(xmax-xmin)/(Nx-1)
    X=np.linspace(xmin,xmax,Nx)
    return X,xmin,xmax,Nx,dx

def material(mat):
    Nmat=int(mat["nombre de milieux"])
    rho=np.zeros(Nmat)
    cm=np.zeros(Nmat)
    for i in range(Nmat):
        nameRho="Rho"+'_'+str(i)
        rho="rho"+'_'+str(i)
        rho=float(mat["Rho"])
        nameCm="Cel"+'_'+str(i)
        cm="cm"+'_'+str(i)
        cm=float(mat["Cel"])
        E=cm**2*rho
        A=np.array([[0,-1/rho],[-E,0]])
    return Nmat,rho,cm,E,A

# def frontiere(config,mat,front):
#     for i in range(Nx):
#         test

def source(sourcef):
    fs=float(sourcef["Frequence"])
    CFL=float(sourcef["CFL"])
    tshift=float(sourcef["Tshift"])
    return fs,CFL,tshift