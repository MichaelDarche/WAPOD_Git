#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:50:43 2023

@author: michael
"""
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("../")

import scheme1D
import source
import inputReading



################# Input files 
#% Demarrer file
configuration=inputReading.readfile('Demarrer.txt')
mat=inputReading.readfile('Milieu.txt')
frontiere=inputReading.readfile('Frontiere.txt')     
sourcef=inputReading.readfile('Source.txt')   
################# Parameters of the simulation   
X,xmin,xmax,Nx,dx=inputReading.geometric(configuration)
fs,CFL,tshift=inputReading.source(sourcef)
Nmat,rho,cm,E,A=inputReading.material(mat)
dt=scheme1D.timeStep(CFL,dx,cm)
lambdaOnde=cm/fs
xsmin=cm*tshift-lambdaOnde
xsmax=cm*tshift

#%

#% Initial conditions
V0=scheme1D.cauchyProblem(configuration,frontiere,mat,sourcef)

###### Mise en oeuvre
Nt=300
Tmax=Nt*dt
K=2
Res=scheme1D.FD(V0,Nt,Nx,K,A,dt,dx)

plt.figure()
plt.plot(X,V0[0,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$v$ (m/s)")
plt.title("LW : $t=$0 s")
plt.savefig('Figures/LW-V_0.eps', format='eps')

plt.figure()
plt.plot(X,V0[1,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$\sigma$ (Pa)")
plt.title("LW : $t=$0 s")
plt.savefig('Figures/LW-S_0.eps', format='eps')

plt.figure()
plt.plot(X,Res[0,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$v$ (m/s)")
plt.title("LW : $t=T_{max}$")
plt.savefig('Figures/LW-V_300.eps', format='eps')

plt.figure()
plt.plot(X,Res[1,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$\sigma$ (Pa)")
plt.title("LW : $t=T_{max}$")
plt.savefig('Figures/LW-S_300.eps', format='eps')


###############################

with open("DataBruno/LW-S-Iter0.txt", "r") as fichier:
    # Lire le contenu du fichier ligne par ligne
    lignes = fichier.readlines()
    
colonne1 = []
colonne2 = []

# Parcourir les lignes et extraire les donnees des deux colonnes
for ligne in lignes:
    # Divisez la ligne en mots en utilisant l'espace comme separateur
    mots = ligne.split()
    
    # Verifiez s'il y a au moins deux mots dans la ligne
    if len(mots) >= 2:
        try:
            # Convertissez les deux mots en nombres (entiers ou flottants)
            valeur_colonne1 = float(mots[0])
            valeur_colonne2 = float(mots[1])
            
            # Ajoutez les valeurs aux listes correspondantes
            colonne1.append(valeur_colonne1)
            colonne2.append(valeur_colonne2)
        except ValueError:
            # Ignorez les lignes qui ne contiennent pas deux nombres
            pass
        
LWS0x = np.array(colonne1)
LWS0y = np.array(colonne2)        


with open("DataBruno/LW-V-Iter0.txt", "r") as fichier:
    # Lire le contenu du fichier ligne par ligne
    lignes = fichier.readlines()

colonne1 = []
colonne2 = []

# Parcourir les lignes et extraire les donnees des deux colonnes
for ligne in lignes:
    # Divisez la ligne en mots en utilisant l'espace comme separateur
    mots = ligne.split()
    
    # Verifiez s'il y a au moins deux mots dans la ligne
    if len(mots) >= 2:
        try:
            # Convertissez les deux mots en nombres (entiers ou flottants)
            valeur_colonne1 = float(mots[0])
            valeur_colonne2 = float(mots[1])
            
            # Ajoutez les valeurs aux listes correspondantes
            colonne1.append(valeur_colonne1)
            colonne2.append(valeur_colonne2)
        except ValueError:
            # Ignorez les lignes qui ne contiennent pas deux nombres
            pass
LWV0x = np.array(colonne1)
LWV0y = np.array(colonne2)

with open("DataBruno/LW-V-Iter300.txt", "r") as fichier:
    # Lire le contenu du fichier ligne par ligne
    lignes = fichier.readlines()

colonne1 = []
colonne2 = []

# Parcourir les lignes et extraire les donnees des deux colonnes
for ligne in lignes:
    # Divisez la ligne en mots en utilisant l'espace comme separateur
    mots = ligne.split()
    
    # Verifiez s'il y a au moins deux mots dans la ligne
    if len(mots) >= 2:
        try:
            # Convertissez les deux mots en nombres (entiers ou flottants)
            valeur_colonne1 = float(mots[0])
            valeur_colonne2 = float(mots[1])
            
            # Ajoutez les valeurs aux listes correspondantes
            colonne1.append(valeur_colonne1)
            colonne2.append(valeur_colonne2)
        except ValueError:
            # Ignorez les lignes qui ne contiennent pas deux nombres
            pass
LWV300x = np.array(colonne1)
LWV300y = np.array(colonne2)

plt.figure()
plt.plot(X,V0[0,:])
plt.plot(X,LWV0y.T,'--')
plt.xlabel("$x$ (m)")
plt.ylabel("$v$ (m/s)")
plt.title("LW : $t=0 s$")
plt.savefig('Figures/Comp_LW-V_0.eps', format='eps')

plt.figure()
plt.plot(X,V0[0,:]-LWV0y.T)
plt.xlabel("$x$ (m)")
plt.ylabel("$v-v_{Bruno}$ (m/s)")
plt.title("LW : $t=0 s$")
plt.savefig('Figures/Diff_LW-V_0.eps', format='eps')

plt.figure()
plt.plot(X,V0[1,:])
plt.plot(X,LWS0y.T,'--')
plt.xlabel("$x$ (m)")
plt.ylabel("$\sigma$ (Pa)")
plt.title("LW : $t=0 s$")
plt.savefig('Figures/Comp_LW-S_0.eps', format='eps')

plt.figure()
plt.plot(X,V0[1,:]-LWS0y.T)
plt.xlabel("$x$ (m)")
plt.ylabel("$\sigma-\sigma_{Bruno}$ (Pa)")
plt.title("LW : $t=0 s$")
plt.savefig('Figures/Diff_LW-S_0.eps', format='eps')

plt.figure()
plt.plot(X,Res[0,:])
plt.plot(X,LWV300y.T,'--')
plt.xlabel("$x$ (m)")
plt.ylabel("$v$ (m/s)")
plt.title("LW : $t=T_{max}$")
plt.savefig('Figures/Comp_LW-V_300.eps', format='eps')

plt.figure()
plt.plot(X,Res[0,:]-LWV300y.T)
plt.xlabel("$x$ (m)")
plt.ylabel("$v-v_{Bruno}$ (m/s)")
plt.title("LW : $t=T_{max}$")
plt.savefig('Figures/Diff_LW-V_300.eps', format='eps')


K=4
Res=scheme1D.FD(V0,Nt,Nx,K,A,dt,dx)

plt.figure()
plt.plot(X,V0[0,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$v$ (m/s)")
plt.title("ADER4 : $t=$0 s")
plt.savefig('Figures/ADER4-V_0.eps', format='eps')

plt.figure()
plt.plot(X,V0[1,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$\sigma$ (Pa)")
plt.title("ADER4 : $t=$0 s")
plt.savefig('Figures/ADER4-S_0.eps', format='eps')


plt.figure()
plt.plot(X,Res[0,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$v$ (m/s)")
plt.title("ADER4 : $t=T_{max}$")
plt.savefig('Figures/ADER4-V_300.eps', format='eps')

plt.figure()
plt.plot(X,Res[1,:])
plt.xlabel("$x$ (m)")
plt.ylabel("$\sigma$ (Pa)")
plt.title("ADER4 : $t=T_{max}$")
plt.savefig('Figures/ADER4-S_300.eps', format='eps')

with open("DataBruno/ADER4-V-Iter300.txt", "r") as fichier:
    # Lire le contenu du fichier ligne par ligne
    lignes = fichier.readlines()
    
colonne1 = []
colonne2 = []

# Parcourir les lignes et extraire les donnees des deux colonnes
for ligne in lignes:
    # Divisez la ligne en mots en utilisant l'espace comme separateur
    mots = ligne.split()
    
    # Verifiez s'il y a au moins deux mots dans la ligne
    if len(mots) >= 2:
        try:
            # Convertissez les deux mots en nombres (entiers ou flottants)
            valeur_colonne1 = float(mots[0])
            valeur_colonne2 = float(mots[1])
            
            # Ajoutez les valeurs aux listes correspondantes
            colonne1.append(valeur_colonne1)
            colonne2.append(valeur_colonne2)
        except ValueError:
            # Ignorez les lignes qui ne contiennent pas deux nombres
            pass
A4V300x = np.array(colonne1)
A4V300y = np.array(colonne2)

plt.figure()
plt.plot(X,Res[0,:])
plt.plot(X,A4V300y.T,'--')
plt.xlabel("$x$ (m)")
plt.ylabel("$v$ (m/s)")
plt.title("ADER4 : $t=T_{max}$")
plt.savefig('Figures/Comp_ADER4-V_300.eps', format='eps')

plt.figure()
plt.plot(X,Res[0,:]-A4V300y.T)
plt.xlabel("$x$ (m)")
plt.ylabel("$v-v_{Bruno}$ (m/s)")
plt.title("ADER4 : $t=T_{max}$")
plt.savefig('Figures/Diff_ADER4-V_300.eps', format='eps')
