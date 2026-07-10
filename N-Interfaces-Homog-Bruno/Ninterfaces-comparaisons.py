#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 13:52:51 2025

@author: michael
"""

import os
#os.chdir("/home/darche/Documents/ModulatedMetamaterials/WAPRO_TD/1-Interface_Imerged_PS_modulated2/")
os.chdir("/home/michael/Documents/ModulatedMetamaterials/WAPRO_TD/N-Interfaces-Homog-Bruno/")

import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

import sys
sys.path.append("../WAPOD/")

import inputReading
import scheme1D_TDold as scheme1D
import figures as ffg

configuration=inputReading.readfile('Demarrer.txt')
mat=inputReading.readfile('Milieu.txt')
frontiere=inputReading.readfile('Frontiere.txt')     
sourcef=inputReading.readfile('Source.txt')   
modulation=inputReading.readfile('Modulation.txt')
#frontierestot=inputReading.createfront(frontiere,modulation)
################ Parameters of the simulation   
X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(configuration)
fs,CFL,x_0=inputReading.source(sourcef)
Nmat,rho,cm=inputReading.material(mat)
alpha=inputReading.frontiere(frontiere)
rhox,Celx=inputReading.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
dt=min(scheme1D.timeStep(CFL,dx,cm))
fm=modulation["Freq"]


def nomsfichiersU(fc,fm,T,Nx):
    nomUmicro="Configs/Direct/Direct-U-Fc"+str(fc)+"-Fm"+str(fm)+"-T"+str(T)+"-Nx"+str(Nx)+".txt"
    nomUhomog0="Configs/SMmod/Ordre0/SMmod-U0-Fc"+str(fc)+"-Fm"+str(fm)+"-T"+str(T)+"-Nx"+str(Nx)+".txt"
    nomUhomog2="Configs/SMmod/Ordre2/SMmod-U2corr-Fc"+str(fc)+"-Fm"+str(fm)+"-T"+str(T)+"-Nx"+str(Nx)+".txt"
    return nomUmicro,nomUhomog0,nomUhomog2

def nomsfichiersV(fc,fm,T,Nx):
    nomVmicro="Configs/Direct/Direct-V-Fc"+str(fc)+"-Fm"+str(fm)+"-T"+str(T)+"-Nx"+str(Nx)+".txt"
    nomVhomog0="Configs/SMmod/Ordre0/SMmod-V0-Fc"+str(fc)+"-Fm"+str(fm)+"-T"+str(T)+"-Nx"+str(Nx)+".txt"
    nomVhomog2="Configs/SMmod/Ordre2/SMmod-V2corr-Fc"+str(fc)+"-Fm"+str(fm)+"-T"+str(T)+"-Nx"+str(Nx)+".txt"
    return nomVmicro,nomVhomog0,nomVhomog2


def resultsB(fc,fm,T,Nx):
    nomUmic,nomUhom0,nomUhom2=nomsfichiersU(fc,fm,T,Nx)
    nomVmic,nomVhom0,nomVhom2=nomsfichiersV(fc,fm,T,Nx)
    Xmic,Umic=ffg.txt2np(nomUmic)
    Xmic,Vmic=ffg.txt2np(nomVmic)
    Xhom0,Uhom0=ffg.txt2np(nomUhom0)
    Xhom0,Vhom0=ffg.txt2np(nomVhom0)
    Xhom2,Uhom2=ffg.txt2np(nomUhom2)
    Xhom2,Vhom2=ffg.txt2np(nomVhom2)
    return Xmic,Xhom0,Xhom2,Umic,Vmic,Uhom0,Vhom0,Uhom2,Vhom2


##### fm = 16 Hz #####
Xmic,Xhom0,Xhom2,Umic,Vmic,Uhom0,Vhom0,Uhom2,Vhom2=resultsB(16, 16, "018", Nx-1)

plt.figure()
for i in range(alpha.size):
    plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(Xmic,Umic,marker=r'$\bigcirc$',label=r'$U_h$',color='black',s=60)
plt.plot(Xmic,Umic,linestyle=None,marker='o',color='white',markersize=4)
plt.plot(Xhom0,Uhom0,linewidth=3,label=r'$U_0$')
plt.plot(Xhom2,Uhom2,linewidth=2,label=r'$U_2$')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.title("$t=$ 0.18 s",fontsize=16)
plt.legend(fontsize=14)
plt.tight_layout()
namefig="U_K_16_16.eps"
plt.savefig(str(namefig), format='eps')

plt.figure()
for i in range(alpha.size):
    plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(Xmic,Vmic,marker=r'$\bigcirc$',label=r'$V_h$',color='black',s=60)
plt.plot(Xmic,Vmic,linestyle=None,marker='o',color='white',markersize=4)
plt.plot(Xhom0,Vhom0,linewidth=3,label=r'$V_0$')
plt.plot(Xhom2,Vhom2,linewidth=2,label=r'$V_2$')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$V$ (m/s)',fontsize=16)
plt.title("$t=$ 0.18 s",fontsize=16)
plt.legend(fontsize=14)
plt.tight_layout()
namefig="V_K_16_16.eps"
plt.savefig(str(namefig), format='eps')

##### fm = 32 Hz #####
Xmic,Xhom0,Xhom2,Umic,Vmic,Uhom0,Vhom0,Uhom2,Vhom2=resultsB(16, 32, "018", Nx-1)

plt.figure()
for i in range(alpha.size):
    plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(Xmic,Umic,marker=r'$\bigcirc$',label=r'$U_h$',color='black',s=60)
plt.plot(Xmic,Umic,linestyle=None,marker='o',color='white',markersize=4)
plt.plot(Xhom0,Uhom0,linewidth=3,label=r'$U_0$')
plt.plot(Xhom2,Uhom2,linewidth=2,label=r'$U_2$')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.title("$t=$ 0.18 s",fontsize=16)
plt.legend(fontsize=14)
plt.tight_layout()
namefig="U_K_16_32.eps"
plt.savefig(str(namefig), format='eps')

plt.figure()
for i in range(alpha.size):
    plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(Xmic,Vmic,marker=r'$\bigcirc$',label=r'$V_h$',color='black',s=60)
plt.plot(Xmic,Vmic,linestyle=None,marker='o',color='white',markersize=4)
plt.plot(Xhom0,Vhom0,linewidth=3,label=r'$V_0$')
plt.plot(Xhom2,Vhom2,linewidth=2,label=r'$V_2$')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$V$ (m/s)',fontsize=16)
plt.title("$t=$ 0.18 s",fontsize=16)
plt.legend(fontsize=14)
plt.tight_layout()
namefig="V_K_16_32.eps"
plt.savefig(str(namefig), format='eps')

##### fm = 48 Hz #####
Xmic,Xhom0,Xhom2,Umic,Vmic,Uhom0,Vhom0,Uhom2,Vhom2=resultsB(16, 48, "018", Nx-1)

plt.figure()
for i in range(alpha.size):
    plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(Xmic,Umic,marker=r'$\bigcirc$',label=r'$U_h$',color='black',s=60)
plt.plot(Xmic,Umic,linestyle=None,marker='o',color='white',markersize=4)
plt.plot(Xhom0,Uhom0,linewidth=3,label=r'$U_0$')
plt.plot(Xhom2,Uhom2,linewidth=2,label=r'$U_2$')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.title("$t=$ 0.18 s",fontsize=16)
plt.legend(fontsize=14)
plt.tight_layout()
namefig="U_K_16_48.eps"
plt.savefig(str(namefig), format='eps')

plt.figure()
for i in range(alpha.size):
    plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(Xmic,Vmic,marker=r'$\bigcirc$',label=r'$V_h$',color='black',s=60)
plt.plot(Xmic,Vmic,linestyle=None,marker='o',color='white',markersize=4)
plt.plot(Xhom0,Vhom0,linewidth=3,label=r'$V_0$')
plt.plot(Xhom2,Vhom2,linewidth=2,label=r'$V_2$')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$V$ (m/s)',fontsize=16)
plt.title("$t=$ 0.18 s",fontsize=16)
plt.legend(fontsize=14)
plt.tight_layout()
namefig="V_K_16_48.eps"
plt.savefig(str(namefig), format='eps')
##### fm = 64 Hz #####
Xmic,Xhom0,Xhom2,Umic,Vmic,Uhom0,Vhom0,Uhom2,Vhom2=resultsB(16, 64, "018", Nx-1)

plt.figure()
for i in range(alpha.size):
    plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(Xmic,Umic,marker=r'$\bigcirc$',label=r'$U_h$',color='black',s=60)
plt.plot(Xmic,Umic,linestyle=None,marker='o',color='white',markersize=4)
plt.plot(Xhom0,Uhom0,linewidth=3,label=r'$U_0$')
plt.plot(Xhom2,Uhom2,linewidth=2,label=r'$U_2$')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.title("$t=$ 0.18 s",fontsize=16)
plt.legend(fontsize=14)
plt.tight_layout()
namefig="U_K_16_64.eps"
plt.savefig(str(namefig), format='eps')

plt.figure()
for i in range(alpha.size):
    plt.axvline(x=alpha[i],color='gray',linewidth=0.75)
plt.scatter(Xmic,Vmic,marker=r'$\bigcirc$',label=r'$V_h$',color='black',s=60)
plt.plot(Xmic,Vmic,linestyle=None,marker='o',color='white',markersize=4)
plt.plot(Xhom0,Vhom0,linewidth=3,label=r'$V_0$')
plt.plot(Xhom2,Vhom2,linewidth=2,label=r'$V_2$')
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$V$ (m/s)',fontsize=16)
plt.title("$t=$ 0.18 s",fontsize=16)
plt.legend(fontsize=14)
plt.tight_layout()
namefig="V_K_16_64.eps"
plt.savefig(str(namefig), format='eps')