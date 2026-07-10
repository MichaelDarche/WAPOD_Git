#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 19:18:16 2025

@author: michael
"""

import sympy as syp
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############# Symbols ###################
t,tau=syp.symbols("t tau", real=True)
omega=syp.symbols("omega",nonzero=True)
C0,M0,Omega=syp.symbols("C0 M0 Omega", positive=True)
C1,M1=syp.symbols("C1 M1", positive=True)
delta,Delta=syp.symbols("delta Delta",real=True)
rho0,c0=syp.symbols("rho0 c0", positive=True)
rho1,c1=syp.symbols("rho1 c1", positive=True)
rho2,c2=syp.symbols("rho2 c2", positive=True)
rho,c=syp.symbols("rho c", positive=True)
D0,D1=syp.symbols("D0 D1")
##########################################
E=rho*c**2
Z=rho*c

A=syp.Matrix([[0,-1/rho],[-E,0]])
##########################################
phi=syp.Function('phi')(t)
phiM=0*syp.Function('phi_M')(t)
phiC=0*syp.Function('phi_C')(t)
#phi=(1+delta*syp.sin(t))
Cp=C0*phi
#Cp=phi
Cp1=Cp.diff(t)
Cp2=Cp1.diff(t)
#Mp=M0*phi
Mp=Z**2*Cp
Mp1=Mp.diff(t)
Mp2=Mp1.diff(t)


psiC=syp.Function('psi_C')(t)
psi=syp.Function('psi')(t)
Cd=C1*psi
# Cd=C0*(1-delta*syp.sin(t-tau))
Cd1=Cd.diff(t)
Cd2=Cd1.diff(t)
# psiM=syp.Function('psi_M')(t)
Md=Z**2*Cd#*psi
# Md=M0*(1-delta*syp.sin(t-tau))
Md1=Md.diff(t)
Md2=Md1.diff(t)
# C1=syp.Function('C1')(t)
# Cp1=C1.diff(t)
# M1=syp.Function('M1')(t)
# Mp1=M1.diff(t)
#########################################
# N1P = syp.Matrix([[1, -Cd1/2], 
#                   [-Md1/2, 1]])

# N1M = syp.Matrix([[1, Cd1/2], 
#                   [Md1/2, 1]])

# NP = syp.Matrix([[1, -Cp1/2], 
#                   [-Mp1/2, 1]])

# NM = syp.Matrix([[1, Cp1/2], 
#                   [Mp1/2, 1]])

# # Calcul des inverses
# N1Pm1 = N1P.inv()
# NMm1 = NM.inv()
####################################
# leftN=N1Pm1*N1M
# rightN=NMm1*NP


N1P = syp.Matrix([[1, -Cd1/2, -Cd*E/2, 0], 
                  [-Md1/2, 1, 0, -Md/2/rho],
                  [0, -Cd2/2, -Cd1*E, 1/rho],
                  [-Md2/2, 0, E, -Md1/rho]])

N1M = syp.Matrix([[1, Cd1/2, Cd*E/2, 0], 
                  [Md1/2, 1, 0,Md/2/rho],
                  [0, Cd2/2, Cd1*E, 1/rho],
                  [Md2/2, 0, E, Md1/rho]])

NP = syp.Matrix([[1, -Cp1/2, -Cp*E/2, 0], 
                  [-Mp1/2, 1, 0, -Mp/2/rho],
                  [0, -Cp2/2, -Cp1*E, 1/rho],
                  [-Mp2/2, 0, E, -Mp1/rho]])

NM = syp.Matrix([[1, Cp1/2, Cp*E/2, 0], 
                  [Mp1/2, 1, 0,Mp/2/rho],
                  [0, Cp2/2, Cp1*E, 1/rho],
                  [Mp2/2, 0, E, Mp1/rho]])

ZER = syp.Matrix([[1, 0, 0, 0], 
                  [0, 1, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]])
ZER = syp.Matrix([[0, 0, 0, 0], 
                  [0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [0, 0, 0, 0]])

# Calcul des inverses
N1Pm1 = N1P.inv()
NMm1 = NM.inv()

N1Mm1=N1M.inv()
####################################
leftN=syp.simplify(N1Pm1*N1M)
left2=syp.simplify(N1Mm1*N1P)
rightN=syp.simplify(NMm1*NP)

eq=syp.Eq(leftN,rightN)
Transfer=left2*rightN
eq2=syp.Eq(Transfer.simplify(),ZER)

