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

##########################################
rho=rho0
c=c0
E=rho*c**2
Z=rho*c
##########################################
phi=syp.Function('phi')(t)
#phi=(1+delta*syp.sin(t-tau))
C=C0*phi
Cp=C.diff(t)
M=M0*phi
Mp=M.diff(t)

psiC=syp.Function('psi_C')(t)
Cd=psiC
Cd1=Cd.diff(t)
psiM=syp.Function('psi_M')(t)
Md=psiM
Md1=Md.diff(t)
# C1=syp.Function('C1')(t)
# Cp1=C1.diff(t)
# M1=syp.Function('M1')(t)
# Mp1=M1.diff(t)
#########################################
N1P = syp.Matrix([[1, -Cd1/2, 0, 0], 
                  [-Md1/2, 1, 0, 0],
                  [0, 0, -Cd*E/2, 0],
                  [0, 0, 0, -Md/2/rho]])

N1M = syp.Matrix([[1, Cd1/2, 0, 0], 
                  [Md1/2, 1, 0, 0],
                  [0, 0, Cd*E/2, 0],
                  [0, 0, 0, Md/2/rho]])

NP = syp.Matrix([[1, -Cp/2, 0, 0], 
                 [-Mp/2, 1, 0, 0],
                 [0, 0, -C*E/2, 0],
                 [0, 0, 0, -M/2/rho]])

NM = syp.Matrix([[1, Cp/2, 0, 0], 
                 [Mp/2, 1, 0, 0],
                 [0, 0, C*E/2, 0],
                 [0, 0, 0, M/2/rho]])

# Calcul des inverses
N1Pm1 = N1P.inv()
NMm1 = NM.inv()
####################################
leftN=N1Pm1*N1M
rightN=NMm1*NP

eq1=syp.Eq(leftN[0,0],rightN[0,0])
eqG=syp.Eq(leftN,rightN)
eqGred=syp.Eq(leftN[0:2,0:2],rightN[0:2,0:2])

