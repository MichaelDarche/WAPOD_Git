#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 17:07:27 2024

@author: darche
"""

import numpy as np
import sympy as syp
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############# Symbols ##################
t,x=syp.symbols("t x")
omega=syp.symbols("omega",nonzero=True)
K0,M0,Omega=syp.symbols("K0 M0 Omega", positive=True)
delta,Delta=syp.symbols("delta Delta",real=True)
rho0,c0=syp.symbols("rho0 c0", positive=True)
rm1,rm2,r0,rp1,rp2=syp.symbols("rm1 rm2 r0 rp1 rp2")
tm1,tm2,t0,tp1,tp2=syp.symbols("tm1 tm2 t0 tp1 tp2")
p0,phi=syp.symbols("p0 phi")
Zm0,Zv0,m=syp.symbols("Zm0 Zv0 m")

omegam2=omega-2*Omega
omegam1=omega-1*Omega
omegap1=omega+1*Omega
omegap2=omega+2*Omega

k=omega/c0
km2=omegam2/c0
km1=omegam2/c0
kp1=omegam2/c0
kp2=omegap2/c0
#########################################
pi=p0*syp.exp(syp.I*(omega*t-k*x))
prm2=rm2*p0*syp.exp(syp.I*(omegam2*t+km2*x))
prm1=rm1*p0*syp.exp(syp.I*(omegam1*t+km1*x))
pr0=r0*p0*syp.exp(syp.I*(omega*t+k*x))
prp1=rp1*p0*syp.exp(syp.I*(omegap1*t+kp1*x))
prp2=rp2*p0*syp.exp(syp.I*(omegap2*t+kp2*x))
pr=prm2+prm1+pr0+prp1+prp2
ptm2=tm2*p0*syp.exp(syp.I*(omegam2*t-km2*x))
ptm1=tm1*p0*syp.exp(syp.I*(omegam1*t-km1*x))
pt0=t0*p0*syp.exp(syp.I*(omega*t-k*x))
ptp1=tp1*p0*syp.exp(syp.I*(omegap1*t-kp1*x))
ptp2=tp2*p0*syp.exp(syp.I*(omegap2*t-kp2*x))
pt=ptm2+ptm1+pt0+ptp1+ptp2
############ Constants ####################
Z0=rho0*c0
###########################################
vi=p0/Z0*syp.exp(syp.I*(omega*t-k*x))
vrm2=-rm2*p0/Z0*syp.exp(syp.I*(omegam2*t+km2*x))
vrm1=-rm1*p0/Z0*syp.exp(syp.I*(omegam1*t+km1*x))
vr0=-r0*p0/Z0*syp.exp(syp.I*(omega*t+k*x))
vrp1=-rp1*p0/Z0*syp.exp(syp.I*(omegap1*t+kp1*x))
vrp2=-rp2*p0/Z0*syp.exp(syp.I*(omegap2*t+kp2*x))
vr=vrm2+vrm1+vr0+vrp1+vrp2
vtm2=tm2*p0/Z0*syp.exp(syp.I*(omegam2*t-km2*x))
vtm1=tm1*p0/Z0*syp.exp(syp.I*(omegam1*t-km1*x))
vt0=t0*p0/Z0*syp.exp(syp.I*(omega*t-k*x))
vtp1=tp1*p0/Z0*syp.exp(syp.I*(omegap1*t-kp1*x))
vtp2=tp2*p0/Z0*syp.exp(syp.I*(omegap2*t-kp2*x))
vt=vtm2+vtm1+vt0+vtp1+vtp2
############################################
Zm=Zm0+m*Zv0/2*(syp.exp(syp.I*(Omega*t)))
############################################
csV=-vi-vr+vt
moP=1/2*(pi+pr+pt)
csP=-pi-pr+pt
moV=1/2*(vi+vr+vt)
############################################
h,H=syp.symbols("h H",nonzero=True)
CV=csV.expand().subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(-syp.I*Omega*t),H)
CP=csP.expand().subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(-syp.I*Omega*t),H)
TP=(Zm*moV).expand().subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(-syp.I*Omega*t),H)

eqV=syp.factor(syp.Eq(CV-0,0))
eqP=syp.factor(syp.Eq(CP-TP,0))



PolyHV=syp.poly(((csVn-ccVn)*h*H**(n+1)).expand(),H)
if PolyHV.degree(H)!=2*(n+1):
    PolyHV=syp.poly(((csVn-ccVn)*h*H**(n+2)).expand(),H)
PolyHS=syp.poly(((csSn-ccSn)*h*H**(n+1)).expand(),H)
if PolyHS.degree(H)!=2*(n+1):
    PolyHS=syp.poly(((csSn-ccSn)*h*H**(n+2)).expand(),H)
    
polVCoeffs=PolyHV.all_coeffs()
polSCoeffs=PolyHS.all_coeffs()
eqs=list()
for i in range(2*n+1):
    hV2=syp.poly(polVCoeffs[i+1],h).all_coeffs()
    hS2=syp.poly(polSCoeffs[i+1],h).all_coeffs()
    eqs.append(syp.Eq((hV2[0]).subs(h,syp.exp(syp.I*omega*t)).simplify(),0))
    eqs.append(syp.Eq((hS2[0]).subs(h,syp.exp(syp.I*omega*t)).simplify(),0))
    eqs.append(syp.Eq((hV2[2]).subs(h,syp.exp(syp.I*omega*t)).simplify(),0))
    eqs.append(syp.Eq((hS2[2]).subs(h,syp.exp(syp.I*omega*t)).simplify(),0))


coeff_matrix, constants = syp.linear_eq_to_matrix((eqV1, eqS1), (R1, T1))



vI=syp.diff(uI,t)
vR1=syp.diff(uR1,t)
vT1=syp.diff(uT1,t)
vR2=syp.diff(uR2,t)
vT2=syp.diff(uT2,t)
dSI=syp.diff(sI,t)
dSR1=syp.diff(sR1,t)
dST1=syp.diff(sT1,t)
dSR2=syp.diff(sR2,t)
dST2=syp.diff(sT2,t)

##### 1-order modulation #############
csV1=vT1-vI-vR1
moV1=1/2*(vT1+vI+vR1)
csS1=sT1-sI-sR1
moS1=1/2*(sT1+sI+sR1)
ccV1=dtF*moS1+F*syp.diff(moS1,t)
ccS1=dtM*moV1+M*syp.diff(moV1,t)
## to simplify
## interface en x=0

#syp.exp(-syp.I*omega*t)*
eqV1_trunc=syp.Eq(R1-T1+1,-(R1*Z0-T1*Z1-Z0)/(4*K0)*2*syp.I*omega)
eqS1_trunc=syp.Eq(R1*Z0+T1*Z1-Z0,-M0*(R1+T1+1)/4*2*syp.I*omega)
# coeff_matrix1, constants1 = syp.linear_eq_to_matrix((eqV1_trunc, eqS1_trunc), (R1, T1))
# RT1=syp.linsolve((eqV1_trunc, eqS1_trunc), (R1, T1))
######################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### 2-order modulation #############
csV2=vT2-vI-vR2
moV2=1/2*(vT2+vI+vR2)
csS2=sT2-sI-sR2
moS2=1/2*(sT2+sI+sR2)
ccV2=dtF*moS2+F*syp.diff(moS2,t)
ccS2=dtM*moV2+M*syp.diff(moV2,t)
## to simplify
## interface en x=0
eqV2=syp.Eq(csV2,ccV2).subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(syp.I*Omega*t),H)
eqS2=syp.Eq(csS2,ccS2).subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(syp.I*Omega*t),H)
# coeff_matrix2, constants2 = syp.linear_eq_to_matrix((eqV2, eqS2), (R1,T1,R1p,R1m,T1p,T1m))
# CM2=syp.factor(coeff_matrix2.subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(syp.I*Omega*t),H))
#syp.exp(-syp.I*omega*t)*
eqV20_trunc=syp.Eq(T1-R1-1,delta/(4*K0)*(omegam/omega*(R1m*Z0-T1m*Z1)*(Omega+omegam)+omegap/omega*(R1p*Z0-T1p*Z1)*(Omega-omegap))+syp.I*omega/(2*K0)*(R1*Z0-T1*Z1-Z0))#.subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0)
eqS20_trunc=syp.Eq(Z0-R1*Z0-T1*Z1,Delta*M0/4*((R1m+T1m)*omegam/omega*(Omega+omegam)+(R1p+T1p)*omegap/omega*(Omega-omegap))+syp.I*M0/2*omega*(R1+T1+1))#.subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0)
eqV2p_trunc=syp.Eq((T1m-R1m)*omegam/omega,delta/(4*K0)*(R1*Z0-T1*Z1-Z0)*(Omega-omega)+syp.I*omega/(2*K0)*(R1m*Z0-T1m*Z1)*(omegam/omega)**2)#.subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0)
eqS2p_trunc=syp.Eq(-R1m*Z0*omegam/omega-T1m*Z1*omegam/omega,Delta*M0/4*(Omega-omega)*(R1+T1+1)+syp.I*M0*omega/2*(omegam/omega)**2*(R1m+T1m))#.subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0)
eqV2m_trunc=syp.Eq((T1p-R1p)*omegap/omega,delta/(4*K0)*(R1*Z0-T1*Z1-Z0)*(Omega+omega)+syp.I*omega/(2*K0)*(R1p*Z0-T1p*Z1)*(omegap/omega)**2)#.subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0)
eqS2m_trunc=syp.Eq(-R1p*Z0*omegap/omega-T1p*Z1*omegap/omega,Delta*M0/4*(Omega+omega)*(R1+T1+1)+syp.I*M0*omega/2*(omegap/omega)**2*(R1p+T1p))#.subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0)

coeff_matrix2, constants2 = syp.linear_eq_to_matrix((eqV2p_trunc,eqS2p_trunc,eqV20_trunc,eqS20_trunc,eqV2m_trunc,eqS2m_trunc), (R1m,T1m,R1,T1,R1p,T1p))
# detM=syp.det(coeff_matrix2)
#RT2=syp.linsolve((eqV2p_trunc,eqS2p_trunc,eqV20_trunc,eqS20_trunc,eqV2m_trunc,eqS2m_trunc), (R1m,T1m,R1,T1,T1p,R1p))

#######################################
eqV20_truncN=syp.Eq(T1-R1-1,delta/(4*K0)*(omegam/omega*(R1m*Z0-T1m*Z1)*(Omega+omegam)+omegap/omega*(R1p*Z0-T1p*Z1)*(Omega-omegap))+syp.I*omega/(2*K0)*(R1*Z0-T1*Z1-Z0))#.subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0)#.subs(Delta,0)
eqS20_truncN=syp.Eq(Z0-R1*Z0-T1*Z1,Delta*M0/4*((R1m+T1m)*omegam/omega*(Omega+omegam)+(R1p+T1p)*omegap/omega*(Omega-omegap))+syp.I*M0/2*omega*(R1+T1+1))#.subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0)#.subs(Delta,0)
eqV2m_truncN=syp.Eq((T1m-R1m)*omegam/omega,delta/(4*K0)*(R1*Z0-T1*Z1-Z0)*(Omega-omega)+syp.I*omega/(2*K0)*(R1m*Z0-T1m*Z1)*(omegam/omega)**2).subs(c1,c0)#.subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0)#.subs(Delta,0)
eqS2m_truncN=syp.Eq(-R1m*Z0*omegam/omega-T1m*Z1*omegam/omega,Delta*M0/4*(Omega-omega)*(R1+T1+1)+syp.I*M0*omega/2*(omegam/omega)**2*(R1m+T1m)).subs(c1,c0)#.subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0)#.subs(Delta,0)
eqV2p_truncN=syp.Eq((T1p-R1p)*omegap/omega,delta/(4*K0)*(R1*Z0-T1*Z1-Z0)*(Omega+omega)+syp.I*omega/(2*K0)*(R1p*Z0-T1p*Z1)*(omegap/omega)**2).subs(c1,c0)#.subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0)#.subs(Delta,0)
eqS2p_truncN=syp.Eq(-R1p*Z0*omegap/omega-T1p*Z1*omegap/omega,Delta*M0/4*(Omega+omega)*(R1+T1+1)+syp.I*M0*omega/2*(omegap/omega)**2*(R1p+T1p)).subs(c1,c0)#.subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0)#.subs(Delta,0)

param={omega:(200*2*syp.pi).evalf(),Omega:(17*2*syp.pi).evalf(),rho0:1200,rho1:1200,c0:2800,c1:2800,Delta:0.9,delta:0.9,M0:20000,K0:2.45e9}
eqV20Num=eqV20_truncN.subs(param)
eqS20Num=eqS20_truncN.subs(param)
eqV2pNum=eqV2p_truncN.subs(param)
eqS2pNum=eqS2p_truncN.subs(param)
eqV2mNum=eqV2m_truncN.subs(param)
eqS2mNum=eqS2m_truncN.subs(param)

coeff_matrix2N, constants2N= syp.linear_eq_to_matrix((eqV2mNum,eqS2mNum,eqV20Num,eqS20Num,eqV2pNum,eqS2pNum), (R1m,T1m,R1,T1,R1p,T1p))
RT2Num=syp.linsolve((eqV2pNum,eqS2pNum,eqV20Num,eqS20Num,eqV2mNum,eqS2mNum), (R1m,T1m,R1,T1,R1p,T1p))

# T1N=RTNum.args[0][3]
# T1pN=RTNum.args[0][4]
# T1mN=RTNum.args[0][1]
# Norm=T1N/T1N
# Normp=T1pN/T1N
# Normm=T1mN/T1N
