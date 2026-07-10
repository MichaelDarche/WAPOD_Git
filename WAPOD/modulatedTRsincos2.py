#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 15:41:38 2024

@author: michael
"""
import sympy as syp
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############# Symbols ###################
t,x=syp.symbols("t x")
omega=syp.symbols("omega",nonzero=True)
K0,M0,Omega=syp.symbols("K0 M0 Omega", positive=True)
delta,Delta=syp.symbols("delta Delta",real=True)
rho0,c0=syp.symbols("rho0 c0", positive=True)
rho1,c1=syp.symbols("rho1 c1", positive=True)
R0,R1,r0,r1=syp.symbols("R0 R1 r0 r1")
T0,T1,t0,t1=syp.symbols("T0 T1 t0 t1")
#********* Polynomes *****************
h,H=syp.symbols("h H",nonzero=True)
############ Constants ####################
Z0=rho0*c0
Z1=rho1*c1
omegap=omega+Omega
omegam=omega-Omega
# omegap=syp.Abs(omega+Omega)
# omegam=syp.Abs(omega-Omega)
############# Variables ####################
# M=M0*(1+Delta/(2*syp.I)*(syp.exp(syp.I*Omega*t)-syp.exp(-syp.I*Omega*t)))
# F=1/K0*(1+delta/(2*syp.I)*(syp.exp(syp.I*Omega*t)-syp.exp(-syp.I*Omega*t)))

M=M0*(1+Delta*(syp.sin(Omega*t)))
F=1/K0*(1+delta*(syp.sin(Omega*t)))
########## Derivatives #####################
dtM=syp.diff(M,t)
dtF=syp.diff(F,t)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############# Fields ######################
uI=syp.sin(omega*(t-x/c0))
sI=-syp.I*omega*Z0*syp.cos(omega*(t-x/c0))
#%%%%%%%% Order n %%%%%%%%%%%%%%%%%%%
def somme_variables_symboliques(n):
    # Declarer une variable symbolique i
    l = syp.symbols('l')

    # Creer une somme symbolique avec 2n+1 termes
    uRpn = syp.Sum(syp.Indexed('Rp', l)*syp.exp(syp.I*(omega+l*Omega)*(t+x/c0))+syp.Indexed('rp', l)*syp.exp(-syp.I*(omega+l*Omega)*(t+x/c0)), (l, 1, n))
    uRmn = syp.Sum(syp.Indexed('Rm', l)*syp.exp(syp.I*(omega-l*Omega)*(t+x/c0))+syp.Indexed('rm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t+x/c0)), (l, 1, n))
    uRn=R0+R1*syp.exp(syp.I*(omega*(t+x/c0)))+r1*syp.exp(-syp.I*(omega*(t+x/c0)))+uRpn+uRmn
    
    sRpn = syp.Sum(Z0*syp.I*(omega+l*Omega)*(syp.Indexed('Rp', l)*syp.exp(syp.I*(omega+l*Omega)*(t+x/c0))-syp.Indexed('rp', l)*syp.exp(-syp.I*(omega+l*Omega)*(t+x/c0))), (l, 1, n))
    sRmn = syp.Sum(Z0*syp.I*(omega-l*Omega)*(syp.Indexed('Rm', l)*syp.exp(syp.I*(omega-l*Omega)*(t+x/c0))-syp.Indexed('rm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t+x/c0))), (l, 1, n))
    sRn=Z0*syp.I*omega*R1*syp.exp(syp.I*(omega*(t+x/c0)))-Z0*syp.I*omega*r1*syp.exp(-syp.I*(omega*(t+x/c0)))+sRpn+sRmn
    
    uTpn = syp.Sum(syp.Indexed('Tp', l)*syp.exp(syp.I*(omega+l*Omega)*(t-x/c1))+syp.Indexed('tp', l)*syp.exp(-syp.I*(omega+l*Omega)*(t-x/c1)), (l, 1, n))
    uTmn = syp.Sum(syp.Indexed('Tm', l)*syp.exp(syp.I*(omega-l*Omega)*(t-x/c1))+syp.Indexed('tm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t-x/c1)), (l, 1, n))
    uTn=T0+T1*syp.exp(syp.I*(omega*(t-x/c1)))+t1*syp.exp(-syp.I*(omega*(t-x/c1)))+uTpn+uTmn
    
    sTpn = -syp.Sum(Z1*syp.I*(omega+l*Omega)*(syp.Indexed('Tp', l)*syp.exp(syp.I*(omega+l*Omega)*(t-x/c1))-syp.Indexed('tp', l)*syp.exp(-syp.I*(omega+l*Omega)*(t-x/c1))), (l, 1, n))
    sTmn = -syp.Sum(Z1*syp.I*(omega-l*Omega)*(syp.Indexed('Tm', l)*syp.exp(syp.I*(omega-l*Omega)*(t-x/c1))-syp.Indexed('tm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t-x/c1))), (l, 1, n))
    sTn=-Z1*syp.I*omega*T1*syp.exp(syp.I*(omega*(t-x/c1)))+Z1*syp.I*omega*t1*syp.exp(-syp.I*(omega*(t-x/c1)))+sTpn+sTmn
    # Simplifier l'expression
    matCoeffs=syp.zeros(4*(2*n+1),1)
    matCoeffs[4*n]=R1
    matCoeffs[4*n+1]=r1
    matCoeffs[4*n+2]=T1
    matCoeffs[4*n+3]=t1
    for l in range(1,n+1):
        matCoeffs[4*(n-l)]=(syp.Indexed('Rm', l))
        matCoeffs[4*(n+l)]=(syp.Indexed('Rp', l))
        matCoeffs[4*(n-l)+1]=(syp.Indexed('rm', l))
        matCoeffs[4*(n+l)+1]=(syp.Indexed('rp', l))
        matCoeffs[4*(n-l)+2]=(syp.Indexed('Tm', l))
        matCoeffs[4*(n+l)+2]=(syp.Indexed('Tp', l))
        matCoeffs[4*(n-l)+3]=(syp.Indexed('tm', l))
        matCoeffs[4*(n+l)+3]=(syp.Indexed('tp', l))
    return uRn.simplify(),sRn.simplify(),uTn.simplify(),sTn.simplify(),matCoeffs
def somme_CS(n):
    # Declarer une variable symbolique i
    l = syp.symbols('l')
    # Creer une somme symbolique avec 2n+1 termes
    uRpn = syp.Sum(syp.Indexed('Rp', l)*syp.sin((omega+l*Omega)*(t+x/c0))+syp.Indexed('rp', l)*syp.cos((omega+l*Omega)*(t+x/c0)), (l, 1, n))
    uRmn = syp.Sum(syp.Indexed('Rm', l)*syp.sin((omega-l*Omega)*(t+x/c0))+syp.Indexed('rm', l)*syp.cos((omega-l*Omega)*(t+x/c0)), (l, 1, n))
    uRn=R0+R1*syp.sin((omega*(t+x/c0)))+r1*syp.cos((omega*(t+x/c0)))+uRpn+uRmn
    
    sRpn = syp.Sum(Z0*syp.I*(omega+l*Omega)*(syp.Indexed('Rp', l)*syp.cos((omega+l*Omega)*(t+x/c0))-syp.Indexed('rp', l)*syp.sin((omega+l*Omega)*(t+x/c0))), (l, 1, n))
    sRmn = syp.Sum(Z0*syp.I*(omega-l*Omega)*(syp.Indexed('Rm', l)*syp.cos((omega-l*Omega)*(t+x/c0))-syp.Indexed('rm', l)*syp.sin((omega-l*Omega)*(t+x/c0))), (l, 1, n))
    sRn=Z0*syp.I*omega*R1*syp.cos((omega*(t+x/c0)))-Z0*syp.I*omega*r1*syp.sin((omega*(t+x/c0)))+sRpn+sRmn
    
    uTpn = syp.Sum(syp.Indexed('Tp', l)*syp.sin((omega+l*Omega)*(t-x/c1))+syp.Indexed('tp', l)*syp.cos((omega+l*Omega)*(t-x/c1)), (l, 1, n))
    uTmn = syp.Sum(syp.Indexed('Tm', l)*syp.sin((omega-l*Omega)*(t-x/c1))+syp.Indexed('tm', l)*syp.cos((omega-l*Omega)*(t-x/c1)), (l, 1, n))
    uTn=T0+T1*syp.sin((omega*(t-x/c1)))+t1*syp.cos((omega*(t-x/c1)))+uTpn+uTmn
    
    sTpn = -syp.Sum(Z1*syp.I*(omega+l*Omega)*(syp.Indexed('Tp', l)*syp.cos((omega+l*Omega)*(t-x/c1))-syp.Indexed('tp', l)*syp.sin((omega+l*Omega)*(t-x/c1))), (l, 1, n))
    sTmn = -syp.Sum(Z1*syp.I*(omega-l*Omega)*(syp.Indexed('Tm', l)*syp.cos((omega-l*Omega)*(t-x/c1))-syp.Indexed('tm', l)*syp.sin((omega-l*Omega)*(t-x/c1))), (l, 1, n))
    sTn=-Z1*syp.I*omega*T1*syp.cos((omega*(t-x/c1)))+Z1*syp.I*omega*t1*syp.sin((omega*(t-x/c1)))+sTpn+sTmn
    # Simplifier l'expression
    matCoeffs=syp.zeros(4*(2*n+1),1)
    matCoeffs[4*n]=R1
    matCoeffs[4*n+1]=r1
    matCoeffs[4*n+2]=T1
    matCoeffs[4*n+3]=t1
    for l in range(1,n+1):
        matCoeffs[4*(n-l)]=(syp.Indexed('Rm', l))
        matCoeffs[4*(n+l)]=(syp.Indexed('Rp', l))
        matCoeffs[4*(n-l)+1]=(syp.Indexed('rm', l))
        matCoeffs[4*(n+l)+1]=(syp.Indexed('rp', l))
        matCoeffs[4*(n-l)+2]=(syp.Indexed('Tm', l))
        matCoeffs[4*(n+l)+2]=(syp.Indexed('Tp', l))
        matCoeffs[4*(n-l)+3]=(syp.Indexed('tm', l))
        matCoeffs[4*(n+l)+3]=(syp.Indexed('tp', l))
    return uRn.simplify(),sRn.simplify(),uTn.simplify(),sTn.simplify(),matCoeffs
#% n-Displacement
n=1
uRn,sRn,uTn,sTn,matCoeffs=somme_variables_symboliques(n)
uRn,sRn,uTn,sTn,matCoeffs=somme_CS(n)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vI=syp.diff(uI,t)
vRn=syp.diff(uRn,t)
vTn=syp.diff(uTn,t)
dSRn=syp.diff(sRn,t)
dSTn=syp.diff(sTn,t)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%% Analytique %%%%%%%%%%%%%%%
# exAV=(1+syp.I*omega*Z0/(2*K0))*R1-(1+syp.I*omega*Z1/(2*K0))*T1+1-syp.I*omega*Z0/(2*K0)
# eqVA1=syp.Eq(exAV,0).subs(x,0)
# exAS=(Z0+syp.I*omega*M0/(2))*R1-(Z1+syp.I*omega*M0/(2))*T1-Z0+syp.I*omega*M0/(2)
# eqSA1=syp.Eq(exAS,0).subs(x,0)
# coeff_matrixA, constantsA = syp.linear_eq_to_matrix((eqVA1, eqSA1), (R1, T1))
######################################
param={omega:(200*2*syp.pi).evalf(),Omega:(35*2*syp.pi).evalf(),rho0:1200,rho1:1200,c0:2800,c1:2800,Delta:0.9,delta:0.9,M0:20000,K0:2.45e9}

# T1N=RTNum.args[0][3]
# T1pN=RTNum.args[0][4]
# T1mN=RTNum.args[0][1]
# Norm=T1N/T1N
# Normp=T1pN/T1N
# Normm=T1mN/T1N

######################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### n-order modulation #############
csVn=(vTn-vI-vRn).subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(-syp.I*Omega*t),H)
moVn=1/2*(vTn+vI+vRn).subs(x,0)
csSn=(sTn-sI-sRn).subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(-syp.I*Omega*t),H)
moSn=1/2*(sTn+sI+sRn).subs(x,0)
ccVn=(dtF*moSn+F*syp.diff(moSn,t)).expand().subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(-syp.I*Omega*t),H)
ccSn=(dtM*moVn+M*syp.diff(moVn,t)).expand().subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(-syp.I*Omega*t),H)
## to simplify
## interface en x=0

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
    eqs.append(syp.Eq((hV2[0].expand()/syp.I).subs(h,syp.exp(syp.I*omega*t)).simplify(),0))
    eqs.append(syp.Eq((hS2[0]/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).simplify(),0))
    eqs.append(syp.Eq((hV2[2].expand()/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).simplify(),0))
    eqs.append(syp.Eq((hS2[2]/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).simplify(),0))
    
coeff_matrixN, constantsN = syp.linear_eq_to_matrix(eqs,list(matCoeffs))

eqsRL=list()
for i in range(2*n+1):
    hV2=syp.poly(polVCoeffs[i+1],h).all_coeffs()
    hS2=syp.poly(polSCoeffs[i+1],h).all_coeffs()
    eqsRL.append(syp.Eq((hV2[0].expand()/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))
    eqsRL.append(syp.Eq((hS2[0]/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))
    eqsRL.append(syp.Eq((hV2[2].expand()/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))
    eqsRL.append(syp.Eq((hS2[2]/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))
    
    # eqsRL.append(syp.Eq((polVCoeffs[i+1]/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))
    # eqsRL.append(syp.Eq((polSCoeffs[i+1]/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))
    # eqsRL.append(syp.Eq((polVCoeffs[i+1]/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))
    # eqsRL.append(syp.Eq((polSCoeffs[i+1]/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))

coeff_matrixRL, constantsRL = syp.linear_eq_to_matrix(eqsRL,list(matCoeffs))

pRL={c0:2800,rho0:1200,Delta:9/10,M0:20000,omega:30,Omega:100}
eqsRLN=list()
for i in range(2*n+1):    
    hV2=syp.poly(polVCoeffs[i+1],h).all_coeffs()
    hS2=syp.poly(polSCoeffs[i+1],h).all_coeffs()
    eqsRLN.append(syp.Eq((hV2[0].expand()/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).subs(pRL).simplify(),0))
    eqsRLN.append(syp.Eq((hS2[0]/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).subs(pRL).simplify(),0))
    eqsRLN.append(syp.Eq((hV2[2].expand()/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).subs(pRL).simplify(),0))
    eqsRLN.append(syp.Eq((hS2[2]/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).subs(pRL).simplify(),0))

coeff_matrixRLN, constantsRLN = syp.linear_eq_to_matrix(eqsRLN,list(matCoeffs))
RTRLNum=syp.linsolve(eqsRLN,list(matCoeffs))

eqsNum=list()
for i in range(2*n+1):
    hV2=syp.poly(polVCoeffs[i+1],h).all_coeffs()
    hS2=syp.poly(polSCoeffs[i+1],h).all_coeffs()
    eqsNum.append(syp.Eq((hV2[0].expand()/omega/syp.I).subs(param).simplify(),0))
    eqsNum.append(syp.Eq((hS2[0]/omega/syp.I).subs(param).simplify(),0))
    eqsNum.append(syp.Eq((hV2[2].expand()/omega/syp.I).subs(param).simplify(),0))
    eqsNum.append(syp.Eq((hS2[2]/omega/syp.I).subs(param).simplify(),0))
    # eqsNum.append(syp.Eq((polVCoeffs[i+1].expand()/h/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(param).simplify(),0))
    # eqsNum.append(syp.Eq((polSCoeffs[i+1].expand()/h/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(param).simplify(),0))

coeff_matrixNum, constantsNum = syp.linear_eq_to_matrix(eqsNum,list(matCoeffs))
RTNum=syp.linsolve(eqsNum,list(matCoeffs))


coefsR=[]
coefsRNorm=[1]
coefsR.append(syp.Abs(RTNum.args[0][4*n]))
coefsR.append(syp.Abs(RTNum.args[0][4*n+1]))
for i in range(1,n+1):
    coefsR.append(syp.Abs(RTNum.args[0][4*(n-i)]))#/coefsR[0])
    coefsR.append(syp.Abs(RTNum.args[0][4*(n-i)+1]))
    coefsRNorm.append(syp.Abs(RTNum.args[0][4*(n-i)])/coefsR[0])
    coefsRNorm.append(syp.Abs(RTNum.args[0][4*(n-i)+1])/coefsR[0])
    coefsR.append(syp.Abs(RTNum.args[0][4*(n+i)]))#/coefsR[0])
    coefsR.append(syp.Abs(RTNum.args[0][4*(n+i)+1]))#/coefsR[0])
    coefsRNorm.append(syp.Abs(RTNum.args[0][4*(n+i)])/coefsR[0])
    coefsRNorm.append(syp.Abs(RTNum.args[0][4*(n+i)+1])/coefsR[0])
    
    
coefsT=[]
coefsTNorm=[1]
coefsT.append(syp.Abs(RTNum.args[0][4*n+2]))
coefsT.append(syp.Abs(RTNum.args[0][4*n+3]))
for i in range(1,n+1):
    coefsT.append(syp.Abs(RTNum.args[0][4*(n-i)]))#/coefsR[0])
    coefsT.append(syp.Abs(RTNum.args[0][4*(n-i)+1]))
    coefsTNorm.append(syp.Abs(RTNum.args[0][4*(n-i)])/coefsT[0])
    coefsTNorm.append(syp.Abs(RTNum.args[0][4*(n-i)+1])/coefsT[0])
    coefsT.append(syp.Abs(RTNum.args[0][4*(n+i)]))#/coefsR[0])
    coefsT.append(syp.Abs(RTNum.args[0][4*(n+i)+1]))#/coefsR[0])
    coefsTNorm.append(syp.Abs(RTNum.args[0][4*(n+i)])/coefsT[0])
    coefsTNorm.append(syp.Abs(RTNum.args[0][4*(n+i)+1])/coefsT[0])

VTnO=(vTn+syp.I*vTn).subs({omega:23*2*syp.pi.evalf(),Omega:200*2*syp.pi.evalf(),c0:2800,c1:2800,T1:coefsT[0]})
#VTnm=vTn.subs({omega:-200*2*syp.pi.evalf(),Omega:-17*2*syp.pi.evalf(),c0:2800,c1:2800,T1:0})

for j in range(n):
    #VTnm=VTnm.subs(syp.Indexed('Tm', j+1),coefsT[2*j+1]).subs(syp.Indexed('Tp', j+1),coefsT[2*j+2])
    VTnO=VTnO.subs(syp.Indexed('Tm', j+1),coefsT[2*j+1]).subs(syp.Indexed('Tp', j+1),coefsT[2*j+2])

Nt=4096


syp.plot(syp.re(VTnO).subs(x,300), (t, 0, Nt*dt),adaptive=False,nb_of_points=Nt)


# # import numpy as np
# invTF=syp.zeros(Nx,Nt)
# for xi in range(Nx):
#     for ti in range(Nt):
#         Vinv=(VTnO+VTnm)
#         val=Vinv.subs(x,xi*200/Nx).subs(t,ti*dt).evalf()
#         invTF[xi,ti]=val
        
