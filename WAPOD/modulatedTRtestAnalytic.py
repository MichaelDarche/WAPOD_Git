#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 09:09:26 2024

@author: darche
"""
import sympy as syp
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############# Symbols ###################
t,x=syp.symbols("t x", real=True)
omega=syp.symbols("omega",nonzero=True)
K0,M0,Omega=syp.symbols("K0 M0 Omega", positive=True)
delta,Delta=syp.symbols("delta Delta",real=True)
rho0,c0=syp.symbols("rho0 c0", positive=True)
rho1,c1=syp.symbols("rho1 c1", positive=True)
R0,R1,r0,r1=syp.symbols("R0 R1 r0 r1")
T0,T1,t0,t1=syp.symbols("T0 T1 t0 t1")
#********* Polynomes *****************
h,H,Z=syp.symbols("h H Z",nonzero=True)
############ Constants ####################
Z0=rho0*c0
Z1=rho1*c1
omegap=omega+Omega
omegam=omega-Omega
# omegap=syp.Abs(omega+Omega)
# omegam=syp.Abs(omega-Omega)
############# Variables ####################
M=M0*(1+Delta/(2*syp.I)*(syp.exp(syp.I*Omega*t)-syp.exp(-syp.I*Omega*t)))
F=1/K0*(1+delta/(2*syp.I)*(syp.exp(syp.I*Omega*t)-syp.exp(-syp.I*Omega*t)))

# M=M0*(1+Delta*(syp.sin(Omega*t)))
# F=1/K0*(1+delta*(syp.sin(Omega*t)))
########## Derivatives #####################
dtM=syp.diff(M,t)
dtF=syp.diff(F,t)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############# Fields ######################
uI=syp.exp(syp.I*(omega*(t-x/c0)))
sI=-syp.I*omega*Z0*syp.exp(syp.I*(omega*(t-x/c0)))

# Nfft=4096*2
# df=10*omega/Nt
# Nf=10*Nt
# dtf=Tmax/Nf
# TFt=[]
# for j in range(Nf):
#     f=j*df/(2*syp.pi)
#     TF=0
#     for k in range(Nf):
#         if k==0 or k==int(Nf):
#             TF=TF+1/2*syp.exp(syp.I*(omega*(dtf*k)))*syp.exp(-2j*syp.pi*f*df*k*dtf)*dtf
#         else:
#             TF=TF+syp.exp(syp.I*(omega*(dtf*k)))*syp.exp(-2j*syp.pi*f*df*k*dtf)*dtf
#     TFt.append(TF)


#%%%%%%%% Order n %%%%%%%%%%%%%%%%%%%
def somme_variables_symboliques(n):
    # Declarer une variable symbolique i
    l = syp.symbols('l')

    # Creer une somme symbolique avec 2n+1 termes
    uRpn = syp.Sum(syp.Indexed('Rp', l)*syp.exp(syp.I*(omega+l*Omega)*(t+x/c0)), (l, 1, n))#+syp.Indexed('rp', l)*syp.exp(-syp.I*(omega+l*Omega)*(t+x/c0))
    uRmn = syp.Sum(syp.Indexed('rm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t+x/c0)), (l, 1, n))#+syp.Indexed('rm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t+x/c0))
    uRn=R0+R1*syp.exp(syp.I*(omega*(t+x/c0)))+uRpn+uRmn
    
    sRpn = syp.Sum(Z0*syp.I*(omega+l*Omega)*(syp.Indexed('Rp', l)*syp.exp(syp.I*(omega+l*Omega)*(t+x/c0))), (l, 1, n))
    sRmn = -syp.Sum(Z0*syp.I*(omega-l*Omega)*(syp.Indexed('rm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t+x/c0))), (l, 1, n))
    sRn=Z0*syp.I*omega*R1*syp.exp(syp.I*(omega*(t+x/c0)))+sRpn+sRmn
    
    uTpn = syp.Sum(syp.Indexed('Tp', l)*syp.exp(syp.I*(omega+l*Omega)*(t-x/c1)), (l, 1, n))
    uTmn = syp.Sum(syp.Indexed('tm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t-x/c1)), (l, 1, n))
    uTn=T0+T1*syp.exp(syp.I*(omega*(t-x/c1)))+uTpn+uTmn
    
    sTpn = -syp.Sum(Z1*syp.I*(omega+l*Omega)*(syp.Indexed('Tp', l)*syp.exp(syp.I*(omega+l*Omega)*(t-x/c1))), (l, 1, n))
    sTmn = +syp.Sum(Z1*syp.I*(omega-l*Omega)*(syp.Indexed('tm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t-x/c1))), (l, 1, n))
    sTn=-Z1*syp.I*omega*T1*syp.exp(syp.I*(omega*(t-x/c1)))+sTpn+sTmn
    # Simplifier l'expression
    matCoeffs=syp.zeros(2*(2*n+1),1)
    matCoeffs[2*n]=R1
    matCoeffs[2*n+1]=T1
    for l in range(1,n+1):
        matCoeffs[2*(n-l)]=(syp.Indexed('rm', l))
        matCoeffs[2*(n+l)]=(syp.Indexed('Rp', l))
        matCoeffs[2*(n-l)+1]=(syp.Indexed('tm', l))
        matCoeffs[2*(n+l)+1]=(syp.Indexed('Tp', l))
    return uRn.simplify(),sRn.simplify(),uTn.simplify(),sTn.simplify(),matCoeffs
def somme_CS(n):
    # Declarer une variable symbolique i
    l = syp.symbols('l')

    # Creer une somme symbolique avec 2n+1 termes
    uRpn = syp.Sum(syp.Indexed('Rp', l)*syp.cos((omega+l*Omega)*(t+x/c0))+syp.Indexed('rp', l)*syp.sin((omega+l*Omega)*(t+x/c0)), (l, 1, n))
    uRmn = syp.Sum(syp.Indexed('Rm', l)*syp.cos((omega-l*Omega)*(t+x/c0))+syp.Indexed('rm', l)*syp.sin((omega-l*Omega)*(t+x/c0)), (l, 1, n))
    uRn=R0+R1*syp.cos((omega*(t+x/c0)))+r1*syp.sin((omega*(t+x/c0)))+uRpn+uRmn
    
    sRpn = syp.Sum(Z0*(omega+l*Omega)*(syp.Indexed('Rp', l)*syp.sin((omega+l*Omega)*(t+x/c0))-syp.Indexed('rp', l)*syp.cos((omega+l*Omega)*(t+x/c0))), (l, 1, n))
    sRmn = syp.Sum(Z0*(omega-l*Omega)*(syp.Indexed('Rm', l)*syp.exp(syp.I*(omega-l*Omega)*(t+x/c0))-syp.Indexed('rm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t+x/c0))), (l, 1, n))
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
#% n-Displacement
n=1
uRn,sRn,uTn,sTn,matCoeffs=somme_variables_symboliques(n)
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

# T1N=RTNum.args[0][3]
# T1pN=RTNum.args[0][4]
# T1mN=RTNum.args[0][1]
# Norm=T1N/T1N
# Normp=T1pN/T1N
# Normm=T1mN/T1N

######################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### n-order modulation #############
csVn=(vTn-vI-vRn).expand().subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(syp.I*Omega*t),H)
moVn=1/2*(vTn+vI+vRn).expand().subs(x,0)
csSn=(sTn-sI-sRn).expand().subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(syp.I*Omega*t),H)
moSn=1/2*(sTn+sI+sRn).expand().subs(x,0)
ccVn=(dtF*moSn+F*syp.diff(moSn,t)).expand().subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(syp.I*Omega*t),H)
ccSn=(dtM*moVn+M*syp.diff(moVn,t)).expand().subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(syp.I*Omega*t),H)
## to simplify
## interface en x=0

PolyHV=syp.poly(((csVn-ccVn)*h*H).expand(),H)
# if PolyHV.degree(H)!=2*(n+1):
#     PolyHV=syp.poly(((csVn-ccVn)*h*H**(n+1)).expand(),H)
PolyHS=syp.poly(((csSn-ccSn)*h*H).expand(),H)
# if PolyHS.degree(H)!=2*(n+1):
#     PolyHS=syp.poly(((csSn-ccSn)*h*H**(n+1)).expand(),H)
    
polVCoeffs=PolyHV.all_coeffs()
polSCoeffs=PolyHS.all_coeffs()


eqsA=list()
for i in range(n+1):
    print(i+1)
    hV2A=syp.poly(polVCoeffs[i+1],h).all_coeffs()
    hS2A=syp.poly(polSCoeffs[i+1],h).all_coeffs()
    if i==n+1:
        eqsA.append(syp.Eq((hV2A[0].expand()).simplify(),0))
        eqsA.append(syp.Eq((hS2A[0]).simplify(),0))
    else:
        eqsA.append(syp.Eq((hV2A[0].expand()).simplify(),0))
        eqsA.append(syp.Eq((hS2A[0]).simplify(),0))
        eqsA.append(syp.Eq((hV2A[2].expand()).simplify(),0))
        eqsA.append(syp.Eq((hS2A[2]).simplify(),0))
    # eqsA.append(syp.Eq((hV2A[2].expand()/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(c1,c0).subs(rho1,rho0).simplify(),0))
    # eqsA.append(syp.Eq((hS2A[2]/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(c1,c0).subs(rho1,rho0).simplify(),0))
coeff_matrixNA, constantsNA = syp.linear_eq_to_matrix(eqsA,list(matCoeffs))

param={omega:(30*2*syp.pi).evalf(),Omega:(100*2*syp.pi).evalf(),rho0:1200,rho1:1200,c0:2800,c1:2800,Delta:0.9,delta:0.9,M0:10000,K0:2.45e9}

eqsANum=list()
for i in range(n+1):
    hV2A=syp.poly(polVCoeffs[i+1],h).all_coeffs()
    hS2A=syp.poly(polSCoeffs[i+1],h).all_coeffs()
    if i==n+1:
        eqsANum.append(syp.Eq((hV2A[0].expand()).subs(param).simplify(),0))
        eqsANum.append(syp.Eq((hS2A[0]).subs(param).simplify(),0))
    else:
        eqsANum.append(syp.Eq((hV2A[0].expand()).subs(param).simplify(),0))
        eqsANum.append(syp.Eq((hS2A[0]).subs(param).simplify(),0))
        eqsANum.append(syp.Eq((hV2A[2].expand()).subs(param).simplify(),0))
        eqsANum.append(syp.Eq((hS2A[2]).subs(param).simplify(),0))
    # eqsA.append(syp.Eq((hV2A[2].expand()/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(c1,c0).subs(rho1,rho0).simplify(),0))
    # eqsA.append(syp.Eq((hS2A[2]/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(c1,c0).subs(rho1,rho0).simplify(),0))
coeff_matrixNANum, constantsNANum = syp.linear_eq_to_matrix(eqsANum,list(matCoeffs))
RTN2=syp.linsolve((coeff_matrixNANum, constantsNANum),list(matCoeffs))

coefsR=[]
coefsRNorm=[1]
coefsR.append((RTN2.args[0][2*n]))
for i in range(1,n+1):
    coefsR.append((RTN2.args[0][2*(n-i)]))#/coefsR[0])
    coefsRNorm.append((RTN2.args[0][2*(n-i)])/coefsR[0])
    coefsR.append((RTN2.args[0][2*(n+i)]))#/coefsR[0])
    coefsRNorm.append((RTN2.args[0][2*(n+i)])/coefsR[0])
    
    
coefsT=[]
coefsTNorm=[1]
coefsT.append(RTN2.args[0][2*n+1])
for i in range(1,n+1):
    coefsT.append(RTN2.args[0][2*(n-i)+1].evalf())#/coefsR[0])
    coefsTNorm.append((RTN2.args[0][2*(n-i)+1])/coefsT[0])
    coefsT.append((RTN2.args[0][2*(n+i)+1]))#/coefsR[0])
    coefsTNorm.append((RTN2.args[0][1*(n+i)+1])/coefsT[0])

x0=syp.symbols("x0",nonzero=True)
VTnO=(vTn).subs({T1:syp.exp(syp.I*(omega)/c0*(x-x0))*coefsT[0],t1:syp.exp(-syp.I*(omega)/c0*(x-x0))*coefsT[1]}).subs(param)
#VTnm=vTn.subs({omega:-200*2*syp.pi.evalf(),Omega:-17*2*syp.pi.evalf(),c0:2800,c1:2800,T1:0})
