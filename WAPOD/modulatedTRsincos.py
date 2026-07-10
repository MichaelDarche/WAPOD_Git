#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 15:41:38 2024

@author: michael
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
h,H=syp.symbols("h H",nonzero=True)
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
uI=syp.exp(syp.I*(omega*(t-x/c0)))#+syp.exp(-syp.I*(omega*(t-x/c0)))
sI=-syp.I*omega*Z0*(syp.exp(syp.I*(omega*(t-x/c0))))#+syp.I*omega*Z0*(-syp.exp(syp.I*(omega*(t-x/c0))))

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
    uRmn = syp.Sum(syp.Indexed('Rm', l)*syp.exp(syp.I*(omega-l*Omega)*(t+x/c0)), (l, 1, n))#+syp.Indexed('rm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t+x/c0))
    uRn=R0+R1*syp.exp(syp.I*(omega*(t+x/c0)))+uRpn+uRmn#+r1*syp.exp(-syp.I*(omega*(t+x/c0)))
    
    sRpn = syp.Sum(Z0*syp.I*(omega+l*Omega)*(syp.Indexed('Rp', l)*syp.exp(syp.I*(omega+l*Omega)*(t+x/c0))), (l, 1, n))#-syp.Indexed('rp', l)*syp.exp(-syp.I*(omega+l*Omega)*(t+x/c0))
    sRmn = syp.Sum(Z0*syp.I*(omega-l*Omega)*(syp.Indexed('Rm', l)*syp.exp(syp.I*(omega-l*Omega)*(t+x/c0))), (l, 1, n))#-syp.Indexed('rm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t+x/c0))
    sRn=Z0*syp.I*omega*R1*syp.exp(syp.I*(omega*(t+x/c0)))+sRpn+sRmn#-Z0*syp.I*omega*r1*syp.exp(-syp.I*(omega*(t+x/c0)))
    
    uTpn = syp.Sum(syp.Indexed('Tp', l)*syp.exp(syp.I*(omega+l*Omega)*(t-x/c1)), (l, 1, n))#+syp.Indexed('tp', l)*syp.exp(-syp.I*(omega+l*Omega)*(t-x/c1))
    uTmn = syp.Sum(syp.Indexed('Tm', l)*syp.exp(syp.I*(omega-l*Omega)*(t-x/c1)), (l, 1, n))#+syp.Indexed('tm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t-x/c1))
    uTn=T0+T1*syp.exp(syp.I*(omega*(t-x/c1)))+uTpn+uTmn#+t1*syp.exp(-syp.I*(omega*(t-x/c1)))
    
    sTpn = -syp.Sum(Z1*syp.I*(omega+l*Omega)*(syp.Indexed('Tp', l)*syp.exp(syp.I*(omega+l*Omega)*(t-x/c1))), (l, 1, n))#-syp.Indexed('tp', l)*syp.exp(-syp.I*(omega+l*Omega)*(t-x/c1))
    sTmn = -syp.Sum(Z1*syp.I*(omega-l*Omega)*(syp.Indexed('Tm', l)*syp.exp(syp.I*(omega-l*Omega)*(t-x/c1))), (l, 1, n))#-syp.Indexed('tm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t-x/c1))
    sTn=-Z1*syp.I*omega*T1*syp.exp(syp.I*(omega*(t-x/c1)))+sTpn+sTmn#+Z1*syp.I*omega*t1*syp.exp(-syp.I*(omega*(t-x/c1)))
    # Simplifier l'expression
    matCoeffs=syp.zeros(2*(2*n+1),1)
    matCoeffs[2*n]=R1
    #matCoeffs[4*n+1]=r1
    matCoeffs[2*n+1]=T1
    #matCoeffs[4*n+3]=t1
    for l in range(1,n+1):
        matCoeffs[2*(n-l)]=(syp.Indexed('Rm', l))
        matCoeffs[2*(n+l)]=(syp.Indexed('Rp', l))
        #matCoeffs[4*(n-l)+1]=(syp.Indexed('rm', l))
        #matCoeffs[4*(n+l)+1]=(syp.Indexed('rp', l))
        matCoeffs[2*(n-l)+1]=(syp.Indexed('Tm', l))
        matCoeffs[2*(n+l)+1]=(syp.Indexed('Tp', l))
        #matCoeffs[4*(n-l)+3]=(syp.Indexed('tm', l))
        #matCoeffs[4*(n+l)+3]=(syp.Indexed('tp', l))
    return uRn.simplify(),sRn.simplify(),uTn.simplify(),sTn.simplify(),matCoeffs

def somme_variables_symboliques_tot(n):
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

def somme_variables_symboliquesAbs(n):
    # Declarer une variable symbolique i
    l = syp.symbols('l')

    # Creer une somme symbolique avec 2n+1 termes
    uRpn = syp.Sum(syp.Indexed('Rp', l)*syp.exp(syp.I*(omega+l*Omega)*(t+x/c0))+syp.Indexed('rm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t+x/c0)), (l, 1, n))
    uRmn = 0#syp.Sum(syp.Indexed('Rm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t+x/c0))+syp.Indexed('rm', l)*syp.exp(+syp.I*(omega-l*Omega)*(t+x/c0)), (l, 1, n))
    uRn=R0+R1*syp.exp(syp.I*(omega*(t+x/c0)))+uRpn+uRmn
    
    sRpn = syp.Sum(Z0*syp.I*((omega+l*Omega)*syp.Indexed('Rp', l)*syp.exp(syp.I*(omega+l*Omega)*(t+x/c0))-(omega-l*Omega)*syp.Indexed('rm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t+x/c0))), (l, 1, n))
    sRmn = 0#syp.Sum(-Z0*syp.I*(omega-l*Omega)*(syp.Indexed('Rm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t+x/c0))-syp.Indexed('rm', l)*syp.exp(+syp.I*(omega-l*Omega)*(t+x/c0))), (l, 1, n))
    sRn=Z0*syp.I*omega*R1*syp.exp(syp.I*(omega*(t+x/c0)))-Z0*syp.I*omega*r1*syp.exp(-syp.I*(omega*(t+x/c0)))+sRpn+sRmn
    
    uTpn = syp.Sum(syp.Indexed('Tp', l)*syp.exp(syp.I*(omega+l*Omega)*(t-x/c1))+syp.Indexed('tm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t-x/c1)), (l, 1, n))
    uTmn = 0#syp.Sum(syp.Indexed('Tm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t-x/c1))+syp.Indexed('tm', l)*syp.exp(+syp.I*(omega-l*Omega)*(t-x/c1)), (l, 1, n))
    uTn=T0+T1*syp.exp(syp.I*(omega*(t-x/c1)))+uTpn+uTmn
    
    sTpn = -syp.Sum(Z1*syp.I*((omega+l*Omega)*syp.Indexed('Tp', l)*syp.exp(syp.I*(omega+l*Omega)*(t-x/c1))-(omega-l*Omega)*syp.Indexed('tm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t-x/c1))), (l, 1, n))
    sTmn = -0#syp.Sum(-Z1*syp.I*(omega-l*Omega)*(syp.Indexed('Tm', l)*syp.exp(-syp.I*(omega-l*Omega)*(t-x/c1))-syp.Indexed('tm', l)*syp.exp(+syp.I*(omega-l*Omega)*(t-x/c1))), (l, 1, n))
    sTn=-Z1*syp.I*omega*T1*syp.exp(syp.I*(omega*(t-x/c1)))+sTpn+sTmn
    # Simplifier l'expression
    matCoeffs=syp.zeros(2*(2*n+1),1)
    matCoeffs[2*n]=R1
    #matCoeffs[2*n+1]=r1
    matCoeffs[2*n+1]=T1
    #matCoeffs[2*n+3]=t1
    for l in range(1,n+1):
        #matCoeffs[2*(n-l)]=(syp.Indexed('Rm', l))
        matCoeffs[2*(n+l)]=(syp.Indexed('Rp', l))
        matCoeffs[2*(n-l)]=(syp.Indexed('rm', l))
        #matCoeffs[2*(n+l)+1]=(syp.Indexed('rp', l))
        #matCoeffs[2*(n-l)+2]=(syp.Indexed('Tm', l))
        matCoeffs[2*(n+l)+1]=(syp.Indexed('Tp', l))
        matCoeffs[2*(n-l)+1]=(syp.Indexed('tm', l))
        #matCoeffs[2*(n+l)+3]=(syp.Indexed('tp', l))
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
n=10
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
csVn=(vTn-vI-vRn).expand().subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(-syp.I*Omega*t),H)
moVn=1/2*(vTn+vI+vRn).expand().subs(x,0)
csSn=(sTn-sI-sRn).expand().subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(-syp.I*Omega*t),H)
moSn=1/2*(sTn+sI+sRn).expand().subs(x,0)
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
    #eqs.append(syp.Eq((hV2[2].expand()/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).simplify(),0))
    #eqs.append(syp.Eq((hS2[2]/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).simplify(),0))
    
coeff_matrixN, constantsN = syp.linear_eq_to_matrix(eqs,list(matCoeffs))


eqsA=list()
for i in range(2*n+1):
    hV2A=syp.poly(polVCoeffs[i+1],h).all_coeffs()
    hS2A=syp.poly(polSCoeffs[i+1],h).all_coeffs()
    eqsA.append(syp.Eq((hV2A[0].expand()/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(c1,c0).subs(rho1,rho0).simplify(),0))
    eqsA.append(syp.Eq((hS2A[0]/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(c1,c0).subs(rho1,rho0).simplify(),0))
    #eqsA.append(syp.Eq((hV2A[2].expand()/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(c1,c0).subs(rho1,rho0).simplify(),0))
    #eqsA.append(syp.Eq((hS2A[2]/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(c1,c0).subs(rho1,rho0).simplify(),0))
coeff_matrixNA, constantsNA = syp.linear_eq_to_matrix(eqsA,list(matCoeffs))

#RTana=syp.linsolve((coeff_matrixN, constantsN),list(matCoeffs))

eqsRL=list()
for i in range(2*n+1):
    hV2=syp.poly(polVCoeffs[i+1],h).all_coeffs()
    hS2=syp.poly(polSCoeffs[i+1],h).all_coeffs()
    eqsRL.append(syp.Eq((hV2[0].expand()/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))
    eqsRL.append(syp.Eq((hS2[0]/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))
    #eqsRL.append(syp.Eq((hV2[2].expand()/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))
    #eqsRL.append(syp.Eq((hS2[2]/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))
    
    # eqsRL.append(syp.Eq((polVCoeffs[i+1]/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))
    # eqsRL.append(syp.Eq((polSCoeffs[i+1]/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))
    # eqsRL.append(syp.Eq((polVCoeffs[i+1]/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))
    # eqsRL.append(syp.Eq((polSCoeffs[i+1]/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))

coeff_matrixRL, constantsRL = syp.linear_eq_to_matrix(eqsRL,list(matCoeffs))

pRL={c0:2800,rho0:1200,Delta:9/10,M0:4000,omega:200,Omega:35}
eqsRLN=list()
for i in range(2*n+1):    
    hV2=syp.poly(polVCoeffs[i+1],h).all_coeffs()
    hS2=syp.poly(polSCoeffs[i+1],h).all_coeffs()
    eqsRLN.append(syp.Eq((hV2[0].expand()/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).subs(pRL).simplify(),0))
    eqsRLN.append(syp.Eq((hS2[0]/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).subs(pRL).simplify(),0))
    #eqsRLN.append(syp.Eq((hV2[2].expand()/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).subs(pRL).simplify(),0))
    #eqsRLN.append(syp.Eq((hS2[2]/omega/syp.I).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).subs(pRL).simplify(),0))

coeff_matrixRLN, constantsRLN = syp.linear_eq_to_matrix(eqsRLN,list(matCoeffs))
RTRLNum=syp.linsolve(eqsRLN,list(matCoeffs))

coefsRLR=[]
coefsRLRNorm=[1]
coefsRLR.append((RTRLNum.args[0][2*n]))
#coefsRLR.append(syp.Abs(RTRLNum.args[0][2*n+1]))
for i in range(1,n+1):
    coefsRLR.append((RTRLNum.args[0][2*(n-i)]))#/coefsR[0])
    #coefsRLR.append(syp.Abs(RTRLNum.args[0][2*(n-i)+1]))
    coefsRLRNorm.append((RTRLNum.args[0][2*(n-i)])/coefsRLR[0])
    #coefsRLRNorm.append(syp.Abs(RTRLNum.args[0][2*(n-i)+1])/coefsRLR[0])
    coefsRLR.append((RTRLNum.args[0][2*(n+i)]))#/coefsR[0])
    #coefsRLR.append(syp.Abs(RTRLNum.args[0][2*(n+i)+1]))#/coefsR[0])
    coefsRLRNorm.append((RTRLNum.args[0][2*(n+i)])/coefsRLR[0])
    #coefsRLRNorm.append(syp.Abs(RTRLNum.args[0][2*(n+i)+1])/coefsRLR[0])
    
    
coefsRLT=[]
coefsRLTNorm=[1]
coefsRLT.append((RTRLNum.args[0][2*n+1]))
#coefsRLT.append(syp.Abs(RTRLNum.args[0][2*n+1]))
for i in range(1,n+1):
    coefsRLT.append((RTRLNum.args[0][2*(n-i)+1]))#/coefsR[0])
    #coefsRLT.append(syp.Abs(RTRLNum.args[0][2*(n-i)+3]))
    coefsRLTNorm.append((RTRLNum.args[0][2*(n-i)+1])/coefsRLT[0])
    #coefsRLTNorm.append(syp.Abs(RTRLNum.args[0][2*(n-i)+3])/coefsRLT[0])
    coefsRLT.append((RTRLNum.args[0][2*(n+i)+1]))#/coefsR[0])
    #coefsRLT.append(syp.Abs(RTRLNum.args[0][2*(n+i)+3]))#/coefsR[0])
    coefsRLTNorm.append((RTRLNum.args[0][2*(n+i)+1])/coefsRLT[0])
    #coefsRLTNorm.append(syp.Abs(RTRLNum.args[0][2*(n+i)+3])/coefsRLT[0])

x0=syp.symbols("x0",nonzero=True)
VTnORL=(vTn).subs({omega:200*2*syp.pi.evalf(),Omega:35*2*syp.pi.evalf(),c0:2800,c1:2800,T1:syp.exp(syp.I*(omega)/c0*(x-x0))*coefsRLT[0],t1:syp.exp(-syp.I*(omega)/c0*(x-x0))*coefsRLT[1]})
#VTnm=vTn.subs({omega:-200*2*syp.pi.evalf(),Omega:-17*2*syp.pi.evalf(),c0:2800,c1:2800,T1:0})

for j in range(n):
    #VTnm=VTnm.subs(syp.Indexed('Tm', j+1),coefsT[2*j+1]).subs(syp.Indexed('Tp', j+1),coefsT[2*j+2])
    VTnORL=VTnORL.subs(syp.Indexed('Tm', j+1),syp.exp(syp.I*(omega-j*Omega)/c0*(x-x0))*coefsRLT[2*j+1]).subs(syp.Indexed('Tp', j+1),syp.exp(syp.I*(omega+j*Omega)/c0*(x-x0))*coefsRLT[2*j+2])
    #VTnORL=VTnORL.subs(syp.Indexed('tm', j+1),syp.exp(-syp.I*(omega-j*Omega)/c0*(x-x0))*coefsRLT[4*j+3]).subs(syp.Indexed('tp', j+1),syp.exp(-syp.I*(omega+j*Omega)/c0*(x-x0))*coefsRLT[4*j+5])
Nt=4096
CFL=0.95
dx=0.4
cm=2800
dt=CFL*dx/cm

VTnORL=(VTnORL/c0/2/omega).subs({omega:200*2*syp.pi.evalf(),Omega:35*2*syp.pi.evalf(),c0:2800,x0:-75})
syp.plot(syp.re(VTnORL).subs(x,100), (t, 0, Nt*dt),adaptive=False,nb_of_points=Nt)




param={omega:(500*2*syp.pi).evalf(),Omega:(40*2*syp.pi).evalf(),rho0:1200,rho1:1200,c0:2800,c1:2800,Delta:0.9,delta:0.9,M0:10000,K0:2.45e9}


eqsNum=list()
for i in range(2*n+1):    
    hV2=syp.poly(polVCoeffs[i+1],h).all_coeffs()
    hS2=syp.poly(polSCoeffs[i+1],h).all_coeffs()
    eqsNum.append(syp.Eq((hV2[0].expand()/omega/syp.I).subs(param).simplify(),0))
    eqsNum.append(syp.Eq((hS2[0]/omega/syp.I).subs(param).simplify(),0))
    #eqsNum.append(syp.Eq((hV2[2].expand()/omega/syp.I).subs(param).simplify(),0))
    #eqsNum.append(syp.Eq((hS2[2]/omega/syp.I).subs(c1,c0).subs(param).simplify(),0))

coeff_matrixNum, constantsNum = syp.linear_eq_to_matrix(eqsNum,list(matCoeffs))


def NormalizationParam(Mat,Vec):
    size=Vec.shape[0]
    for i in range(size):
        maxV=max(max(abs(Mat[i,:])),abs(Vec[i]))
        Vec[i]=Vec[i]/maxV
        Mat[i,:]=Mat[i,:]/maxV
    return Mat,Vec

coeff_matrixNum, constantsNum=NormalizationParam(coeff_matrixNum, constantsNum)
RTNum=syp.linsolve((coeff_matrixNum,constantsNum),list(matCoeffs))



coefsR=[]
coefsRNorm=[1,1]
coefsR.append((RTNum.args[0][2*n]))
#coefsR.append(syp.Abs(RTNum.args[0][2*n+1]))
for i in range(1,n+1):
    coefsR.append((RTNum.args[0][2*(n-i)]))#/coefsR[0])
    #coefsR.append(syp.Abs(RTNum.args[0][4*(n-i)+1]))
    coefsRNorm.append((RTNum.args[0][2*(n-i)])/coefsR[0])
    #coefsRNorm.append(syp.Abs(RTNum.args[0][4*(n-i)+1])/coefsR[0])
    coefsR.append((RTNum.args[0][2*(n+i)]))#/coefsR[0])
    #coefsR.append(syp.Abs(RTNum.args[0][4*(n+i)+1]))#/coefsR[0])
    coefsRNorm.append((RTNum.args[0][2*(n+i)])/coefsR[0])
    #coefsRNorm.append(syp.Abs(RTNum.args[0][4*(n+i)+1])/coefsR[0])
    
    
coefsT=[]
coefsTNorm=[1.]
coefsT.append(RTNum.args[0][2*n+1])
#coefsT.append(RTNum.args[0][4*n+3])
for i in range(1,n+1):
    coefsT.append(RTNum.args[0][2*(n-i)+1].evalf())#/coefsR[0])
    #coefsT.append((RTNum.args[0][4*(n-i)+3]))
    coefsTNorm.append((RTNum.args[0][2*(n-i)+1])/coefsT[0])
    #coefsTNorm.append((RTNum.args[0][4*(n-i)+3])/coefsT[0])
    coefsT.append((RTNum.args[0][2*(n+i)+1]))#/coefsR[0])
    #coefsT.append((RTNum.args[0][4*(n+i)+3]))#/coefsR[0])
    coefsTNorm.append((RTNum.args[0][2*(n+i)+1])/coefsT[0])
    #coefsTNorm.append((RTNum.args[0][4*(n+i)+3])/coefsT[0])
coefsTr=[]
coefsTrNorm=[1.]
coefsTr.append(syp.re(RTNum.args[0][2*n+1]))
#coefsT.append(RTNum.args[0][4*n+3])
for i in range(1,n+1):
    coefsTr.append(syp.re(RTNum.args[0][2*(n-i)+1].evalf()))#/coefsR[0])
    #coefsT.append((RTNum.args[0][4*(n-i)+3]))
    coefsTrNorm.append(syp.re((RTNum.args[0][2*(n-i)+1])/coefsTr[0]))
    #coefsTNorm.append((RTNum.args[0][4*(n-i)+3])/coefsT[0])
    coefsTr.append(syp.re((RTNum.args[0][2*(n+i)+1])))#/coefsR[0])
    #coefsT.append((RTNum.args[0][4*(n+i)+3]))#/coefsR[0])
    coefsTrNorm.append(syp.re((RTNum.args[0][2*(n+i)+1])/coefsTr[0]))
    #coefsTNorm.append((RTNum.args[0][4*(n+i)+3])/coefsT[0])
coefsTi=[]
coefsTiNorm=[1.]
coefsTi.append(RTNum.args[0][2*n+1])
#coefsT.append(RTNum.args[0][4*n+3])
for i in range(1,n+1):
    coefsTi.append(syp.im(RTNum.args[0][2*(n-i)+1].evalf()))#/coefsR[0])
    #coefsT.append((RTNum.args[0][4*(n-i)+3]))
    coefsTiNorm.append(syp.im((RTNum.args[0][2*(n-i)+1])/coefsT[0]))
    #coefsTNorm.append((RTNum.args[0][4*(n-i)+3])/coefsT[0])
    coefsTi.append(syp.im((RTNum.args[0][2*(n+i)+1])))#/coefsR[0])
    #coefsT.append((RTNum.args[0][4*(n+i)+3]))#/coefsR[0])
    coefsTiNorm.append(syp.im((RTNum.args[0][2*(n+i)+1])/coefsT[0]))
    #coefsTNorm.append((RTNum.args[0][4*(n+i)+3])/coefsT[0])
coefsTa = []
coefsTaNorm = [1.]
coefsTa.append(syp.Abs(RTNum.args[0][2*n+1]))
#coefsT.append(RTNum.args[0][4*n+3])
for i in range(1, n+1):
    coefsTa.append(syp.Abs(RTNum.args[0][2*(n-i)+1].evalf()))  # /coefsR[0])
    #coefsT.append((RTNum.args[0][4*(n-i)+3]))
    coefsTaNorm.append(syp.Abs((RTNum.args[0][2*(n-i)+1])/coefsTa[0]))
    #coefsTNorm.append((RTNum.args[0][4*(n-i)+3])/coefsT[0])
    coefsTa.append(syp.Abs((RTNum.args[0][2*(n+i)+1])))  # /coefsR[0])
    #coefsT.append((RTNum.args[0][4*(n+i)+3]))#/coefsR[0])
    coefsTaNorm.append(syp.Abs((RTNum.args[0][2*(n+i)+1])/coefsTa[0]))
    #coefsTNorm.append((RTNum.args[0][4*(n+i)+3])/coefsT[0])

x0 = syp.symbols("x0", nonzero=True)
VTnO=(vTn).subs({T1:syp.exp(syp.I*(omega)/c0*(x-x0))*coefsT[0]}).subs(param)
#VTnm=vTn.subs({omega:-200*2*syp.pi.evalf(),Omega:-17*2*syp.pi.evalf(),c0:2800,c1:2800,T1:0})

### Dephasage dependant de Omega
# for j in range(n):
#     #VTnm=VTnm.subs(syp.Indexed('Tm', j+1),coefsT[2*j+1]).subs(syp.Indexed('Tp', j+1),coefsT[2*j+2])
#     VTnO=VTnO.subs(param).subs(syp.Indexed('Tm', j+1),syp.exp(syp.I*(omega-j*Omega)/c0*(x-x0))*coefsT[4*j+2]).subs(syp.Indexed('Tp', j+1),syp.exp(syp.I*(omega+j*Omega)/c0*(x-x0))*coefsT[4*j+4])
#     VTnO=VTnO.subs(param).subs(syp.Indexed('tm', j+1),syp.exp(-syp.I*(omega-j*Omega)/c0*(x-x0))*coefsT[4*j+3]).subs(syp.Indexed('tp', j+1),syp.exp(-syp.I*(omega+j*Omega)/c0*(x-x0))*coefsT[4*j+5])
# Nt=4096*2
### Dephasage independant de Omega
for j in range(n):
    VTnO=VTnO.subs(syp.Indexed('Tm', j+1),coefsT[2*j+1]).subs(syp.Indexed('Tp', j+1),coefsT[2*j+2])
    #VTnO=VTnO.subs(param).subs(syp.Indexed('Tm', j+1),syp.exp(syp.I*(omega)/c0*(x-x0))*coefsT[4*j+2]).subs(syp.Indexed('Tp', j+1),syp.exp(syp.I*(omega)/c0*(x-x0))*coefsT[4*j+4])
    #VTnO=VTnO.subs(param).subs(syp.Indexed('tm', j+1),syp.exp(-syp.I*(omega)/c0*(x-x0))*coefsT[4*j+3]).subs(syp.Indexed('tp', j+1),syp.exp(-syp.I*(omega)/c0*(x-x0))*coefsT[4*j+5])
Nt=4096*4

VTnO=(VTnO/c0/2/omega).subs({x0:100}).subs(param)
syp.plot(syp.re(VTnO).subs(x,100), (t, 0, Nt*dt),adaptive=False,nb_of_points=Nt)




paramS={omega:(30*2*syp.pi).evalf(),Omega:(40*2*syp.pi).evalf(),rho0:1200,rho1:1200,c0:2800,c1:2800,Delta:0.9,delta:0.9,M0:10000,K0:2.45e9}


eqsNumS=list()
for i in range(2*n+1):    
    hV2=syp.poly(polVCoeffs[i+1],h).all_coeffs()
    hS2=syp.poly(polSCoeffs[i+1],h).all_coeffs()
    eqsNumS.append(syp.Eq((hV2[0].expand()/omega/syp.I).subs(paramS).simplify(),0))
    eqsNumS.append(syp.Eq((hS2[0]/omega/syp.I).subs(paramS).simplify(),0))
    #eqsNum.append(syp.Eq((hV2[2].expand()/omega/syp.I).subs(param).simplify(),0))
    #eqsNum.append(syp.Eq((hS2[2]/omega/syp.I).subs(c1,c0).subs(param).simplify(),0))

coeff_matrixNumS, constantsNumS = syp.linear_eq_to_matrix(eqsNumS,list(matCoeffs))


coeff_matrixNumS, constantsNumS=NormalizationParam(coeff_matrixNumS, constantsNumS)
RTNumS=syp.linsolve((coeff_matrixNumS,constantsNumS),list(matCoeffs))



coefsRS=[]
coefsRNormS=[1,1]
coefsRS.append((RTNumS.args[0][2*n]))
#coefsR.append(syp.Abs(RTNum.args[0][2*n+1]))
for i in range(1,n+1):
    coefsRS.append((RTNumS.args[0][2*(n-i)]))#/coefsR[0])
    #coefsR.append(syp.Abs(RTNum.args[0][4*(n-i)+1]))
    coefsRNormS.append((RTNumS.args[0][2*(n-i)])/coefsRS[0])
    #coefsRNorm.append(syp.Abs(RTNum.args[0][4*(n-i)+1])/coefsR[0])
    coefsRS.append((RTNumS.args[0][2*(n+i)]))#/coefsR[0])
    #coefsR.append(syp.Abs(RTNum.args[0][4*(n+i)+1]))#/coefsR[0])
    coefsRNormS.append((RTNumS.args[0][2*(n+i)])/coefsRS[0])
    #coefsRNorm.append(syp.Abs(RTNum.args[0][4*(n+i)+1])/coefsR[0])
    
    
coefsTS=[]
coefsTNormS=[1.]
coefsTS.append(RTNumS.args[0][2*n+1])
#coefsT.append(RTNum.args[0][4*n+3])
for i in range(1,n+1):
    coefsTS.append(RTNumS.args[0][2*(n-i)+1].evalf())#/coefsR[0])
    #coefsT.append((RTNum.args[0][4*(n-i)+3]))
    coefsTNormS.append((RTNumS.args[0][2*(n-i)+1])/coefsTS[0])
    #coefsTNorm.append((RTNum.args[0][4*(n-i)+3])/coefsT[0])
    coefsTS.append((RTNumS.args[0][2*(n+i)+1]))#/coefsR[0])
    #coefsT.append((RTNum.args[0][4*(n+i)+3]))#/coefsR[0])
    coefsTNormS.append((RTNumS.args[0][2*(n+i)+1])/coefsTS[0])
    #coefsTNorm.append((RTNum.args[0][4*(n+i)+3])/coefsT[0])
coefsTrS=[]
coefsTrNormS=[1.]
coefsTrS.append(syp.re(RTNumS.args[0][2*n+1]))
#coefsT.append(RTNum.args[0][4*n+3])
for i in range(1,n+1):
    coefsTrS.append(syp.re(RTNumS.args[0][2*(n-i)+1].evalf()))#/coefsR[0])
    #coefsT.append((RTNum.args[0][4*(n-i)+3]))
    coefsTrNormS.append(syp.re((RTNumS.args[0][2*(n-i)+1])/coefsTrS[0]))
    #coefsTNorm.append((RTNum.args[0][4*(n-i)+3])/coefsT[0])
    coefsTrS.append(syp.re((RTNumS.args[0][2*(n+i)+1])))#/coefsR[0])
    #coefsT.append((RTNum.args[0][4*(n+i)+3]))#/coefsR[0])
    coefsTrNormS.append(syp.re((RTNumS.args[0][2*(n+i)+1])/coefsTrS[0]))
    #coefsTNorm.append((RTNum.args[0][4*(n+i)+3])/coefsT[0])
coefsTiS=[]
coefsTiNormS=[1.]
coefsTiS.append(RTNumS.args[0][2*n+1])
#coefsT.append(RTNum.args[0][4*n+3])
for i in range(1,n+1):
    coefsTiS.append(syp.im(RTNumS.args[0][2*(n-i)+1].evalf()))#/coefsR[0])
    #coefsT.append((RTNum.args[0][4*(n-i)+3]))
    coefsTiNormS.append(syp.im((RTNumS.args[0][2*(n-i)+1])/coefsTS[0]))
    #coefsTNorm.append((RTNum.args[0][4*(n-i)+3])/coefsT[0])
    coefsTiS.append(syp.im((RTNumS.args[0][2*(n+i)+1])))#/coefsR[0])
    #coefsT.append((RTNum.args[0][4*(n+i)+3]))#/coefsR[0])
    coefsTiNormS.append(syp.im((RTNumS.args[0][2*(n+i)+1])/coefsTS[0]))
    #coefsTNorm.append((RTNum.args[0][4*(n+i)+3])/coefsT[0])
coefsTaS=[]
coefsTaNormS=[1.]
coefsTaS.append(RTNumS.args[0][2*n+1])
#coefsT.append(RTNum.args[0][4*n+3])
for i in range(1,n+1):
    coefsTaS.append(syp.Abs(RTNumS.args[0][2*(n-i)+1].evalf()))#/coefsR[0])
    #coefsT.append((RTNum.args[0][4*(n-i)+3]))
    coefsTaNormS.append(syp.Abs((RTNumS.args[0][2*(n-i)+1])/coefsTaS[0]))
    #coefsTNorm.append((RTNum.args[0][4*(n-i)+3])/coefsT[0])
    coefsTaS.append(syp.Abs((RTNumS.args[0][2*(n+i)+1])))#/coefsR[0])
    #coefsT.append((RTNum.args[0][4*(n+i)+3]))#/coefsR[0])
    coefsTaNormS.append(syp.Abs((RTNumS.args[0][2*(n+i)+1])/coefsTaS[0]))
    #coefsTNorm.append((RTNum.args[0][4*(n+i)+3])/coefsT[0])

x0=syp.symbols("x0",nonzero=True)
VTnOS=(vTn).subs({T1:syp.exp(syp.I*(omega)/c0*(x-x0))*coefsT[0]}).subs(paramS)
#VTnm=vTn.subs({omega:-200*2*syp.pi.evalf(),Omega:-17*2*syp.pi.evalf(),c0:2800,c1:2800,T1:0})

### Dephasage dependant de Omega
# for j in range(n):
#     #VTnm=VTnm.subs(syp.Indexed('Tm', j+1),coefsT[2*j+1]).subs(syp.Indexed('Tp', j+1),coefsT[2*j+2])
#     VTnO=VTnO.subs(param).subs(syp.Indexed('Tm', j+1),syp.exp(syp.I*(omega-j*Omega)/c0*(x-x0))*coefsT[4*j+2]).subs(syp.Indexed('Tp', j+1),syp.exp(syp.I*(omega+j*Omega)/c0*(x-x0))*coefsT[4*j+4])
#     VTnO=VTnO.subs(param).subs(syp.Indexed('tm', j+1),syp.exp(-syp.I*(omega-j*Omega)/c0*(x-x0))*coefsT[4*j+3]).subs(syp.Indexed('tp', j+1),syp.exp(-syp.I*(omega+j*Omega)/c0*(x-x0))*coefsT[4*j+5])
# Nt=4096*2
### Dephasage independant de Omega
for j in range(n):
    VTnOS=VTnOS.subs(syp.Indexed('Tm', j+1),coefsTS[2*j+1]).subs(syp.Indexed('Tp', j+1),coefsTS[2*j+2])
    #VTnO=VTnO.subs(param).subs(syp.Indexed('Tm', j+1),syp.exp(syp.I*(omega)/c0*(x-x0))*coefsT[4*j+2]).subs(syp.Indexed('Tp', j+1),syp.exp(syp.I*(omega)/c0*(x-x0))*coefsT[4*j+4])
    #VTnO=VTnO.subs(param).subs(syp.Indexed('tm', j+1),syp.exp(-syp.I*(omega)/c0*(x-x0))*coefsT[4*j+3]).subs(syp.Indexed('tp', j+1),syp.exp(-syp.I*(omega)/c0*(x-x0))*coefsT[4*j+5])
Nt=4096*2

VTnOS=(VTnOS/c0/2/omega).subs({x0:100}).subs(paramS)
syp.plot(syp.re(VTnOS).subs(x,100), (t, 0, Nt*dt),adaptive=False,nb_of_points=Nt)






# eqsNum=list()
# for i in range(2*n+1):
#     hV2=syp.poly(polVCoeffs[i+1],h).all_coeffs()
#     hS2=syp.poly(polSCoeffs[i+1],h).all_coeffs()
#     eqsNum.append(syp.Eq((hV2[0].expand()/omega/syp.I).subs(param).simplify(),0))
#     eqsNum.append(syp.Eq((hS2[0]/omega/syp.I).subs(param).simplify(),0))
#     eqsNum.append(syp.Eq((hV2[2].expand()/omega/syp.I).subs(param).simplify(),0))
#     eqsNum.append(syp.Eq((hS2[2]/omega/syp.I).subs(param).simplify(),0))
#     # eqsNum.append(syp.Eq((polVCoeffs[i+1].expand()/h/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(param).simplify(),0))
#     # eqsNum.append(syp.Eq((polSCoeffs[i+1].expand()/h/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(param).simplify(),0))

# coeff_matrixNum, constantsNum = syp.linear_eq_to_matrix(eqsNum,list(matCoeffs))
# RTNum=syp.linsolve(eqsNum,list(matCoeffs))


# coefsR=[]
# coefsRNorm=[1]
# coefsR.append(syp.Abs(RTNum.args[0][4*n]))
# coefsR.append(syp.Abs(RTNum.args[0][4*n+1]))
# for i in range(1,n+1):
#     coefsR.append(syp.Abs(RTNum.args[0][4*(n-i)]))#/coefsR[0])
#     coefsR.append(syp.Abs(RTNum.args[0][4*(n-i)+1]))
#     coefsRNorm.append(syp.Abs(RTNum.args[0][4*(n-i)])/coefsR[0])
#     coefsRNorm.append(syp.Abs(RTNum.args[0][4*(n-i)+1])/coefsR[0])
#     coefsR.append(syp.Abs(RTNum.args[0][4*(n+i)]))#/coefsR[0])
#     coefsR.append(syp.Abs(RTNum.args[0][4*(n+i)+1]))#/coefsR[0])
#     coefsRNorm.append(syp.Abs(RTNum.args[0][4*(n+i)])/coefsR[0])
#     coefsRNorm.append(syp.Abs(RTNum.args[0][4*(n+i)+1])/coefsR[0])
    
    
# coefsT=[]
# coefsTNorm=[1]
# coefsT.append(syp.Abs(RTNum.args[0][4*n+2]))
# coefsT.append(syp.Abs(RTNum.args[0][4*n+3]))
# for i in range(1,n+1):
#     coefsT.append(syp.Abs(RTNum.args[0][4*(n-i)+2]))#/coefsR[0])
#     coefsT.append(syp.Abs(RTNum.args[0][4*(n-i)+3]))
#     coefsTNorm.append(syp.Abs(RTNum.args[0][4*(n-i)+2])/coefsT[0])
#     coefsTNorm.append(syp.Abs(RTNum.args[0][4*(n-i)+3])/coefsT[0])
#     coefsT.append(syp.Abs(RTNum.args[0][4*(n+i)+2]))#/coefsR[0])
#     coefsT.append(syp.Abs(RTNum.args[0][4*(n+i)+3]))#/coefsR[0])
#     coefsTNorm.append(syp.Abs(RTNum.args[0][4*(n+i)+2])/coefsT[0])
#     coefsTNorm.append(syp.Abs(RTNum.args[0][4*(n+i)+3])/coefsT[0])

# x0=syp.symbols("x0",nonzero=True)
# VTnO=(vTn).subs({omega:200*2*syp.pi.evalf(),Omega:35*2*syp.pi.evalf(),c0:2800,c1:2800,T1:syp.exp(syp.I*(omega)/c0*(x-x0))*coefsT[0],t1:syp.exp(-syp.I*(omega)/c0*(x-x0))*coefsT[1]})
# #VTnm=vTn.subs({omega:-200*2*syp.pi.evalf(),Omega:-17*2*syp.pi.evalf(),c0:2800,c1:2800,T1:0})

# for j in range(n):
#     #VTnm=VTnm.subs(syp.Indexed('Tm', j+1),coefsT[2*j+1]).subs(syp.Indexed('Tp', j+1),coefsT[2*j+2])
#     VTnO=VTnO.subs(syp.Indexed('Tm', j+1),syp.exp(syp.I*(omega-j*Omega)/c0*(x-x0))*coefsT[4*j+2]).subs(syp.Indexed('Tp', j+1),syp.exp(syp.I*(omega+j*Omega)/c0*(x-x0))*coefsT[4*j+3])
#     VTnO=VTnO.subs(syp.Indexed('tm', j+1),syp.exp(-syp.I*(omega-j*Omega)/c0*(x-x0))*coefsT[4*j+4]).subs(syp.Indexed('tp', j+1),syp.exp(-syp.I*(omega+j*Omega)/c0*(x-x0))*coefsT[4*j+5])
# Nt=4096

# VTnO=(VTnO/c0/2/omega).subs({omega:200*2*syp.pi.evalf(),Omega:35*2*syp.pi.evalf(),c0:2800,x0:-75})
# syp.plot(syp.re(VTnO).subs(x,100), (t, 0, Nt*dt),adaptive=False,nb_of_points=Nt)


# # import numpy as np
# invTF=syp.zeros(Nx,Nt)
# for xi in range(Nx):
#     for ti in range(Nt):
#         Vinv=(VTnO+VTnm)
#         val=Vinv.subs(x,xi*200/Nx).subs(t,ti*dt).evalf()
#         invTF[xi,ti]=val
        
