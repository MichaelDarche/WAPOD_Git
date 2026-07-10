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
R0,R1,R1p,R1m=syp.symbols("R0 R1 R1p R1m")
T0,T1,T1p,T1m=syp.symbols("T0 T1 T1p T1m")
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
#%%%%%%%% Order 1 %%%%%%%%%%%%%%%%%%%
#% 1-Displacement
uR1=R0+R1*syp.exp(syp.I*(omega*(t+x/c0)))
uT1=T0+T1*syp.exp(syp.I*(omega*(t-x/c1)))
#% 1-Stress
sR1=+syp.I*omega*Z0*R1*syp.exp(syp.I*(omega*(t+x/c0)))
sT1=-syp.I*omega*Z1*T1*syp.exp(syp.I*(omega*(t-x/c1)))
#%%%%%%%% Order 2 %%%%%%%%%%%%%%%%%%%
#% 2-Displacement
uR2=R0+R1*syp.exp(syp.I*(omega*(t+x/c0)))+R1p*syp.exp(syp.I*(omegap*(t+x/c0)))+R1m*syp.exp(syp.I*(omegam*(t+x/c0)))
uT2=T0+T1*syp.exp(syp.I*(omega*(t-x/c1)))+T1p*syp.exp(syp.I*(omegap*(t-x/c1)))+T1m*syp.exp(syp.I*(omegam*(t-x/c1)))
#% 2-Stress
sR2=+syp.I*omega*Z0*R1*syp.exp(syp.I*(omega*(t+x/c0)))+syp.I*omegap*Z0*R1p*syp.exp(syp.I*(omegap*(t+x/c0)))+syp.I*omegam*Z0*R1m*syp.exp(syp.I*(omegam*(t+x/c0)))
sT2=-syp.I*omega*Z1*T1*syp.exp(syp.I*(omega*(t-x/c1)))-syp.I*omegap*Z1*T1p*syp.exp(syp.I*(omegap*(t-x/c1)))-syp.I*omegam*Z1*T1m*syp.exp(syp.I*(omegam*(t-x/c1)))
#%%%%%%%% Order n %%%%%%%%%%%%%%%%%%%
def somme_variables_symboliques(n):
    # Declarer une variable symbolique i
    l = syp.symbols('l')

    # Creer une somme symbolique avec 2n+1 termes
    uRpn = syp.Sum(syp.Indexed('Rp', l)*syp.exp(syp.I*(omega+l*Omega)*(t+x/c0)), (l, 1, n))
    uRmn = syp.Sum(syp.Indexed('Rm', l)*syp.exp(syp.I*(omega-l*Omega)*(t+x/c0)), (l, 1, n))
    uRn=R0+R1*syp.exp(syp.I*(omega*(t+x/c0)))+uRpn+uRmn
    
    sRpn = syp.Sum(Z0*syp.I*(omega+l*Omega)*syp.Indexed('Rp', l)*syp.exp(syp.I*(omega+l*Omega)*(t+x/c0)), (l, 1, n))
    sRmn = syp.Sum(Z0*syp.I*(omega-l*Omega)*syp.Indexed('Rm', l)*syp.exp(syp.I*(omega-l*Omega)*(t+x/c0)), (l, 1, n))
    sRn=Z0*syp.I*omega*R1*syp.exp(syp.I*(omega*(t+x/c0)))+sRpn+sRmn
    
    uTpn = syp.Sum(syp.Indexed('Tp', l)*syp.exp(syp.I*(omega+l*Omega)*(t-x/c1)), (l, 1, n))
    uTmn = syp.Sum(syp.Indexed('Tm', l)*syp.exp(syp.I*(omega-l*Omega)*(t-x/c1)), (l, 1, n))
    uTn=T0+T1*syp.exp(syp.I*(omega*(t-x/c1)))+uTpn+uTmn
    
    sTpn = -syp.Sum(Z1*syp.I*(omega+l*Omega)*syp.Indexed('Tp', l)*syp.exp(syp.I*(omega+l*Omega)*(t-x/c1)), (l, 1, n))
    sTmn = -syp.Sum(Z1*syp.I*(omega-l*Omega)*syp.Indexed('Tm', l)*syp.exp(syp.I*(omega-l*Omega)*(t-x/c1)), (l, 1, n))
    sTn=-Z1*syp.I*omega*T1*syp.exp(syp.I*(omega*(t-x/c1)))+sTpn+sTmn
    # Simplifier l'expression
    matCoeffs=syp.zeros(2*(2*n+1),1)
    matCoeffs[2*n]=R1
    matCoeffs[2*n+1]=T1
    for l in range(1,n+1):
        matCoeffs[2*(n-l)]=(syp.Indexed('Rm', l))
        matCoeffs[2*(n+l)]=(syp.Indexed('Rp', l))
        matCoeffs[2*(n-l)+1]=(syp.Indexed('Tm', l))
        matCoeffs[2*(n+l)+1]=(syp.Indexed('Tp', l))
    return uRn.simplify(),sRn.simplify(),uTn.simplify(),sTn.simplify(),matCoeffs
#% n-Displacement
n=10
uRn,sRn,uTn,sTn,matCoeffs=somme_variables_symboliques(n)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
vRn=syp.diff(uRn,t)
vTn=syp.diff(uTn,t)
dSRn=syp.diff(sRn,t)
dSTn=syp.diff(sTn,t)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### 1-order modulation #############
csV1=vT1-vI-vR1
moV1=1/2*(vT1+vI+vR1)
csS1=sT1-sI-sR1
moS1=1/2*(sT1+sI+sR1)
ccV1=dtF*moS1+F*syp.diff(moS1,t)
ccS1=dtM*moV1+M*syp.diff(moV1,t)
## to simplify
## interface en x=0
h,H=syp.symbols("h H",nonzero=True)
eqV1=syp.factor(syp.Eq(csV1,ccV1).subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(syp.I*Omega*t),H))
eqS1=syp.factor(syp.Eq(csS1,ccS1).subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(syp.I*Omega*t),H))
coeff_matrix, constants = syp.linear_eq_to_matrix((eqV1, eqS1), (R1, T1))
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
#%%%%%%%%%% Analytique %%%%%%%%%%%%%%%
# exAV=(1+syp.I*omega*Z0/(2*K0))*R1-(1+syp.I*omega*Z1/(2*K0))*T1+1-syp.I*omega*Z0/(2*K0)
# eqVA1=syp.Eq(exAV,0).subs(x,0)
# exAS=(Z0+syp.I*omega*M0/(2))*R1-(Z1+syp.I*omega*M0/(2))*T1-Z0+syp.I*omega*M0/(2)
# eqSA1=syp.Eq(exAS,0).subs(x,0)
# coeff_matrixA, constantsA = syp.linear_eq_to_matrix((eqVA1, eqSA1), (R1, T1))
######################################
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

param={omega:(200*2*syp.pi).evalf(),Omega:(35*2*syp.pi).evalf(),rho0:1200,rho1:1200,c0:2800,c1:2800,Delta:0.0,delta:0.0,M0:20000,K0:2.45e9}
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

######################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### n-order modulation #############
csVn=(vTn-vI-vRn).subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(-syp.I*Omega*t),H)
moVn=1/2*(vTn+vI+vRn).subs(x,0)
csSn=(sTn-sI-sRn).subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(-syp.I*Omega*t),H)
moSn=1/2*(sTn+sI+sRn).subs(x,0)
ccVn=(dtF*moSn+F*syp.diff(moSn,t)).subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(-syp.I*Omega*t),H)
ccSn=(dtM*moVn+M*syp.diff(moVn,t)).subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(-syp.I*Omega*t),H)
## to simplify
## interface en x=0

PolyHV=syp.poly(((csVn-ccVn)*H**(n+1)).expand(),H)
if PolyHV.degree(H)!=2*(n+1):
    PolyHV=syp.poly(((csVn-ccVn)*H**(n+2)).expand(),H)
PolyHS=syp.poly(((csSn-ccSn)*H**(n+1)).expand(),H)
if PolyHS.degree(H)!=2*(n+1):
    PolyHS=syp.poly(((csSn-ccSn)*H**(n+2)).expand(),H)
    
polVCoeffs=PolyHV.all_coeffs()
polSCoeffs=PolyHS.all_coeffs()
eqs=list()
for i in range(2*n+1):
    eqs.append(syp.Eq((polVCoeffs[i+1].expand()/h/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).simplify(),0))
    eqs.append(syp.Eq((polSCoeffs[i+1]/h/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).simplify(),0))

coeff_matrixN, constantsN = syp.linear_eq_to_matrix(eqs,list(matCoeffs))

eqsRL=list()
for i in range(2*n+1):
    eqsRL.append(syp.Eq((polVCoeffs[i+1]/h/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))
    eqsRL.append(syp.Eq((polSCoeffs[i+1]/h/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))

coeff_matrixRL, constantsRL = syp.linear_eq_to_matrix(eqsRL,list(matCoeffs))

pRL={c0:2800,rho0:1200,Delta:0/10,M0:20000,omega:30,Omega:100}
eqsRLN=list()
for i in range(2*n+1):
    eqsRLN.append(syp.Eq((polVCoeffs[i+1]/h/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).subs(pRL).simplify(),0))
    eqsRLN.append(syp.Eq((polSCoeffs[i+1]/h/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).subs(pRL).simplify(),0))

coeff_matrixRLN, constantsRLN = syp.linear_eq_to_matrix(eqsRLN,list(matCoeffs))
RTRLNum=syp.linsolve(eqsRLN,list(matCoeffs))

eqsNum=list()
for i in range(2*n+1):
    eqsNum.append(syp.Eq((polVCoeffs[i+1].expand()/h/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(param).simplify(),0))
    eqsNum.append(syp.Eq((polSCoeffs[i+1].expand()/h/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(param).simplify(),0))

coeff_matrixNum, constantsNum = syp.linear_eq_to_matrix(eqsNum,list(matCoeffs))
RTNum=syp.linsolve(eqsNum,list(matCoeffs))

coefsR=[]
coefsRNorm=[1]
coefsR.append(syp.Abs(RTNum.args[0][2*n]))
for i in range(1,n+1):
    coefsR.append(syp.Abs(RTNum.args[0][2*(n-i)]))#/coefsR[0])
    coefsRNorm.append(syp.Abs(RTNum.args[0][2*(n-i)])/coefsR[0])
    coefsR.append(syp.Abs(RTNum.args[0][2*(n+i)]))#/coefsR[0])
    coefsRNorm.append(syp.Abs(RTNum.args[0][2*(n+i)])/coefsR[0])
    
coefsT=[]
coefsTNorm=[1]
coefsT.append(syp.Abs(RTNum.args[0][2*n+1]))
for i in range(1,n+1):
    coefsT.append(syp.Abs(RTNum.args[0][2*(n-i)+1]))
    coefsTNorm.append(syp.Abs(RTNum.args[0][2*(n-i)+1])/coefsT[0])
    coefsT.append(syp.Abs(RTNum.args[0][2*(n+i)+1]))
    coefsTNorm.append(syp.Abs(RTNum.args[0][2*(n+i)+1])/coefsT[0])

VTnO=vTn.subs({omega:200*2*syp.pi.evalf(),Omega:17*2*syp.pi.evalf(),c0:2800,c1:2800,T1:coefsT[0]})
VTnm=vTn.subs({omega:-200*2*syp.pi.evalf(),Omega:-17*2*syp.pi.evalf(),c0:2800,c1:2800,T1:0})

for j in range(n):
    VTnm=VTnm.subs(syp.Indexed('Tm', j+1),coefsT[4*j+1]).subs(syp.Indexed('Tp', j+1),coefsT[2*j+2]).subs(syp.Indexed('tm', j+1),coefsT[2*j+1]).subs(syp.Indexed('tp', j+1),coefsT[2*j+2])
    VTnO=VTnO.subs(syp.Indexed('Tm', j+1),coefsT[4*j+1]).subs(syp.Indexed('Tp', j+1),coefsT[2*j+2])

Nt=4096


syp.plot(syp.re(VTnm).subs(x,300),syp.im(VTnm).subs(x,300), (t, 0, Nt*dt),adaptive=False,nb_of_points=Nt)


# # import numpy as np
# invTF=syp.zeros(Nx,Nt)
# for xi in range(Nx):
#     for ti in range(Nt):
#         Vinv=(VTnO+VTnm)
#         val=Vinv.subs(x,xi*200/Nx).subs(t,ti*dt).evalf()
#         invTF[xi,ti]=val
        
