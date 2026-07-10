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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############# Fields ######################
uI=syp.exp(syp.I*(omega*(t-x/c0)))
sI=-syp.I*omega*Z0*syp.exp(syp.I*(omega*(t-x/c0)))
#%%%%%%%% Order 1 %%%%%%%%%%%%%%%%%%%

#%%%%%%%% Order 2 %%%%%%%%%%%%%%%%%%%
#% 2-Displacement
uR2=R1*syp.exp(syp.I*(omega*(t+x/c0)))+R1p*syp.exp(syp.I*(omegap*(t+x/c0)))+R1m*syp.exp(syp.I*(omegam*(t+x/c0)))
uT2=T1*syp.exp(syp.I*(omega*(t-x/c1)))+T1p*syp.exp(syp.I*(omegap*(t-x/c1)))+T1m*syp.exp(syp.I*(omegam*(t-x/c1)))
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
    uRn=R1*syp.exp(syp.I*(omega*(t+x/c0)))+uRpn+uRmn
    
    sRpn = syp.Sum(Z0*syp.I*(omega+l*Omega)*syp.Indexed('Rp', l)*syp.exp(syp.I*(omega+l*Omega)*(t+x/c0)), (l, 1, n))
    sRmn = syp.Sum(Z0*syp.I*(omega-l*Omega)*syp.Indexed('Rm', l)*syp.exp(syp.I*(omega-l*Omega)*(t+x/c0)), (l, 1, n))
    sRn=Z0*syp.I*omega*R1*syp.exp(syp.I*(omega*(t+x/c0)))+sRpn+sRmn
    
    uTpn = syp.Sum(syp.Indexed('Tp', l)*syp.exp(syp.I*(omega+l*Omega)*(t-x/c1)), (l, 1, n))
    uTmn = syp.Sum(syp.Indexed('Tm', l)*syp.exp(syp.I*(omega-l*Omega)*(t-x/c1)), (l, 1, n))
    uTn=T1*syp.exp(syp.I*(omega*(t-x/c1)))+uTpn+uTmn
    
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
dSI=syp.diff(sI,t)

h,H=syp.symbols("h H",nonzero=True)

param={omega:(63*2*syp.pi).evalf(),Omega:(115*2*syp.pi).evalf(),rho0:1200,rho1:1200,c0:2800,c1:2800,Delta:0.9,delta:0.9,M0:10000,K0:2.45e9}
######################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### n-order modulation #############
csVn=(uTn-uI-uRn).subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(-syp.I*Omega*t),H)
moVn=1/2*(uTn+uI+uRn).subs(x,0)
csSn=(sTn-sI-sRn).subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(-syp.I*Omega*t),H)
moSn=1/2*(sTn+sI+sRn).subs(x,0)
ccVn=(F*moSn).subs(x,0).subs(syp.exp(syp.I*omega*t),h).subs(syp.exp(-syp.I*Omega*t),H)
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

pRL={c0:2800,rho0:1200,Delta:9/10,M0:20000,omega:30,Omega:100}
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

VTnO=uTn.subs({omega:200*2*syp.pi.evalf(),Omega:17*2*syp.pi.evalf(),c0:2800,c1:2800,T1:coefsT[0]})
VTnm=uTn.subs({omega:-200*2*syp.pi.evalf(),Omega:-17*2*syp.pi.evalf(),c0:2800,c1:2800,T1:0})

for j in range(n):
    VTnm=VTnm.subs(syp.Indexed('Tm', j+1),coefsT[4*j+1]).subs(syp.Indexed('Tp', j+1),coefsT[2*j+2]).subs(syp.Indexed('tm', j+1),coefsT[2*j+1]).subs(syp.Indexed('tp', j+1),coefsT[2*j+2])
    VTnO=VTnO.subs(syp.Indexed('Tm', j+1),coefsT[4*j+1]).subs(syp.Indexed('Tp', j+1),coefsT[2*j+2])

# Nt=4096




# syp.plot(syp.re(VTnm).subs(x,300),syp.im(VTnm).subs(x,300), (t, 0, Nt*dt),adaptive=False,nb_of_points=Nt)


# # import numpy as np
# invTF=syp.zeros(Nx,Nt)
# for xi in range(Nx):
#     for ti in range(Nt):
#         Vinv=(VTnO+VTnm)
#         val=Vinv.subs(x,xi*200/Nx).subs(t,ti*dt).evalf()
#         invTF[xi,ti]=val
        
