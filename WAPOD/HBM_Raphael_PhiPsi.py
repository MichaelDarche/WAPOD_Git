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
# R0,R1,R1p,R1m=syp.symbols("R0 R1 R1p R1m")
# T0,T1,T1p,T1m=syp.symbols("T0 T1 T1p T1m")
############ Constants ####################
Z0=rho0*c0
Z1=rho1*c1
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
uI=syp.exp(-syp.I*(omega*(t-x/c0)))
sI=Z0*c0*syp.diff(uI,x)
#%%%%%%%% Order n %%%%%%%%%%%%%%%%%%%
def somme_variables_symboliques(n):
    # Declarer une variable symbolique i
    l = syp.symbols('l')
    
    
    
    # Creer une somme symbolique avec 2n+1 termes
    uRpn = syp.Sum(syp.Indexed('Rp', l)*syp.exp(-syp.I*(omega+l*Omega)*(t+x/c0)), (l, -n, n))    
    sRpn=Z0*c0*syp.diff(uRpn,x)    

    uTpn = syp.Sum(syp.Indexed('Tp', l)*syp.exp(-syp.I*(omega+l*Omega)*(t-x/c1)), (l, -n, n))
    sTpn=Z1*c1*syp.diff(uTpn,x)
    # sTpn = -syp.Sum(Z1*syp.I*(omega+l*Omega)*syp.Indexed('Tp', l)*syp.exp(syp.I*(omega+l*Omega)*(t-x/c1)), (l, 1, n))
    # sTmn = -syp.Sum(Z1*syp.I*(omega-l*Omega)*syp.Indexed('Tm', l)*syp.exp(syp.I*(omega-l*Omega)*(t-x/c1)), (l, 1, n))
    # sTn=-Z1*syp.I*omega*T1*syp.exp(syp.I*(omega*(t-x/c1)))+sTpn+sTmn
    # Simplifier l'expression
    matCoeffs=syp.zeros(2*(2*n+1),1)
    matCoeffsPsi=syp.zeros((2*n+1),1)
    matCoeffsPhi=syp.zeros((2*n+1),1)
    # matCoeffs[2*n]=R1
    # matCoeffs[2*n+1]=T1
    for l in range(-n,n+1):
        matCoeffs[2*(l+n)]=(syp.Indexed('Rp', l))
        matCoeffs[2*(l+n)+1]=(syp.Indexed('Tp', l))
        matCoeffsPhi[l+n]=(syp.Indexed('Phi', l))
        matCoeffsPsi[l+n]=(syp.Indexed('Psi', l))
        # matCoeffs[2*(n-l)]=(syp.Indexed('Rm', l))
        # matCoeffs[2*(n+l)]=(syp.Indexed('Rp', l))
        # matCoeffs[2*(n-l)+1]=(syp.Indexed('Tm', l))
        # matCoeffs[2*(n+l)+1]=(sypw.Indexed('Tp', l))
    return uRpn.simplify(),sRpn.simplify(),uTpn.simplify(),sTpn.simplify(),matCoeffs,matCoeffsPsi,matCoeffsPhi
#% n-Displacement
n=2
uRn,sRn,uTn,sTn,matCoeffs,matCoeffsPsi,matCoeffsPhi=somme_variables_symboliques(n)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vI=syp.diff(uI,t)
dSI=syp.diff(sI,t)

h,H=syp.symbols("h H",nonzero=True)

param={omega:(63*2*syp.pi).evalf(),Omega:(115*2*syp.pi).evalf(),rho0:1200,rho1:1200,c0:2800,c1:2800,Delta:0.9,delta:0.9,M0:20000,K0:2.45e9}
######################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### n-order modulation #############
csV=(uTn-uI-uRn).simplify().subs(x,0).subs(syp.exp(-syp.I*omega*t),h).subs(syp.exp(syp.I*Omega*t),H)
moV=1/2*(uTn+uI+uRn).subs(x,0).simplify().subs(syp.exp(-syp.I*omega*t),h).subs(syp.exp(syp.I*Omega*t),H)

for l in range(-n,n+1):
    csV=(csV).simplify().subs(syp.Indexed('Tp', l)-syp.Indexed('Rp', l),syp.Indexed('Psi', l))
    moV=(moV).simplify().subs(syp.Indexed('Tp', l)+syp.Indexed('Rp', l),syp.Indexed('Phi', l))

csSn=(sTn-sI-sRn).subs(rho1,rho0).subs(c1,c0).subs(x,0).subs(syp.exp(-syp.I*omega*t),h).subs(syp.exp(syp.I*Omega*t),H).expand().simplify()
moSn=1/2*(sTn+sI+sRn).subs(rho1,rho0).subs(c1,c0).subs(x,0).subs(syp.exp(-syp.I*omega*t),h).subs(syp.exp(syp.I*Omega*t),H).expand().simplify()

csS=csSn
moS=moSn
for l in range(-n,n+1):
    csS=syp.collect((csS),(omega)).simplify().subs(syp.Indexed('Tp', l)+syp.Indexed('Rp', l),syp.Indexed('Phi', l))
    csS=syp.collect((csS),(Omega)).simplify().subs(syp.Indexed('Tp', l)+syp.Indexed('Rp', l),syp.Indexed('Phi', l))
    csS=syp.collect((csS),(H)).simplify().subs(syp.Indexed('Tp', l)+syp.Indexed('Rp', l),syp.Indexed('Phi', l))
    moS=syp.collect((moS),(omega)).simplify().subs(syp.Indexed('Tp', l)-syp.Indexed('Rp', l),syp.Indexed('Psi', l))
    moS=syp.collect((moS),(Omega)).simplify().subs(syp.Indexed('Tp', l)-syp.Indexed('Rp', l),syp.Indexed('Psi', l))
    moS=syp.collect((moS),(H)).simplify().subs(syp.Indexed('Tp', l)-syp.Indexed('Rp', l),syp.Indexed('Psi', l))

ccV=(F*moS).subs(x,0).subs(syp.exp(-syp.I*omega*t),h).subs(syp.exp(syp.I*Omega*t),H)
ccS=(dtM*moV+M*syp.diff(moV,t)).subs(x,0).subs(syp.exp(-syp.I*omega*t),h).subs(syp.exp(syp.I*Omega*t),H)

## to simplify
## interface en x=0
PolyHPsi=syp.poly(((csV-ccV)*H**(n+1)).expand(),H)
if PolyHPsi.degree(H)!=2*(n+1):
    PolyHPsi=syp.poly(((csV-ccV)*H**(n+2)).expand(),H)
PolyHPhi=syp.poly(((csS-ccS)*H**(n+1)).expand(),H)
if PolyHPhi.degree(H)!=2*(n+1):
    PolyHPsi=syp.poly(((csS-ccS)*H**(n+2)).expand(),H)

polVCoeffs=PolyHPsi.all_coeffs()
polSCoeffs=PolyHPhi.all_coeffs()



eqsPsi=list()
eqsPhi=list()
for i in range(2*n+1):
    eqsPsi.append(syp.Eq((polVCoeffs[i+1].expand()/h).subs(h,syp.exp(-syp.I*omega*t)).simplify(),0))
    eqsPhi.append(syp.Eq((polSCoeffs[i+1]/h/omega).subs(h,syp.exp(-syp.I*omega*t)).simplify(),0))

coeff_matrixPsi, constantsPsi = syp.linear_eq_to_matrix(eqsPsi,list(matCoeffsPsi))
coeff_matrixPhi, constantsPhi = syp.linear_eq_to_matrix(eqsPhi,list(matCoeffsPhi))



eqsRL=list()
for i in range(2*n+1):
    eqsRL.append(syp.Eq((polVCoeffs[i+1]/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))
    eqsRL.append(syp.Eq((polSCoeffs[i+1]/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).simplify(),0))

coeff_matrixRL, constantsRL = syp.linear_eq_to_matrix(eqsRL,list(matCoeffs))

pRL={c0:2800,rho0:1200,Delta:0/10,M0:20000,omega:30,Omega:100}
eqsRLN=list()
for i in range(2*n+1):
    eqsRLN.append(syp.Eq((polVCoeffs[i+1]/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).subs(pRL).simplify(),0))
    eqsRLN.append(syp.Eq((polSCoeffs[i+1]/omega/syp.I).subs(h,syp.exp(syp.I*omega*t)).subs(c1,c0).subs(rho1,rho0).subs(delta,Delta).subs(K0,c0**2*rho0**2/M0).subs(pRL).simplify(),0))

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
        
