#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 17:01:49 2024

@author: michael
"""
import matplotlib.pyplot as plt
import sympy as syp
import numpy as np
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
############# Symbols ###################
# t,x=syp.symbols("t x")
# omega=syp.symbols("omega",nonzero=True)
# K0,M0,Omega=syp.symbols("K0 M0 Omega", positive=True)
# delta,Delta=syp.symbols("delta Delta",real=True)
# rho0,c0=syp.symbols("rho0 c0", positive=True)
# rho1,c1=syp.symbols("rho1 c1", positive=True)
# R0,R1,R1p,R1m=syp.symbols("R0 R1 R1p R1m")
# T0,T1,T1p,T1m=syp.symbols("T0 T1 T1p T1m")
############ Constants ####################
# Z0=rho0*c0
# Z1=rho1*c1
# # omegap=syp.Abs(omega+Omega)
# # omegam=syp.Abs(omega-Omega)
# ############# Variables ####################
# M=M0*(1+Delta/(2*syp.I)*(syp.exp(syp.I*Omega*t)-syp.exp(-syp.I*Omega*t)))
# F=1/K0*(1+delta/(2*syp.I)*(syp.exp(syp.I*Omega*t)-syp.exp(-syp.I*Omega*t)))

# # M=M0*(1+Delta*(syp.sin(Omega*t)))
# # F=1/K0*(1+delta*(syp.sin(Omega*t)))
# ########## Derivatives #####################
# dtM=syp.diff(M,t)
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ############# Fields ######################
# uI=syp.exp(-syp.I*(omega*(t-x/c0)))
# sI=Z0*c0*syp.diff(uI,x)
#%%%%%%%% Order n %%%%%%%%%%%%%%%%%%%
# def somme_variables_symboliques(n):
#     # Declarer une variable symbolique i
#     l = syp.symbols('l')
    
    
    
#     # Creer une somme symbolique avec 2n+1 termes
#     uRpn = syp.Sum(syp.Indexed('Rp', l)*syp.exp(syp.I*(omega+l*Omega)*(t+x/c0))+syp.Indexed('Rm', l)*syp.exp(-syp.I*(omega+l*Omega)*(t+x/c0)), (l, -n, n))    
#     sRpn=Z0*c0*syp.diff(uRpn,x)    

#     uTpn = syp.Sum(syp.Indexed('Tp', l)*syp.exp(syp.I*(omega+l*Omega)*(t-x/c1))+syp.Indexed('Tm', l)*syp.exp(-syp.I*(omega+l*Omega)*(t-x/c1)), (l, -n, n))
#     sTpn=Z1*c1*syp.diff(uTpn,x)
#     # sTpn = -syp.Sum(Z1*syp.I*(omega+l*Omega)*syp.Indexed('Tp', l)*syp.exp(syp.I*(omega+l*Omega)*(t-x/c1)), (l, 1, n))
#     # sTmn = -syp.Sum(Z1*syp.I*(omega-l*Omega)*syp.Indexed('Tm', l)*syp.exp(syp.I*(omega-l*Omega)*(t-x/c1)), (l, 1, n))
#     # sTn=-Z1*syp.I*omega*T1*syp.exp(syp.I*(omega*(t-x/c1)))+sTpn+sTmn
#     # Simplifier l'expression
#     matCoeffs=syp.zeros(2*(2*n+1),1)
#     matCoeffsPsi=syp.zeros((2*n+1),1)
#     matCoeffsPhi=syp.zeros((2*n+1),1)
#     # matCoeffs[2*n]=R1
#     # matCoeffs[2*n+1]=T1
#     for l in range(-n,n+1):
#         matCoeffs[2*(l+n)]=(syp.Indexed('Rp', l))
#         matCoeffs[2*(l+n)+1]=(syp.Indexed('Tp', l))
#         matCoeffsPhi[l+n]=(syp.Indexed('Phi', l))
#         matCoeffsPsi[l+n]=(syp.Indexed('Psi', l))
#         # matCoeffs[2*(n-l)]=(syp.Indexed('Rm', l))
#         # matCoeffs[2*(n+l)]=(syp.Indexed('Rp', l))
#         # matCoeffs[2*(n-l)+1]=(syp.Indexed('Tm', l))
#         # matCoeffs[2*(n+l)+1]=(sypw.Indexed('Tp', l))
#     return uRpn.simplify(),sRpn.simplify(),uTpn.simplify(),sTpn.simplify(),matCoeffs,matCoeffsPsi,matCoeffsPhi
#% n-Displacement
def signalNum(x,t,n,omega,Omega,rho1,c1,T):
    valTp=0.+0.j
    valTm=0.+0.j
    for k in range(-n,n+1):
        valTp=valTp+T[k+n]*np.exp(1j*(omega+k*Omega)*(t-x/c1))
        valTm=valTm+T[k+2*n+1+n]*np.exp(-1j*(omega+k*Omega)*(t-x/c1))

   
    return valTp,valTm#uRpn.simplify(),sRpn.simplify(),uTpn.simplify(),sTpn.simplify(),matCoeffs,matCoeffsPsi,matCoeffsPhi
n=12

# param={omega:(63*2*syp.pi).evalf(),Omega:(115*2*syp.pi).evalf(),rho0:1200,rho1:1200,c0:2800,c1:2800,Delta:0.9,delta:0.9,M0:10000,K0:2.45e9}
# uRn,sRn,uTn,sTn,matCoeffs,matCoeffsPsi,matCoeffsPhi=somme_variables_symboliques(n)

tr=np.random.rand(2*n+1)
phi=np.pi/2*np.random.rand(2*n+1)
T=np.zeros(2*(2*n+1),dtype="complex_")
for i in range(-n,n+1):
    T[i]=tr[i]*np.exp(1j*phi[i])
    T[i+2*n+1]=tr[i]*np.exp(-1j*phi[i])


x=100
f=123
F=350

omega=(f*2*np.pi)
Omega=(F*2*np.pi)
rho1=1200
c1=2800

Nt=2**14
dt=6.785714285714286e-05
sigT=np.zeros([int(Nt)],dtype='complex_')
for k in range(Nt):
    sigTp,sigTm=signalNum(x,k*dt,n,omega,Omega,rho1,c1,T)
    sigT[k]=1/(2.j)*(sigTp-sigTm)
    
    
def cosine_taper(signal, alpha=0.00,zp=20):
    zeropadding=2**zp
    Size=signal.shape[0]
    M=np.floor((Size*alpha)/2+0.5);
    tapered=np.zeros(zeropadding)
    for j in range(int(Size)):
        if j<=M+1:
            tapered[j]=signal[j] * (0.5 * ( 1-np.cos(j*np.pi/(M+1))))
        elif j<Size - M-1:
            tapered[j]=signal[j]
        elif j<=Size:
            tapered[j]=signal[j] * (0.5 * (1-np.cos((Size-j)*np.pi/(M+1))))
    return tapered
N=18
signal=cosine_taper(sigT,zp=N)
# FFTsig=np.zeros(2**N)
#FFTMapSft=np.zeros([int(N),Nx])
FFTSig=1/Nt*np.abs(np.fft.fftshift(np.fft.fft(signal,n=int(2**N))))
freq = np.fft.fftshift(np.fft.fftfreq(2**N,d=dt))
freqp=[]
freqm=[]
for i in range(-n,n+1):
    freqp.append(f+i*F)
    freqm.append(-(f+i*F))
plt.figure()
plt.plot(freq,FFTSig)
plt.scatter(freqp,np.abs((1/2j)*T[0:2*n+1]),c='red')
plt.scatter(freqm,np.abs((1/2j)*T[2*n+1:]),c='purple')


###############################################################
sigSyn=np.zeros(N)
sigSyn2=np.zeros(N)
sigSyn3=np.zeros(N)
cRand=np.random.rand(2*n+1)
A=freqp[0]
B=12/5*freqp[5]
for nt in range(Nt):
    sigSyn[nt]=np.sin(2*np.pi*A*nt*dt)*np.sin(2*np.pi*B*nt*dt)
    sigSyn2[nt]=-1/4*(np.exp(2j*np.pi*(A+B)*nt*dt)+np.exp(-2j*np.pi*(A+B)*nt*dt)-(np.exp(2j*np.pi*(A-B)*nt*dt)+np.exp(-2j*np.pi*(A-B)*nt*dt)))
    sigSyn3[nt]=-1/2*(np.cos(2*np.pi*(A+B)*nt*dt)-np.cos(2*np.pi*(A-B)*nt*dt))
    # for i in range(2*n+1):
    #     sigSyn[nt]=sigSyn[nt]+cRand[i]*np.exp(-1j*2*np.pi*frp[i]*nt*dt)
# coefs=[float(Norm.evalf()),np.abs(float(syp.Abs(Normm.evalf())),float(Normp.evalf()))]
fftSyn=1/Nt*(np.fft.fftshift(np.abs(np.fft.fft(sigSyn,n=N))))
fftSyn2=1/Nt*(np.fft.fftshift(np.abs(np.fft.fft(sigSyn2,n=N))))
fftSyn3=1/Nt*(np.fft.fftshift(np.abs(np.fft.fft(sigSyn3,n=N))))

# freq = np.fft.fftshift(np.fft.fftfreq(N,d=dt))
# plt.figure()
# plt.plot(tim,sigSyn[0:Nt+1],linewidth=3)
# plt.plot(tim,sigSyn2[0:Nt+1],":",linewidth=2)
# plt.plot(tim,sigSyn3[0:Nt+1],":",linewidth=1)
# plt.figure()
# plt.plot(freq,fftSyn,linewidth=2)
# plt.plot(freq,fftSyn2,":",linewidth=2)
# plt.plot(freq,fftSyn3,":",linewidth=1)
# plt.scatter([-A-B,-A+B,A-B,A+B], [0.25,0.25,0.25,0.25], c='red')
# plt.xlabel('f')
# plt.ylabel('|F(u)|')