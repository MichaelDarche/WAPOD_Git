
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 14:01:07 2024

@author: michael
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as spl
import scipy.sparse.linalg as ssl

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

#Parametres materiaux
# rho0=1200
# c0=2800
# E0=rho0*c0*c0
# h=10
# K0=2.45e9
# M0=2e4
# DK=-0.9
# #DM=0.75
# DM=0.9#(M0*E0+E0*rho0*h)/(M0*E0+K0*M0*h)*DK

# rho=rho0+M0/h
# Delta=DM*M0/h/rho
# E=E0*K0*h/(K0*h+E0)
# delta=DK*E/K0/h

# rho=1000
# Delta=0.0
# E=2.5e9
# delta=0.6

rho=1200
Delta=0.0#0.46875
c=2800
E=rho*c**2
delta=0.75
#Periode spatiale
tau=2*np.pi/20
phi=1/2

deltan=[0.199, 0.328, 0.980, 1.112, 2.292, 2.388]
Omegan=2*np.pi/tau
knth=Omegan/c*np.sqrt(deltan)

#rho0=1200



def fctheavyside(t,tau=tau,phi=phi):
    Tphi=tau*phi
    res=t.copy()
    for i in range(t.shape[0]):
        if t[i]<Tphi:
            res[i]=+1
        else:
            res[i]=-1
    return res    


fctM="np.cos"
#fct=
# fctm="heavy"


def rhot(t,tau=tau,rho=rho,Delta=Delta):
    if fctM=="np.sinplus":
        rhot=rho*(1+Delta*np.sin((2*np.pi*t/tau)))
    elif fctM=="np.sin":
        rhot=rho*(1+Delta*np.sin((2*np.pi*t/tau)))
    if fctM=="np.cos":
        rhot=rho*(1+Delta*np.cos((2*np.pi*t/tau)))
    elif fctM=="np.sign":
        rhot0=np.sin(2*np.pi*t/tau)
        rhot=rho*(1+Delta*np.sign(rhot0))
    elif fctM=="doublesin":
        rhot=rho*(1+Delta/2*(np.sin((2*np.pi*t/tau))+np.sin((4*np.pi*t/tau))))
    elif fctM=="heavy":
        rhot=rho*(1+Delta*fctheavyside(t,tau))
    return rhot
#


def Et(t,tau=tau,E=E,delta=delta):
    if fctM=="np.sinplus":
        Et=E*(1+delta*np.sin((2*np.pi*t/tau)))
    elif fctM=="np.sin":
        Et=1/(1/E*(1+delta*np.sin((2*np.pi*t/tau))))
    elif fctM=="np.cos":
        Et=E*(1+delta*np.cos((2*np.pi*t/tau)))
    elif fctM=="np.sign":
        e0=np.sin((2*np.pi*t/tau))
        Et=1/(1/E*(1+delta*np.sign(e0)))
    elif fctM=="doublesin":
        Et=1/(1/E*(1+delta/2*(np.sin((2*np.pi*t/tau))+np.sin((4*np.pi*t/tau)))))    
    elif fctM=="heavy":
        Et=1/(1/E*(1+delta*fctheavyside(t,tau)))
    return Et


# Series de Fourier
def SF(fct,n,tau,M):
    t=np.linspace(0,tau,M)
    dt=t[1]-t[0]
    an=dt/tau*np.sum((fct(t[1:])+fct(t[:-1]))/2*np.exp(-2j*n*np.pi/tau*(t[1:]+t[:-1])/2))   #(fct(t[:-1])*np.exp(-2j*n*np.pi/tau*(t[:-1]+dt/2)))
    return an



# Bloch vector
Nom=50
om=np.linspace(0,2*np.pi/tau,Nom+1)

#Taille de la SF
N=2**4
#Integration pour les coefs
Ni=500
#Nbre de modes
NE=12


# Initialisation des matrices
MatrixA=np.zeros([2*N+1,2*N+1],dtype=complex)
MatrixB=np.zeros([2*N+1,2*N+1],dtype=complex)
VP=np.zeros([NE,Nom+1],dtype=complex)
VPI=np.zeros([NE,Nom+1],dtype=complex)


# Remplissage des matrices
for omeg in range(Nom+1):
    print(omeg)
    for n in range(-N,N+1):
        nn=n
        for p in range(-N,N+1):
            pp=p
            anE=SF(Et,pp-nn,tau,Ni)
            anR=SF(rhot,pp-nn,tau,Ni)
            MatrixA[p,n]=(2*np.pi*nn/tau-om[omeg])*(2*np.pi*pp/tau-om[omeg])*anR
            MatrixB[p,n]=anE
    EigV=ssl.eigsh(MatrixA,k=NE,M=MatrixB,which='SM',return_eigenvectors=False)
    EigVI=ssl.eigs(-MatrixA,k=NE,M=MatrixB,which='SM',return_eigenvectors=False)
    VP[:,omeg]=np.sqrt(EigV)
    VPI[:,omeg]=np.sqrt(EigVI)
    
recoE=np.zeros(Ni)
recoR=np.zeros(Ni)
t=np.linspace(0,tau,Ni)
dt=t[1]-t[0]
for n in range(-N,N+1):
    nn=n
    #pp=p-N
    SFpn=SF(rhot,nn,tau,Ni)
    SFEn=SF(Et,nn,tau,Ni)
    for ti in range(Ni):
        recoE[ti]=recoE[ti]+SFEn*np.exp(2j*nn*np.pi/tau*(t[ti]))
        recoR[ti]=recoR[ti]+SFpn*np.exp(2j*nn*np.pi/tau*(t[ti]))

VP=np.nan_to_num(VP) 
VP.sort(axis=0)
VPI.sort(axis=0)
# plt.figure()
# for i in range(2*N):
#     plt.plot((VP[i,:]),om)

############## To compare : small number of Fourier
# Bloch vector

#Taille de la SF
N2=4
NE2=N2*2+1


# Initialisation des matrices
MatrixA2=np.zeros([2*N2+1,2*N2+1],dtype=complex)
MatrixB2=np.zeros([2*N2+1,2*N2+1],dtype=complex)
VP2=np.zeros([NE2,Nom+1],dtype=complex)
VPI2=np.zeros([NE2,Nom+1],dtype=complex)


# Remplissage des matrices
for omeg in range(Nom+1):
    print(omeg)
    for n in range(-N2,N2+1):
        nn=n
        for p in range(-N2,N2+1):
            pp=p
            anE=SF(Et,pp-nn,tau,Ni)
            anR=SF(rhot,pp-nn,tau,Ni)
            MatrixA2[p,n]=(2*np.pi*nn/tau-om[omeg])*(2*np.pi*pp/tau-om[omeg])*anR
            MatrixB2[p,n]=anE
    EigV2=ssl.eigsh(MatrixA2,k=NE,M=MatrixB2,which='SM',return_eigenvectors=False)
    EigVI2=ssl.eigs(-MatrixA2,k=NE,M=MatrixB2,which='SM',return_eigenvectors=False)
    VP2[:,omeg]=np.sqrt(EigV2)
    VPI2[:,omeg]=np.sqrt(EigVI2)
    
recoE2=np.zeros(Ni)
recoR2=np.zeros(Ni)
t=np.linspace(0,tau,Ni)
dt=t[1]-t[0]
for n in range(-N2,N2+1):
    nn=n
    #pp=p-N
    SFpn2=SF(rhot,nn,tau,Ni)
    SFEn2=SF(Et,nn,tau,Ni)
    for ti in range(Ni):
        recoE2[ti]=recoE2[ti]+SFEn2*np.exp(2j*nn*np.pi/tau*(t[ti]))
        recoR2[ti]=recoR2[ti]+SFpn2*np.exp(2j*nn*np.pi/tau*(t[ti]))

VP2=np.nan_to_num(VP2) 
VP2.sort(axis=0)
VPI2.sort(axis=0)
##################################
#Comparaison

rho1=rho*(1+Delta)
E1=1/(1/E*(1+delta))

rho2=rho*(1-Delta)
E2=1/(1/E*(1-delta))

c1=np.sqrt(E1/rho1)
c2=np.sqrt(E2/rho2)
cr=(phi)*c1+(1-phi)*c2





plt.figure()
plt.plot(t,rhot(t),label=r'$\rho$ (t)')
plt.plot(t,recoR,label=r'PWE - $N=16$')
plt.plot(t,recoR2,label=r'PWE - $N=3$')
# plt.plot(np.linspace(0,phi*tau,Ni),np.linspace(rho1,rho1,Ni),)
# plt.plot(np.linspace(phi*tau,tau,Ni),np.linspace(rho2,rho2,Ni))#,color='gray',linestyle='--')
plt.xlabel(r'$t$ (s)',fontsize=16)
plt.ylabel(r'$\rho(t)$ (kg/$m^{3}$)',fontsize=16)
plt.legend(fontsize=16)
plt.tight_layout()

plt.figure()
plt.plot(t,Et(t),label=r'$E(t)$')
plt.plot(t,recoE,label=r'PWE - $N=16$')
plt.plot(t,recoE2,label=r'PWE - $N=3$')
# plt.plot(np.linspace(0,phi*tau,Ni),np.linspace(E1,E1,Ni))
# plt.plot(np.linspace(phi*tau,tau,Ni),np.linspace(E2,E2,Ni))
plt.xlabel(r'$t$ (s)',fontsize=16)
plt.ylabel(r'$E(t)$ (Pa)',fontsize=16)
plt.legend(fontsize=16)
plt.tight_layout()
#,color='gray',linestyle='--')

plt.figure()
plt.plot(t,np.sqrt(Et(t)/rhot(t)),label=r'$c(t)$')
plt.plot(t,np.sqrt(recoE/recoR),label=r'PWE - $N=16$')
plt.plot(t,np.sqrt(recoE2/recoR2),label=r'PWE - $N=3$')
# plt.plot(np.linspace(0,phi*tau,Ni),np.linspace(E1,E1,Ni))
# plt.plot(np.linspace(phi*tau,tau,Ni),np.linspace(E2,E2,Ni))
plt.xlabel(r'$t$ (s)',fontsize=16)
plt.ylabel(r'$c(t)$ (m/s)',fontsize=16)
plt.legend(fontsize=16)
plt.tight_layout()


Z1=rho1*c1
Z2=rho2*c2

tau1 = phi*tau 
tau2 = (1-phi)*tau 

def D_micro(k,omega):
    Dmicro=np.cos(omega*tau)-np.cos(c1*k*tau1)*np.cos(c2*k*tau2)+1/2*(Z1/Z2+Z2/Z1)*np.sin(c1*k*tau1)*np.sin(c2*k*tau2)
    return Dmicro

Nome=200
Nk=400
fsB=20
fsG=10
c0homo=np.sqrt(E/rho)
lambda0B=c0homo/fsB
lambda0G=c0homo/fsG
k0B=2*np.pi/lambda0B
k0G=2*np.pi/lambda0G

Om_plot = np.linspace(0,2*np.pi/tau,Nome) 
k_map = np.linspace(-25/c0homo/tau,25/c0homo/tau,Nk)


DD_micro = np.zeros([Nome,Nk])
DD_det = np.zeros([Nome,Nk])


for ind in range(Nk):
    for ind2 in range(Nome):
        DD_micro[ind2,ind] = D_micro(k_map[ind],Om_plot[ind2])
        #DD_det(ind2,ind) = Det_micro(k_map(ind),Omega_plot(ind2)






plt.figure()
#plt.pcolor(k_map,Om_plot,np.log10(np.abs(DD_micro)))
for i in range(NE):
    plt.plot((VP[i,:]),om,'red',linewidth=2)
for i in range(knth.shape[0]):
    plt.axvline(knth[i],linewidth=2)
plt.xlabel(r'$k$ (rad/m)',fontsize=16)
plt.ylabel(r'$\omega$ (rad/s)',fontsize=16)
# plt.legend(fontsize=16)
plt.title(r'$N=16$',fontsize=16)
plt.tight_layout()

plt.figure()
plt.pcolor(k_map,Om_plot,np.log10(np.abs(DD_micro)))
# for i in range(NE2):
#     plt.plot((VP2[i,:]),om,'red',linewidth=3)
plt.xlabel(r'$k$ (rad/m)',fontsize=16)
plt.ylabel(r'$\omega$ (rad/s)',fontsize=16)
# plt.legend(fontsize=16)
# plt.title(r'$N=3$',fontsize=16)
plt.tight_layout()
    
    
plt.figure()
# plt.pcolor(k_map*tau*c0homo/2/np.pi,Om_plot*tau,np.log10(np.abs(DD_micro)))
for i in range(NE):
    plt.plot((VP[i,:])*tau*c0homo/2/np.pi,om*tau,'red',linewidth=3)
plt.xlabel(r'$kc_0\tau_m$ ',fontsize=16)
plt.ylabel(r'$\omega\tau_m$ (rad)',fontsize=16)
plt.yticks([-np.pi,0,np.pi,2*np.pi],[r'$-\pi$',0,r'$\pi$',r'$2\pi$'],fontsize=16)
plt.xlim([0,2])
plt.ylim([0,2*np.pi])
# plt.legend(fontsize=16)
# plt.title(r'$N=16$',fontsize=16)
plt.tight_layout()

plt.figure()
plt.pcolor(k_map*tau*cr/2/np.pi,Om_plot*tau/2/np.pi,np.log10(np.abs(DD_micro)))
for i in range(NE2):
    plt.plot((VP2[i,:])*tau*cr/2/np.pi,om*tau/2/np.pi,'red',linewidth=3)
plt.xlabel(r'$k\,c_r\tau$ ',fontsize=16)
plt.ylabel(r'$\omega\tau$ ',fontsize=16)
# plt.legend(fontsize=16)
plt.title(r'$N=3$',fontsize=16)
plt.tight_layout()

plt.figure()
#plt.pcolor(k_map*tau*cr/2/np.pi,Om_plot*tau/2/np.pi,np.log10(np.abs(DD_micro)))
for i in range(NE):
    plt.plot((VP[i,:])*tau*cr/2/np.pi,om*tau,'red')
plt.axvline(x=k0B*cr*tau/2/np.pi,color='gray',linestyle='--')
plt.axvline(x=k0G*cr*tau/2/np.pi,color='gray',linestyle='--')

plt.figure()
#plt.pcolor(k_map*tau*cr/2/np.pi,Om_plot*tau/2/np.pi,np.log10(np.abs(DD_micro)))
for i in range(NE):
    plt.plot((VP[i,:]),om,'red')
plt.axvline(x=k0B,color='gray',linestyle='--')
plt.axvline(x=k0G,color='blue',linestyle='--')

#--------- approximations du k-gap 1 -----------
def FuncK1m(x):
    y =Omega/(2*Cel*np.sqrt(1+x/2))
    return y

def FuncK1p(x):
    y =Omega/(2*Cel*np.sqrt(1-x/2))
    return y

#--------- approximations du k-gap 2 -----------
def FuncK2m(x):
    alpha=6/x/x*(1-np.sqrt(1-1/3*x**2))
    y =Omega/Cel*np.sqrt(alpha)
#   y=Omega/Cel*(1+x*x/24)
    return y

def FuncK2p(x):
    Delta = 6 / 7 / (x*x) * (1 - np.sqrt(1 - 7 * x**2/ 3))
    y =Omega/Cel*np.sqrt(Delta)
#   y=Omega/Cel*(1+7*x*x/24)
    return y

def FuncK2mb(x):
    alpha=6/x/x*(1-np.sqrt(1-1/3*x**2))
    y =Omega/Cel*np.sqrt(alpha)
    y=Omega/Cel*(1+x*x/24)
    return y

def FuncK2pb(x):
    Delta = 6 / 7 / (x*x) * (1 - np.sqrt(1 - 7 * x**2/ 3))
    y =Omega/Cel*np.sqrt(Delta)
    y=Omega/Cel*(1+7*x*x/24)
    return y

#--------- tracer des k-gaps -----------
Omega=20
Cel=2800
Eps=1e-3
Nu=np.linspace(Eps,1-Eps,100)
K1m=FuncK1m(Nu)
K1p=FuncK1p(Nu)
K2m=FuncK2m(Nu)
K2p=FuncK2p(Nu)
K2mb=FuncK2mb(Nu)
K2pb=FuncK2pb(Nu)
X=(0,max(K2p))
Y=(0.75,0.75)


plt.figure()
for i in range(NE):
    plt.plot((VP[i,:]),om,'gray',linewidth=2)
plt.axvline(max(VP[0]),color='r',linewidth=2)
plt.axvline(min(VP[1]),color='r',linewidth=2)
plt.axvline(max(VP[1]),color='r',linewidth=2)
plt.axvline(min(VP[2]),color='r',linewidth=2)
plt.xlabel(r'$k$',fontsize=16)
plt.ylabel(r"$\omega$",fontsize=16)
plt.twinx()
plt.plot(K1m, Nu, linewidth=1.5, color='k',label='k1-')
plt.plot(K1p, Nu, linewidth=1.5, color='k',label='k1+')
plt.plot(K2m, Nu, linewidth=1.5, color='k',label='k2-')
plt.plot(K2p, Nu, linewidth=1.5, color='k',label='k2+')
plt.plot(K2mb, Nu, linewidth=1.5, color='g',label='k2-')
plt.plot(K2pb, Nu, linewidth=1.5, color='g',label='k2+')
plt.plot(X,Y,'--',color='r',linewidth=1.5)
plt.ylabel(r"$\nu$",fontsize=16)
# for i in range(knth.shape[0]):
#     plt.axvline(knth[i],linewidth=2)
# plt.legend()
plt.tight_layout()
plt.show() 

