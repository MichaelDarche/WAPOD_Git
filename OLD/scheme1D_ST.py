# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import inputReading
import source
############ Geometrical configuration ############

###################################################
def timeStep(CFL,dx,cm):
    dt=CFL*dx/cm
    return dt

##### Material properties for Nmat materials #####
#% Construction of the operator A        
####################################################
def initCauchyProblem(configuration,frontiere,mat,sourceparameters):
    X,xmin,xmax,Nx,dx=inputReading.geometric(configuration)
    Nmat,rho,cm=inputReading.material(mat)
    fs,CFL,tshift=inputReading.source(sourceparameters)
    alphaF=inputReading.frontiere(frontiere)
    rhox,Celx=inputReading.affectMaterials(alphaF, rho, cm, xmin,xmax, Nx, dx)
    # dt=timeStep(CFL,dx,cm)
    # lambdaOnde=cm/fs
    # xsmin=cm*tshift-lambdaOnde
    # xsmax=cm*tshift
    # amp=float(sourceparameters["Force"])
    coefsSource=np.zeros(2)
    if sourceparameters["Onde"]=="V":
        coefsSource[0]=1/Celx[0]
        coefsSource[1]=-rhox[0]
    if sourceparameters["Onde"]=="S":
        coefsSource[0]=1/Celx[0]
        coefsSource[1]=-rhox[0]
    V0=np.zeros([2,Nx])
    for i in range(Nx):
        V0[0,i]=coefsSource[0]*source.choice_timefct(sourceparameters, tshift-X[i]/Celx[0])
        V0[1,i]=coefsSource[1]*source.choice_timefct(sourceparameters, tshift-X[i]/Celx[0])
    return V0

def initPointSourceProblem(configuration,frontiere,mat,sourceparameters):
    X,xmin,xmax,Nx,dx=inputReading.geometric(configuration)
    # Nmat,rho,cm=inputReading.material(mat)
    # fs,CFL,tshift=inputReading.source(sourceparameters)
    # dt=timeStep(CFL,dx,cm)
    # lambdaOnde=cm/fs
    # xsmin=cm*tshift-lambdaOnde
    # xsmax=cm*tshift
    V0=np.zeros([2,Nx])
    return V0
############### Simulation parameters ##############

####################################################

#####################################################        
### Settings of the integration scheme ###
def coefScheme(K):
    if K==2:
        # print("Lax-Wendroff")
        gamma=np.zeros([K+1,K])
        gamma[0,0]=-1/2
        gamma[0,1]=-1/2
        gamma[1,0]=0
        gamma[1,1]=1
        gamma[2,0]=1/2
        gamma[2,1]=-1/2
        #print(gamma)
    elif K==4:
        # print("ADER4")
        gamma=np.zeros([K+1,K])
        gamma[0,0]=1/12
        gamma[0,1]=1/24
        gamma[0,2]=-1/12
        gamma[0,3]=-1/24
        gamma[1,0]=-2/3
        gamma[1,1]=-2/3
        gamma[1,2]=1/6
        gamma[1,3]=1/6
        gamma[2,0]=0
        gamma[2,1]=5/4
        gamma[2,2]=0
        gamma[2,3]=-1/4
        gamma[3,0]=2/3
        gamma[3,1]=-2/3
        gamma[3,2]=-1/6
        gamma[3,3]=1/6
        gamma[4,0]=-1/12
        gamma[4,1]=1/24
        gamma[4,2]=1/12
        gamma[4,3]=-1/24
    elif np.floor(K/2)!=K/2:
        print("Odd number not accepted")
    elif K>4:
        gamma=np.zeros([K+1,K])
        for i in range(K+1):
            for j in range(K):
                texte='coef '+str(i)+','+str(j)+":"
                gamma[i,j]=float(input(texte))
    return gamma

#% Computation of CS
def CS(K,s,A,dt,dx):
    gamma=coefScheme(K)
    cs=np.zeros(np.shape(A))
    for m in range(K):
        print(m)
        csm=gamma[s,m]*(dt/dx)**(m+1)*np.linalg.matrix_power(A,m+1)
        cs=cs+csm
    return cs


########### Integration scheme ############
## Finite differences without source term
def FD_cauchy(U0,Nt,Nx,K,rhox,cmx,dt,dx,film):
    U=U0.copy()
    V=np.zeros([int(Nt),int(Nx)])
    Sig=np.zeros([int(Nt),int(Nx)])
    for ts in range(Nt):
        Un = np.zeros([2,Nx])
        for xstep in range(2,Nx-2):
            Us=np.zeros([2])
            Us=Us.T
            for s in range(K+1):
                rhos=rhox[xstep+s-int(K/2)]
                cms=cmx[xstep+s-int(K/2)]
                A=inputReading.comA(rhos,cms)
                cs=CS(K,s,A,dt,dx)
                sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
                Us=Us+sumUs
            Un[:,xstep]=U[:,xstep]-Us
        U=Un.copy()
        V[ts,:]=U[0,:]
        Sig[ts,:]=U[1,:]
    if film=="yes":
        return U,V,Sig
    else:
        return U

## Finite differences with source term

def dirac_source(x,x_0,dx):
    if abs(x-x_0)<dx/2:
        return 1.0/dx
    elif abs(x-x_0)==dx/2:
        return 1.0/dx/2
    else:
        return 0.0

# def sce_dirac(sce_t,x,x_0,dx):
#     def sced(t,x):
#         sce=dirac_source(x, x_0, dx)*sce_t
#         return sce
#     return sced

def FD_sourceV(U0,Nt,Nx,K,rho,cm,dt,dx,x_0,sce_t,film):
    U=U0.copy()
    Utn=np.zeros([4,int(Nx)])
    V=np.zeros([int(Nt),int(Nx)])
    Sig=np.zeros([int(Nt),int(Nx)])
    for ts in range(Nt):
        Un = np.zeros([2,Nx])
        for xstep in range(2,Nx-2):
            Us=np.zeros([2])
            Us=Us.T
            for s in range(K+1):
                A=inputReading.comA(rho[xstep+s-int(K/2)],cm[xstep+s-int(K/2)])
                cs=CS(K,s,A,dt,dx)
                sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
                Us=Us+sumUs
            Un[:,xstep]=U[:,xstep]-Us
            Utn[0:2,xstep]=Un[:,xstep]
            Utn[2:4,xstep]=Un[:,xstep]+np.array([1.,0.])*dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
        U=Utn[2:4,:].copy()
        V[ts,:]=U[0,:]
        Sig[ts,:]=U[1,:]
    if film=="yes":
        return U,V,Sig
    else:
        return U

def FD_sourceS(U0,Nt,Nx,K,rho,cm,dt,dx,x_0,sce_t,film):
    U=U0.copy()
    Utn=np.zeros([4,int(Nx)])
    V=np.zeros([int(Nt),int(Nx)])
    Sig=np.zeros([int(Nt),int(Nx)])
    for ts in range(Nt):
        Un = np.zeros([2,Nx])
        for xstep in range(2,Nx-2):
            Us=np.zeros([2])
            Us=Us.T
            for s in range(K+1):
                A=inputReading.comA(rho[xstep+s-int(K/2)],cm[xstep+s-int(K/2)])
                cs=CS(K,s,A,dt,dx)
                sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
                Us=Us+sumUs
            Un[:,xstep]=U[:,xstep]-Us
            Utn[0:2,xstep]=Un[:,xstep]
            Utn[2:4,xstep]=Un[:,xstep]+np.array([0.,1.])*dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
        U=Utn[0:2,:].copy()
        V[ts,:]=U[0,:]
        Sig[ts,:]=U[1,:]
    if film=="yes":
        return U,V,Sig
    else:
        return U
            

