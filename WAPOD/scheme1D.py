#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import inputReading
import immergedInterfaces
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
    X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(configuration)
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
    X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(configuration)
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
def CSeig(K,s,A,dt,dx):
    eigenvalues, eigenvectors = np.linalg.eig(A)
    V = eigenvectors
    Λ = np.diag(eigenvalues)
    Vm1 = np.linalg.inv(V)
    gamma=coefScheme(K)
    csm=np.zeros(np.shape(A))
    for m in range(K):
        csym=gamma[s,m]*(dt/dx)**(m+1)*np.linalg.matrix_power(Λ,m+1)
        csm=csm+csym
    cs=np.dot(np.dot(V,csm),Vm1)
    return cs

def CS(K,s,A,dt,dx):
    gamma=coefScheme(K)
    csm=np.zeros(np.shape(A))
    for m in range(K):
        csym=gamma[s,m]*(dt/dx)**(m+1)*np.linalg.matrix_power(A,m+1)
        csm=csm+csym
    #cs=np.dot(np.dot(V,csm),Vm1)
    return csm

# def MatVec(A,X):
#     """ calcul le produit Y=A.X de A par X """ 
#     n,p = A.shape
#     Y = np.zeros(n)
#     for i in range(n):
#         S = 0.
#         for j in range(p):
#             S +=  A[i, j] * X[j]
#         Y[i] = S
#     return Y

########### Integration scheme ############
## Finite differences without source term
def FD_cauchy(U0,Nt,Nx,K,rhox,cmx,dt,dx,film):
    U=U0.copy()
    V=np.zeros([int(Nt),int(Nx)])
    Sig=np.zeros([int(Nt),int(Nx)])
    csTot=constructCS(K,Nx,rhox,cmx,dt,dx)
    for ts in range(Nt):
        # print(ts)
        Un = np.zeros([2,Nx])
        for xstep in range(int(K/2),Nx-int(K/2)):
            Us=np.zeros([2,1])
            # rhos=rhox[xstep]
            # cms=cmx[xstep]
            # A=inputReading.comA(rhos,cms)
            for s in range(K+1):
                cs=csTot[s,xstep]
                # cs=CS(K,s,A,dt,dx)
                sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
                Us[0,0]=Us[0,0]+sumUs[0]
                Us[1,0]=Us[1,0]+sumUs[1]
            Un[0,xstep]=U[0,xstep]-Us[0,0]
            Un[1,xstep]=U[1,xstep]-Us[1,0]
        U=Un.copy()
        V[ts,:]=U[0,:]
        Sig[ts,:]=U[1,:]
    if film=="yes":
        return U,V,Sig
    else:
        return U

def preTMat(frontiere,rho,cm,order):
    Ninter=int(frontiere["nombre de frontieres "])
    matCtot=np.empty((4, Ninter), dtype=object)
    for Ni in range(Ninter):
        contacti="Contact_"+str(Ni)
        if frontiere[contacti]=="parfait":
            matCtot[0,Ni]=immergedInterfaces.matC(rho[Ni], cm[Ni], order)
            matCtot[1,Ni]=immergedInterfaces.matC(rho[Ni+1], cm[Ni+1], order)
            matCtot[2,Ni]=immergedInterfaces.matCInter(matCtot[0,Ni], matCtot[1,Ni])
        if frontiere[contacti]=="masse-ressort":
            nameKni="Kn_"+str(Ni)
            nameMni="Mn_"+str(Ni)
            F=immergedInterfaces.SMMatrix(float(frontiere[nameKni]),float(frontiere[nameMni]))
            matCtot[0,Ni]=immergedInterfaces.matC_SM(rho[Ni], cm[Ni], order,F)
            matCtot[1,Ni]=immergedInterfaces.matC_SM(rho[Ni+1], cm[Ni+1], order,-F)
            matCtot[2,Ni]=immergedInterfaces.matCInter(matCtot[0,Ni], matCtot[1,Ni])
    return matCtot

def constructCS(K,Nx,rho,c,dt,dx):
    csTot=np.empty((K+1,Nx), dtype=object)
    for x in range(int(K/2),Nx-int(K/2)):
        for s in range(K+1):
            A=inputReading.comA(rho[x],c[x])
            csTot[s,x]=CS(K,s,A,dt,dx)
    return csTot

def FD_cauchyII(U0,Nt,Nx,K,rhox,cmx,dt,dx,config,mat,frontiere,film):
    X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(config)
    order=2*ESIM-1
    Nmat,rho,cm=inputReading.material(mat)
    alpha=inputReading.frontiere(frontiere)
    alpha[-1]=xmax+0.0001*xmax
    matCtot=preTMat(frontiere,rho,cm,order)
    print(matCtot)
    dim=2
    U=U0.copy()
    V=np.zeros([int(Nt),int(Nx)])
    Sig=np.zeros([int(Nt),int(Nx)])
    csTot=constructCS(K,Nx,rhox,cmx,dt,dx)
    for ts in range(Nt):
        # print(ts)
        interS=0
        Un = np.zeros([2,Nx])
        for xstep in range(int(K/2),Nx-int(K/2)):
            if (xstep*dx+xmin)<=alpha[interS]-K/2*dx:
                Us=np.zeros([2,1])
                # rhos=rhox[xstep]
                # cms=cmx[xstep]
                # A=inputReading.comA(rhos,cms)
                for s in range(K+1):
                    cs=csTot[s,xstep]
                    #cs=CS(K,s,A,dt,dx)
                    sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
                    Us[0,0]=Us[0,0]+sumUs[0]
                    Us[1,0]=Us[1,0]+sumUs[1]
                Un[0,xstep]=U[0,xstep]-Us[0,0]
                Un[1,xstep]=U[1,xstep]-Us[1,0]
            else:
                xj=int(alpha[interS]/dx)
                if xstep==xj-int(K/2)+1:
                    matCL=matCtot[0,interS]
                    matCR=matCtot[1,interS]
                    x=np.linspace(xj-Npt+1, xj+Npt,2*Npt)*dx
                    mR=immergedInterfaces.matR(dim,order,alpha[interS],x,Npt,matCL,matCR)
                    mL=np.dot(matCtot[2,interS],mR)
                    print(mL)
                    Uint=np.reshape(U[:,xj-Npt+1:xj+Npt+1],(dim*2*Npt,1),order='F')
                    UmR=immergedInterfaces.modifiedValuesR(U,Uint,mR,xj,Npt,dim,order,alpha[interS],x)
                    UmL=immergedInterfaces.modifiedValuesL(U,Uint,mL,xj,Npt,dim,order,alpha[interS],x)
                if K/2*dx>alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin<=alpha[interS]:
                    Us=np.zeros([2,1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs=csTot[s,xstep]
                        #cs=CS(K,s,A,dt,dx)
                        sumUs=np.dot(cs,UmR[:,xstep+s-int(K/2)])
                        Us[0,0]=Us[0,0]+sumUs[0]
                        Us[1,0]=Us[1,0]+sumUs[1]        
                    Un[0,xstep]=U[0,xstep]-Us[0,0]
                    Un[1,xstep]=U[1,xstep]-Us[1,0]
                elif xstep*dx+xmin>alpha[interS] and K/2*dx+alpha[interS]>=(xstep*dx+xmin):
                    Us=np.zeros([2,1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs=csTot[s,xstep]
                        #cs=CS(K,s,A,dt,dx)
                        sumUs=np.dot(cs,UmL[:,xstep+s-int(K/2)])
                        Us[0,0]=Us[0,0]+sumUs[0]
                        Us[1,0]=Us[1,0]+sumUs[1]
                    Un[0,xstep]=U[0,xstep]-Us[0,0]
                    Un[1,xstep]=U[1,xstep]-Us[1,0]
                else:
                    Us=np.zeros([2,1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    #A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs=csTot[s,xstep]
                        #cs=CS(K,s,A,dt,dx)
                        sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
                        Us[0,0]=Us[0,0]+sumUs[0]
                        Us[1,0]=Us[1,0]+sumUs[1]
                    Un[0,xstep]=U[0,xstep]-Us[0,0]
                    Un[1,xstep]=U[1,xstep]-Us[1,0]
                    if xstep == xj+int(K/2)+1:
                        interS=interS+1
        U=Un.copy()
        V[ts,:]=U[0,:] 
        Sig[ts,:]=U[1,:]
    if film=="no":
        return U
    else:
        return U,V[::film],Sig[::film]

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
    Utn=np.zeros([2,int(Nx)])
    V=np.zeros([int(Nt),int(Nx)])
    Sig=np.zeros([int(Nt),int(Nx)])
    csTot=constructCS(K,Nx,rho,cm,dt,dx)
    for ts in range(Nt):
        # print(ts)
        Un = np.zeros([2,Nx])
        for xstep in range(int(K/2),Nx-int(K/2)):
            Us=np.zeros([2,1])
            for s in range(K+1):
                cs=csTot[s,xstep]
                sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
                Us[0,0]=Us[0,0]+sumUs[0]
                Us[1,0]=Us[1,0]+sumUs[1]
            Un[0,xstep]=U[0,xstep]-Us[0,0]
            Un[1,xstep]=U[1,xstep]-Us[1,0]
            Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
            Utn[1,xstep]=Un[1,xstep]
        U=Utn[0:2,:].copy()
        V[ts,:]=U[0,:]
        Sig[ts,:]=U[1,:]
    if film=="no":
        return U
    else:
        return U,V[::film],Sig[::film]

def FD_sourceS(U0,Nt,Nx,K,rho,cm,dt,dx,x_0,sce_t,film):
    U=U0.copy()
    Utn=np.zeros([2,int(Nx)])
    V=np.zeros([int(Nt),int(Nx)])
    Sig=np.zeros([int(Nt),int(Nx)])
    csTot=constructCS(K,Nx,rho,cm,dt,dx)
    for ts in range(Nt):
        # print(ts)
        Un = np.zeros([2,Nx])
        for xstep in range(int(K/2),Nx-int(K/2)):
            Us=np.zeros([2,1])
            for s in range(K+1):
                cs=csTot[s,xstep]
                sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
                Us[0,0]=Us[0,0]+sumUs[0]
                Us[1,0]=Us[1,0]+sumUs[1]
            Un[0,xstep]=U[0,xstep]-Us[0,0]
            Un[1,xstep]=U[1,xstep]-Us[1,0]
            Utn[0,xstep]=Un[0,xstep]
            Utn[1,xstep]=Un[1,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
        U=Utn[0:2,:].copy()
        V[ts,:]=U[0,:]
        Sig[ts,:]=U[1,:]
    if film=="no":
        return U
    else:
        return U,V[::film],Sig[::film]
                
def FD_sourceVII(U0,Nt,Nx,K,rhox,cmx,dt,dx,x_0,sce_t,config,mat,frontiere,film):
    X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(config)
    order=2*ESIM-1
    Nmat,rho,cm=inputReading.material(mat)
    alpha=inputReading.frontiere(frontiere)
    alpha[-1]=xmax+0.0001*xmax
    matCtot=preTMat(frontiere,rho,cm,order)
    # print(matCtot)
    dim=2
    U=U0.copy()
    Utn=np.zeros([2,int(Nx)])
    V=np.zeros([int(Nt),int(Nx)])
    Sig=np.zeros([int(Nt),int(Nx)])
    csTot=constructCS(K,Nx,rhox,cmx,dt,dx)
    for ts in range(Nt):
        print(ts)
        interS=0
        Un = np.zeros([2,Nx])
        for xstep in range(int(K/2),Nx-int(K/2)):
            if (xstep*dx+xmin)<=alpha[interS]-K/2*dx:
                Us=np.zeros([2,1])
                # rhos=rhox[xstep]
                # cms=cmx[xstep]
                #A=inputReading.comA(rhos,cms)
                for s in range(K+1):
                    cs=csTot[s,xstep]
                    #cs=CS(K,s,A,dt,dx)
                    sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
                    Us[0,0]=Us[0,0]+sumUs[0]
                    Us[1,0]=Us[1,0]+sumUs[1]
                Un[0,xstep]=U[0,xstep]-Us[0,0]
                Un[1,xstep]=U[1,xstep]-Us[1,0]
                Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                Utn[1,xstep]=Un[1,xstep]
            else:
                xj=int(alpha[interS]/dx)
                if xstep==xj-int(K/2)+1:
                    matCL=matCtot[0,interS]
                    matCR=matCtot[1,interS]
                    x=np.linspace(xj-Npt+1, xj+Npt,2*Npt)*dx
                    mR=immergedInterfaces.matR(dim,order,alpha[interS],x,Npt,matCL,matCR)
                    mL=np.dot(matCtot[2,interS],mR)
                    # print(mL)
                    Uint=np.reshape(U[:,xj-Npt+1:xj+Npt+1],(dim*2*Npt,1),order='F')
                    UmR=immergedInterfaces.modifiedValuesR(U,Uint,mR,xj,Npt,dim,order,alpha[interS],x)
                    UmL=immergedInterfaces.modifiedValuesL(U,Uint,mL,xj,Npt,dim,order,alpha[interS],x)
                if K/2*dx>alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin<=alpha[interS]:
                    Us=np.zeros([2,1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    #A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs=csTot[s,xstep]
                        #cs=CS(K,s,A,dt,dx)
                        sumUs=np.dot(cs,UmR[:,xstep+s-int(K/2)])
                        Us[0,0]=Us[0,0]+sumUs[0]
                        Us[1,0]=Us[1,0]+sumUs[1]        
                    Un[0,xstep]=U[0,xstep]-Us[0,0]
                    Un[1,xstep]=U[1,xstep]-Us[1,0]
                    Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1,xstep]=Un[1,xstep]
                elif xstep*dx+xmin>alpha[interS] and K/2*dx+alpha[interS]>=(xstep*dx+xmin):
                    Us=np.zeros([2,1])
                    #rhos=rhox[xstep]
                    #cms=cmx[xstep]
                    #A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs=csTot[s,xstep]
                        #cs=CS(K,s,A,dt,dx)
                        sumUs=np.dot(cs,UmL[:,xstep+s-int(K/2)])
                        Us[0,0]=Us[0,0]+sumUs[0]
                        Us[1,0]=Us[1,0]+sumUs[1]
                    Un[0,xstep]=U[0,xstep]-Us[0,0]
                    Un[1,xstep]=U[1,xstep]-Us[1,0]
                    Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1,xstep]=Un[1,xstep]
                else:
                    Us=np.zeros([2,1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    #A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs=csTot[s,xstep]
                        #cs=CS(K,s,A,dt,dx)
                        sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
                        Us[0,0]=Us[0,0]+sumUs[0]
                        Us[1,0]=Us[1,0]+sumUs[1]
                    Un[0,xstep]=U[0,xstep]-Us[0,0]
                    Un[1,xstep]=U[1,xstep]-Us[1,0]
                    Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1,xstep]=Un[1,xstep]
                    if xstep == xj+int(K/2)+1:
                        interS=interS+1
        U=Utn[0:2,:].copy()
        V[ts,:]=U[0,:]
        Sig[ts,:]=U[1,:]
    if film=="no":
        return U
    else:
        return U,V[::film],Sig[::film]
    
def FD_sourceVH(U0,Nt,Nx,K,rhox,cmx,dt,dx,x_0,sce_t,config,mat,frontiere,film):
    X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(config)
    order=2*ESIM-1
    Nmat,rho,cm=inputReading.material(mat)
    alpha=inputReading.frontiere(frontiere)
    alpha[-1]=xmax+0.0001*xmax
    matCtot=preTMat(frontiere,rho,cm,order)
    dim=2
    U=U0.copy()
    Utn=np.zeros([2,int(Nx)])
    V=np.zeros([int(Nt),int(Nx)])
    Sig=np.zeros([int(Nt),int(Nx)])
    csTot=constructCS(K,Nx,rhox,cmx,dt,dx)
    for ts in range(Nt):
        # print(ts)
        interS=0
        Un = np.zeros([2,Nx])
        for xstep in range(int(K/2),Nx-int(K/2)):
            if (xstep*dx+xmin)<=alpha[interS]-K/2*dx:
                Us=np.zeros([2,1])
                # rhos=rhox[xstep]
                # cms=cmx[xstep]
                #A=inputReading.comA(rhos,cms)
                for s in range(K+1):
                    cs=csTot[s,xstep]
                    #cs=CS(K,s,A,dt,dx)
                    sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
                    Us[0,0]=Us[0,0]+sumUs[0]
                    Us[1,0]=Us[1,0]+sumUs[1]
                Un[0,xstep]=U[0,xstep]-Us[0,0]
                Un[1,xstep]=U[1,xstep]-Us[1,0]
                Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                Utn[1,xstep]=Un[1,xstep]
            else:
                xj=int(alpha[interS]/dx)
                if xstep==xj-int(K/2)+1:
                    matCL=matCtot[0,interS]
                    matCR=matCtot[1,interS]
                    x=np.linspace(xj-Npt+1, xj+Npt,2*Npt)*dx
                    mR=immergedInterfaces.matR(dim,order,alpha[interS],x,Npt,matCL,matCR)
                    mL=np.dot(matCtot[2,interS],mR)
                    Uint=np.reshape(U[:,xj-Npt+1:xj+Npt+1],(dim*2*Npt,1),order='F')
                    UmR=immergedInterfaces.modifiedValuesR(U,Uint,mR,xj,Npt,dim,order,alpha[interS],x)
                    UmL=immergedInterfaces.modifiedValuesL(U,Uint,mL,xj,Npt,dim,order,alpha[interS],x)
                if K/2*dx>alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin<=alpha[interS]:
                    Us=np.zeros([2,1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    #A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs=csTot[s,xstep]
                        #cs=CS(K,s,A,dt,dx)
                        sumUs=np.dot(cs,UmR[:,xstep+s-int(K/2)])
                        Us[0,0]=Us[0,0]+sumUs[0]
                        Us[1,0]=Us[1,0]+sumUs[1]        
                    Un[0,xstep]=U[0,xstep]-Us[0,0]
                    Un[1,xstep]=U[1,xstep]-Us[1,0]
                    Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1,xstep]=Un[1,xstep]
                elif xstep*dx+xmin>alpha[interS] and K/2*dx+alpha[interS]>=(xstep*dx+xmin):
                    Us=np.zeros([2,1])
                    #rhos=rhox[xstep]
                    #cms=cmx[xstep]
                    #A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs=csTot[s,xstep]
                        #cs=CS(K,s,A,dt,dx)
                        sumUs=np.dot(cs,UmL[:,xstep+s-int(K/2)])
                        Us[0,0]=Us[0,0]+sumUs[0]
                        Us[1,0]=Us[1,0]+sumUs[1]
                    Un[0,xstep]=U[0,xstep]-Us[0,0]
                    Un[1,xstep]=U[1,xstep]-Us[1,0]
                    Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1,xstep]=Un[1,xstep]
                    
                else:
                    Us=np.zeros([2,1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    #A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs=csTot[s,xstep]
                        #cs=CS(K,s,A,dt,dx)
                        sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
                        Us[0,0]=Us[0,0]+sumUs[0]
                        Us[1,0]=Us[1,0]+sumUs[1]
                    Un[0,xstep]=U[0,xstep]-Us[0,0]
                    Un[1,xstep]=U[1,xstep]-Us[1,0]
                    Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1,xstep]=Un[1,xstep]
                    if xstep == xj+int(K/2)+1:
                        interS=interS+1
        U=Utn[0:2,:].copy()
        V[ts,:]=U[0,:]
        Sig[ts,:]=U[1,:]
    if film=="no":
        return U
    else:
        return U,V[::film],Sig[::film]


def FD_sourceSII(U0,Nt,Nx,K,rhox,cmx,dt,dx,x_0,sce_t,config,mat,frontiere,film):
    X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(config)
    order=2*ESIM-1
    Nmat,rho,cm=inputReading.material(mat)
    alpha=inputReading.frontiere(frontiere)
    alpha[-1]=xmax+0.0001*xmax
    matCtot=preTMat(frontiere,rho,cm,order)
    dim=2
    U=U0.copy()
    Utn=np.zeros([2,int(Nx)])
    V=np.zeros([int(Nt),int(Nx)])
    Sig=np.zeros([int(Nt),int(Nx)])
    csTot=constructCS(K,Nx,rhox,cmx,dt,dx)
    for ts in range(Nt):
        # print(ts)
        interS=0
        Un = np.zeros([2,Nx])
        for xstep in range(int(K/2),Nx-int(K/2)):
            if (xstep*dx+xmin)<=alpha[interS]-K/2*dx:
                Us=np.zeros([2,1])
                # rhos=rhox[xstep]
                # cms=cmx[xstep]
                #A=inputReading.comA(rhos,cms)
                for s in range(K+1):
                    cs=csTot[s,xstep]
                    #cs=CS(K,s,A,dt,dx)
                    sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
                    Us[0,0]=Us[0,0]+sumUs[0]
                    Us[1,0]=Us[1,0]+sumUs[1]
                Un[0,xstep]=U[0,xstep]-Us[0,0]
                Un[1,xstep]=U[1,xstep]-Us[1,0]
                Utn[0,xstep]=Un[0,xstep]
                Utn[1,xstep]=Un[1,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
            else:
                xj=int(alpha[interS]/dx)
                if xstep==xj-int(K/2)+1:
                    matCL=matCtot[0,interS]
                    matCR=matCtot[1,interS]
                    x=np.linspace(xj-Npt+1, xj+Npt,2*Npt)*dx
                    mR=immergedInterfaces.matR(dim,order,alpha[interS],x,Npt,matCL,matCR)
                    mL=np.dot(matCtot[2,interS],mR)
                    Uint=np.reshape(U[:,xj-Npt+1:xj+Npt+1],(dim*2*Npt,1),order='F')
                    UmR=immergedInterfaces.modifiedValuesR(U,Uint,mR,xj,Npt,dim,order,alpha[interS],x)
                    UmL=immergedInterfaces.modifiedValuesL(U,Uint,mL,xj,Npt,dim,order,alpha[interS],x)
                if K/2*dx>alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin<=alpha[interS]:
                    Us=np.zeros([2,1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    #A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs=csTot[s,xstep]
                        #cs=CS(K,s,A,dt,dx)
                        sumUs=np.dot(cs,UmR[:,xstep+s-int(K/2)])
                        Us[0,0]=Us[0,0]+sumUs[0]
                        Us[1,0]=Us[1,0]+sumUs[1]        
                    Un[0,xstep]=U[0,xstep]-Us[0,0]
                    Un[1,xstep]=U[1,xstep]-Us[1,0]
                    Utn[0,xstep]=Un[0,xstep]
                    Utn[1,xstep]=Un[1,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                elif xstep*dx+xmin>alpha[interS] and K/2*dx+alpha[interS]>=(xstep*dx+xmin):
                    Us=np.zeros([2,1])
                    #rhos=rhox[xstep]
                    #cms=cmx[xstep]
                    #A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs=csTot[s,xstep]
                        #cs=CS(K,s,A,dt,dx)
                        sumUs=np.dot(cs,UmL[:,xstep+s-int(K/2)])
                        Us[0,0]=Us[0,0]+sumUs[0]
                        Us[1,0]=Us[1,0]+sumUs[1]
                    Un[0,xstep]=U[0,xstep]-Us[0,0]
                    Un[1,xstep]=U[1,xstep]-Us[1,0]
                    Utn[0,xstep]=Un[0,xstep]
                    Utn[1,xstep]=Un[1,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    
                else:
                    Us=np.zeros([2,1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    #A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs=csTot[s,xstep]
                        #cs=CS(K,s,A,dt,dx)
                        sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
                        Us[0,0]=Us[0,0]+sumUs[0]
                        Us[1,0]=Us[1,0]+sumUs[1]
                    Un[0,xstep]=U[0,xstep]-Us[0,0]
                    Un[1,xstep]=U[1,xstep]-Us[1,0]
                    Utn[0,xstep]=Un[0,xstep]
                    Utn[1,xstep]=Un[1,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    if xstep == xj+int(K/2)+1:
                        interS=interS+1
        U=Utn[0:2,:].copy()
        V[ts,:]=U[0,:]
        Sig[ts,:]=U[1,:]
    if film=="no":
        return U
    else:
        return U,V[::film],Sig[::film]