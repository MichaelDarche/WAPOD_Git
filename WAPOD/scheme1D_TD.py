#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import sympy as sp
import inputReading
import immergedInterfaces
import source
import time
# import sumUs_c

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

def preTMat(frontiere,rho,cm,order,modulation):
    Ninter=int(frontiere["nombre de frontieres "])
    matCtot=np.empty((4, Ninter), dtype=object)
    for Ni in range(Ninter):
        nalpha="Alpha_"+str(Ni)
        alpha=float(frontiere[nalpha])
        contacti="Contact_"+str(Ni)
        if frontiere[contacti]=="parfait":
            matCtot[0,Ni]=immergedInterfaces.matC(rho[Ni], cm[Ni], order)
            matCtot[1,Ni]=immergedInterfaces.matC(rho[Ni+1], cm[Ni+1], order)
            matCtot[2,Ni]=immergedInterfaces.matCInter(matCtot[0,Ni], matCtot[1,Ni])
        if frontiere[contacti]=="masse-ressort":
            nameKni="Kn_"+str(Ni)
            nameMni="Mn_"+str(Ni)
            F=immergedInterfaces.SMMatrix(float(frontiere[nameKni]),float(frontiere[nameMni]))
            matCtot[0,Ni]=immergedInterfaces.matC_SM(rho[Ni], cm[Ni], order,-F)
            matCtot[1,Ni]=immergedInterfaces.matC_SM(rho[Ni+1], cm[Ni+1], order,F)
            matCtot[2,Ni]=immergedInterfaces.matCInter(matCtot[0,Ni], matCtot[1,Ni])
        if frontiere[contacti]=="masse-ressort module":
            if modulation["Synchro"]=="No":
            else:
                fctModul,t,x,Omega,phi=immergedInterfaces.fctMod(modulation)
                ampF=sp.Float(modulation["DeltaF"])
                ampM=sp.Float(modulation["DeltaM"])
                FreqM=sp.Float(modulation["FreqM"])
                celM=sp.Float(modulation["ModC"])
            if celM==0:
                phiM=0
            else:
                phiM=FreqM*2*np.pi*alpha/celM
            DF,DM= sp.symbols('DF DM')
            nameKni="Kn_"+str(Ni)
            nameMni="Mn_"+str(Ni)
            Fmod=(1+DF*fctModul)/sp.Float(frontiere[nameKni])
            Mmod=(1+DM*fctModul)*sp.Float(frontiere[nameMni])
            formalB=immergedInterfaces.CMMatrix(Fmod,Mmod)
            matCtot[0,Ni]=immergedInterfaces.formatC_SM_modulated(rho[Ni], cm[Ni], order,formalB,x,t,Omega,phi,DF,DM)
            matCtot[1,Ni]=immergedInterfaces.formatC_SM_modulated(rho[Ni+1], cm[Ni+1], order,-formalB,x,t,Omega,phi,DF,DM)
            matCtot[0,Ni]=matCtot[0,Ni].subs(Omega,FreqM*2*np.pi)
            matCtot[0,Ni]=matCtot[0,Ni].subs(DF,ampF)
            matCtot[0,Ni]=matCtot[0,Ni].subs(DM,ampM)
            matCtot[0,Ni]=matCtot[0,Ni].subs(phi,phiM)
            matCtot[1,Ni]=matCtot[1,Ni].subs(Omega,FreqM*2*np.pi)
            matCtot[1,Ni]=matCtot[1,Ni].subs(DF,ampF)
            matCtot[1,Ni]=matCtot[1,Ni].subs(DM,ampM)
            matCtot[1,Ni]=matCtot[1,Ni].subs(phi,phiM)
            #matCtot[2,Ni]=immergedInterfaces.matCInter(matCtot[0,Ni], matCtot[1,Ni])
    return matCtot

def constructCS(K,Nx,rho,c,dt,dx):
    csTot=np.empty((K+1,Nx),dtype=object)
    for x in range(int(K/2),Nx-int(K/2)):
        for s in range(K+1):
            A=inputReading.comA(rho[x],c[x])
            csTot[s,x]=CS(K,s,A,dt,dx)       
    return csTot

def FD_cauchyII(U0,Nt,Nx,K,rhox,cmx,dt,dx,config,mat,frontiere,modulation,film):
    X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(config)
    order=2*ESIM-1
    Nmat,rho,cm=inputReading.material(mat)
    Ninter=int(frontiere["nombre de frontieres "])
    alpha=inputReading.frontiere(frontiere)
    alpha[-1]=xmax+0.0001*xmax
    dim=2
    U=U0.copy()
    V=np.zeros([int(Nt),int(Nx)])
    Sig=np.zeros([int(Nt),int(Nx)])
    csTot=constructCS(K,Nx,rhox,cmx,dt,dx)
    formalCMat=preTMat(frontiere,rho,cm,order,modulation)
    t= sp.symbols('t')
    FMatrix=np.empty((4, Ninter), dtype=object)
    for Ni in range(Ninter):
        contacti="Contact_"+str(Ni)
        if frontiere[contacti]=="parfait" or frontiere[contacti]=="masse-ressort":
            FMatrix[:,Ni]=formalCMat[:,Ni]
        else:
            FMatrix[0,Ni]=sp.lambdify(t, formalCMat[0,Ni],modules="numpy")
            FMatrix[1,Ni]=sp.lambdify(t, formalCMat[1,Ni],modules="numpy")
    matCtot=np.empty((4, Ninter), dtype=object)
    for ts in range(Nt):
        print(ts)
        interS=0
        # Un = np.zeros([2,Nx])
        for Ni in range(Ninter):
            contacti="Contact_"+str(Ni)
            if frontiere[contacti]=="parfait" or frontiere[contacti]=="masse-ressort":
                matCtot[:,Ni]=formalCMat[:,Ni]
            else:
                matCtot[0,Ni]=immergedInterfaces.valueModMatrix(FMatrix[0,Ni],t,ts*dt)
                matCtot[1,Ni]=immergedInterfaces.valueModMatrix(FMatrix[1,Ni],t,ts*dt)
                matCtot[2,Ni]=immergedInterfaces.matCInter(matCtot[0,Ni], matCtot[1,Ni])
        # print(matCtot)
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
                
def FD_sourceVII(U0,Nt,Nx,K,rhox,cmx,dt,dx,x_0,sce_t,config,mat,frontiere,modulation,film):
    X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(config)
    order=2*ESIM-1
    Nmat,rho,cm=inputReading.material(mat)
    Ninter=int(frontiere["nombre de frontieres "])
    alpha=inputReading.frontiere(frontiere)
    alpha[-1]=xmax+0.0001*xmax
    dim=2
    U=U0.copy()
    Utn=np.zeros([2,int(Nx)])
    V=np.zeros([int(Nt),int(Nx)])
    Sig=np.zeros([int(Nt),int(Nx)])
    csTot=constructCS(K,Nx,rhox,cmx,dt,dx)
    formalCMat=preTMat(frontiere,rho,cm,order,modulation)
    t= sp.symbols('t')
    FMatrix=np.empty((4, Ninter), dtype=object)
    for Ni in range(Ninter):
        contacti="Contact_"+str(Ni)
        if frontiere[contacti]=="parfait" or frontiere[contacti]=="masse-ressort":
            FMatrix[:,Ni]=formalCMat[:,Ni]
        else:
            FMatrix[0,Ni]=sp.lambdify(t, formalCMat[0,Ni],modules="numpy")
            FMatrix[1,Ni]=sp.lambdify(t, formalCMat[1,Ni],modules="numpy")
    matCtot=np.empty((4, Ninter), dtype=object)
    #formalCMat=immergedInterfaces.formatC_SM_modulated(rho,cm,order,formalB,x,t,Lambda,Omega,DF,DM)#,alpha,ts,LambdaM,OmegaM,ampF,ampM)
    for ts in range(Nt):
        print(ts)
        interS=0
        # Un = np.zeros([2,Nx])
        for Ni in range(Ninter):
            contacti="Contact_"+str(Ni)
            if frontiere[contacti]=="parfait" or frontiere[contacti]=="masse-ressort":
                matCtot[:,Ni]=formalCMat[:,Ni]
            else:
                matCtot[0,Ni]=immergedInterfaces.valueModMatrix(FMatrix[0,Ni],t,ts*dt)
                matCtot[1,Ni]=immergedInterfaces.valueModMatrix(FMatrix[1,Ni],t,ts*dt)
                matCtot[2,Ni]=immergedInterfaces.matCInter(matCtot[0,Ni], matCtot[1,Ni])
        # print(matCtot)
        Un = np.zeros([2,Nx])
        for xstep in range(int(K/2),Nx-int(K/2)):
            if (xstep*dx+xmin)<=alpha[interS]-K/2*dx:
                # Us=np.zeros([2,1])
                # # rhos=rhox[xstep]
                # # cms=cmx[xstep]
                # #A=inputReading.comA(rhos,cms)
                Uc = U.astype(np.double)
                Uc = U.astype(np.double)
                Unc = Un.astype(np.double)
                Us=sumUs_c.compute_sumUs(Uc, Uc, Unc, csTot, xstep, K)
                # for s in range(K+1):
                #     cs=csTot[s,xstep]
                #     #cs=CS(K,s,A,dt,dx)
                #     sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
                #     Us[0,0]=Us[0,0]+sumUs[0]
                #     Us[1,0]=Us[1,0]+sumUs[1]
                # Un[0,xstep]=U[0,xstep]-Us[0,0]
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
                    # Us=np.zeros([2,1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    #A=inputReading.comA(rhos,cms)
                    Uc = U.astype(np.double)
                    UmRc = UmR.astype(np.double)
                    Unc = Un.astype(np.double)
                    csTotc = csTot.astype(np.double)
                    Us=sumUs_c.compute_sumUs(Uc, UmRc, Unc, csTotc, xstep, K)
                    # for s in range(K+1):
                    #     cs=csTot[s,xstep]
                    #     #cs=CS(K,s,A,dt,dx)
                    #     sumUs=np.dot(cs,UmR[:,xstep+s-int(K/2)])
                    #     Us[0,0]=Us[0,0]+sumUs[0]
                    #     Us[1,0]=Us[1,0]+sumUs[1]        
                    Un[0,xstep]=U[0,xstep]-Us[0,0]
                    Un[1,xstep]=U[1,xstep]-Us[1,0]
                    Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1,xstep]=Un[1,xstep]
                elif xstep*dx+xmin>alpha[interS] and K/2*dx+alpha[interS]>=(xstep*dx+xmin):
                    # Us=np.zeros([2,1])
                    #rhos=rhox[xstep]
                    #cms=cmx[xstep]
                    #A=inputReading.comA(rhos,cms)
                    Uc = U.astype(np.double)
                    UmLc = UmL.astype(np.double)
                    Unc = Un.astype(np.double)
                    csTotc = csTot.astype(np.double)
                    Us=sumUs_c.compute_sumUs(Uc, UmLc, Unc, csTotc, xstep, K)
                    # for s in range(K+1):
                    #     cs=csTot[s,xstep]
                    #     #cs=CS(K,s,A,dt,dx)
                    #     sumUs=np.dot(cs,UmL[:,xstep+s-int(K/2)])
                    #     Us[0,0]=Us[0,0]+sumUs[0]
                    #     Us[1,0]=Us[1,0]+sumUs[1]
                    # Un[0,xstep]=U[0,xstep]-Us[0,0]
                    # Un[1,xstep]=U[1,xstep]-Us[1,0]
                    Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1,xstep]=Un[1,xstep]
                else:
                    # Us=np.zeros([2,1])
                    # # rhos=rhox[xstep]
                    # # cms=cmx[xstep]
                    # #A=inputReading.comA(rhos,cms)
                    Uc = U.astype(np.double)
                    Unc = Un.astype(np.double)
                    csTotc = csTotc = np.array(csTot, dtype=np.double)#csTot.astype(np.double)
                    Us=sumUs_c.compute_sumUs(Uc, Uc, Unc, csTotc, xstep, K)
                    # for s in range(K+1):
                    #     cs=csTot[s,xstep]
                    #     #cs=CS(K,s,A,dt,dx)
                    #     sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
                    #     Us[0,0]=Us[0,0]+sumUs[0]
                    #     Us[1,0]=Us[1,0]+sumUs[1]
                    # Un[0,xstep]=U[0,xstep]-Us[0,0]
                    # Un[1,xstep]=U[1,xstep]-Us[1,0]
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

def FD_sourceVII_opti(U0,Nt,Nx,K,rhox,cmx,dt,dx,x_0,sce_t,config,mat,frontiere,modulation,film):
    X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(config)
    order=2*ESIM-1
    Nmat,rho,cm=inputReading.material(mat)
    Ninter=int(frontiere["nombre de frontieres "])
    alpha=inputReading.frontiere(frontiere)
    alpha[-1]=xmax+0.0001*xmax
    dim=2
    U=U0.copy()
    Utn=np.zeros([2,int(Nx)])
    V=np.zeros([int(Nt),int(Nx)])
    Sig=np.zeros([int(Nt),int(Nx)])
    csTot=constructCS(K,Nx,rhox,cmx,dt,dx)
    formalCMat=preTMat(frontiere,rho,cm,order,modulation,ts=0,dt=1)
    t= sp.symbols('t')
    FMatrix=np.empty((4, Ninter), dtype=object)
    for Ni in range(Ninter):
        contacti="Contact_"+str(Ni)
        if frontiere[contacti]=="parfait" or frontiere[contacti]=="masse-ressort":
            FMatrix[:,Ni]=formalCMat[:,Ni]
        else:
            FMatrix[0,Ni]=sp.lambdify(t, formalCMat[0,Ni],modules="numpy")
            FMatrix[1,Ni]=sp.lambdify(t, formalCMat[1,Ni],modules="numpy")
    matCtot=np.empty((4, Ninter), dtype=object)
    #formalCMat=immergedInterfaces.formatC_SM_modulated(rho,cm,order,formalB,x,t,Lambda,Omega,DF,DM)#,alpha,ts,LambdaM,OmegaM,ampF,ampM)
    for ts in range(Nt):
        print(ts)
        interS=0
        # Un = np.zeros([2,Nx])
        for Ni in range(Ninter):
            contacti="Contact_"+str(Ni)
            if frontiere[contacti]=="parfait" or frontiere[contacti]=="masse-ressort":
                matCtot[:,Ni]=formalCMat[:,Ni]
            else:
                matCtot[0,Ni]=immergedInterfaces.valueModMatrix(FMatrix[0,Ni],t,ts*dt)
                matCtot[1,Ni]=immergedInterfaces.valueModMatrix(FMatrix[1,Ni],t,ts*dt)
                matCtot[2,Ni]=immergedInterfaces.matCInter(matCtot[0,Ni], matCtot[1,Ni])
        for xstep in range(int(K/2),Nx-int(K/2)):
            if (xstep*dx+xmin)<=alpha[interS]-K/2*dx:
                # time0=time.time()
                Us=np.zeros([2,1])
                # rhos=rhox[xstep]
                # cms=cmx[xstep]
                #A=inputReading.comA(rhos,cms)
                time0=time.time()
                s_values = np.arange(K + 1)
                cst=np.stack(csTot[s_values,xstep])
                cst=cst.reshape(K+1, 2, 2)
                Ured=U[:,xstep-int(K/2)+s_values]
                print("Forme de cst:", cst.shape)
                print("Forme de Ured:", Ured.shape)
                # Ured=Ured.reshape(-1, 2)
                # if cst.shape[2] != Ured.shape[0]:
                #     raise ValueError("Les dimensions ne correspondent pas.")
                Us_sum =np.matmul(cst.transpose(1, 2, 0), Ured)#.sum(axis=0)#np.tensordot(cst, Ured, axes=([0], [0]))# np.einsum('ijk,ki->ij', cst, Ured).sum(axis=0)
                print(time.time()-time0)
                #np.einsum('ijk,ik->ij', cst, Ured).sum(axis=0, keepdims=True).T
                # #Us = np.dot(Ured,cst)
                time0=time.time()
                Us=np.zeros([2,1])
                for s in range(K+1):
                    cs=csTot[s,xstep]
                    #cs=CS(K,s,A,dt,dx)
                    sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
                    #print(sumUs)
                    Us[0,0]=Us[0,0]+sumUs[0]
                    Us[1,0]=Us[1,0]+sumUs[1]
                # print(Us)
                print(time.time()-time0,'time')
                # Un[0,xstep]=U[0,xstep]-Us[0,0]
                # Un[1,xstep]=U[1,xstep]-Us[1,0]
                Utn[0,xstep]=U[0,xstep]-Us_sum[0]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                Utn[1,xstep]=U[1,xstep]-Us_sum[1]
                # timef=time.time()-time0
                # print(timef)
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
                    s_values = np.arange(K + 1)
                    cst=np.stack(csTot[s_values,xstep])
                    cst=cst.reshape(K+1, 2, 2)
                    Ured=UmR[:,xstep-int(K/2)+s_values]
                    # Ured=Ured.reshape(-1, 2)
                    # if cst.shape[2] != Ured.shape[0]:
                    #     raise ValueError("Les dimensions ne correspondent pas.")
                    Us = np.einsum('ijk,ki->ij', cst, Ured)#np.einsum('ijk,ik->ij', cst, Ured).sum(axis=0, keepdims=True).T
                    Us_sum = np.sum(Us, axis=0)
                    # for s in range(K+1):
                    #     cs=csTot[s,xstep]
                    #     #cs=CS(K,s,A,dt,dx)
                    #     sumUs=np.dot(cs,UmR[:,xstep+s-int(K/2)])
                    #     Us[0,0]=Us[0,0]+sumUs[0]
                    #     Us[1,0]=Us[1,0]+sumUs[1]        
                    # Un[0,xstep]=U[0,xstep]-Us[0,0]
                    # Un[1,xstep]=U[1,xstep]-Us[1,0]
                    Utn[0,xstep]=U[0,xstep]-Us_sum[0]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1,xstep]=U[1,xstep]-Us_sum[1]
                elif xstep*dx+xmin>alpha[interS] and K/2*dx+alpha[interS]>=(xstep*dx+xmin):
                    Us=np.zeros([2,1])
                    #rhos=rhox[xstep]
                    #cms=cmx[xstep]
                    #A=inputReading.comA(rhos,cms)
                    s_values = np.arange(K + 1)
                    cst=np.stack(csTot[s_values,xstep])
                    cst=cst.reshape(K+1, 2, 2)
                    Ured=UmL[:,xstep-int(K/2)+s_values]
                    # Ured=Ured.reshape(-1, 2)
                    # if cst.shape[2] != Ured.shape[0]:
                    #     raise ValueError("Les dimensions ne correspondent pas.")
                    Us = np.einsum('ijk,ki->ij', cst, Ured)#np.einsum('ijk,ik->ij', cst, Ured).sum(axis=0, keepdims=True).T
                    Us_sum = np.sum(Us, axis=0)
                    # for s in range(K+1):
                    #     cs=csTot[s,xstep]
                    #     #cs=CS(K,s,A,dt,dx)
                    #     sumUs=np.dot(cs,UmL[:,xstep+s-int(K/2)])
                    #     Us[0,0]=Us[0,0]+sumUs[0]
                    #     Us[1,0]=Us[1,0]+sumUs[1]
                    # Un[0,xstep]=U[0,xstep]-Us[0,0]
                    # Un[1,xstep]=U[1,xstep]-Us[1,0]
                    Utn[0,xstep]=U[0,xstep]-Us_sum[0]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1,xstep]=U[1,xstep]-Us_sum[1]
                else:
                    Us=np.zeros([2,1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    #A=inputReading.comA(rhos,cms)
                    s_values = np.arange(K + 1)
                    cst=np.stack(csTot[s_values,xstep])
                    cst=cst.reshape(K+1, 2, 2)
                    Ured=U[:,xstep-int(K/2)+s_values]
                    # Ured=Ured.reshape(-1, 2)
                    Us = np.einsum('ijk,ki->ij', cst, Ured)
                    Us_sum = np.sum(Us, axis=0)
                    # for s in range(K+1):
                    #     cs=csTot[s,xstep]
                    #     #cs=CS(K,s,A,dt,dx)
                    #     sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
                    #     Us[0,0]=Us[0,0]+sumUs[0]
                    #     Us[1,0]=Us[1,0]+sumUs[1]
                    # Un[0,xstep]=U[0,xstep]-Us[0,0]
                    # Un[1,xstep]=U[1,xstep]-Us[1,0]
                    Utn[0,xstep]=U[0,xstep]-Us_sum[0]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1,xstep]=U[1,xstep]-Us_sum[1]
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


def FD_sourceSII(U0,Nt,Nx,K,rhox,cmx,dt,dx,x_0,sce_t,config,mat,frontiere,modulation,film):
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
        print(ts)
        interS=0
        Un = np.zeros([2,Nx])
        matCtot=preTMat(frontiere,rho,cm,order,modulation,ts*dt,dt)
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