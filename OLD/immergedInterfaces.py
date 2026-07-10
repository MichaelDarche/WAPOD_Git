#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 17:28:26 2023

@author: michael
"""
import inputReading
import numpy as np

# Cas parfait
def matC(rho,c,order):
    A=inputReading.comA(rho,c)
    matrice=np.zeros([2*(order+1),2*(order+1)])
    matrice[0:2,0:2]=np.identity(2)
    for m in range(1,order+1):
        matrice[2*m:2*m+2,2*m:2*m+2]=np.linalg.matrix_power(A,m)
    return matrice 

def SMMatrix(spring,mass):
    # B=np.zeros([2,2])
    B=np.array([[0,1/spring],[mass,0]])
    return B

def matC_SM(rho,c,order,B):
    A=inputReading.comA(rho,c)
    F=np.dot(B,A)
    matrice=np.zeros([2*(order+1),2*(order+1)])
    matrice[0:2,0:2]=np.identity(2)
    for m in range(1,order+1):
        matrice[2*m:2*m+2,2*m:2*m+2]=np.linalg.matrix_power(A,m)
        matrice[2*m-2:2*m,2*m:2*m+2]=np.dot(F,np.linalg.matrix_power(A,m-1))/2
    return matrice

def normalizeMatrix(Mat0,Mat1):
    maxi=0.
    Mat0N=np.zeros(Mat0.shape)
    Mat1N=np.zeros(Mat1.shape)
    maxi0=np.abs(Mat0.max())
    maxi1=np.abs(Mat1.max())
    mini0=np.abs(Mat0.min())
    mini1=np.abs(Mat1.min())
    maxi=max(maxi0,maxi1,mini0,mini1)
    # for i in range(Mat0.shape[0]):
        # for j in range(Mat0.shape[0]):
        #     if np.abs(Mat0[i,j])>maxi:
        #         maxi=np.abs(Mat0[i,j])
        #     if np.abs(Mat1[i,j])>maxi:
        #         maxi=np.abs(Mat1[i,j])
    if maxi !=0.:
        Mat0N=Mat0/maxi
        Mat1N=Mat1/maxi
    return Mat0N, Mat1N

def invertMatrix(A):
    sizeOfMatrix=np.shape(A)
    if sizeOfMatrix[0]==sizeOfMatrix[1]:
        return np.linalg.inv(A)
    else:
        return np.linalg.pinv(A)
    
def matTaylor(dim,order,alpha,x):
    dist=x-alpha
    matrixT=np.zeros([dim,dim*(order+1)])
    for i in range(order+1):
        for j in range(dim):
            coef=1/np.math.factorial(i)*dist**i
            matrixT[j,dim*i+j]=coef
    return matrixT

def matCInter(mat0,mat1):
    mat0N, Mat1N=normalizeMatrix(mat0,mat1)
    mat1Inv=invertMatrix(Mat1N)
    return np.dot(mat1Inv,mat0N)

def matMR(dim, order,alpha,x,Npt,matCL,matCR):
    M=np.zeros([2*2*Npt,2*(order+1)])
    for xi in range(Npt):
        matT=matTaylor(dim, order, alpha, x[xi])
        M[2*xi:2*xi+2,:]=matT
    for xi in range(Npt,2*Npt): 
        matT=matTaylor(dim, order, alpha, x[xi])
        MCC=matCInter(matCL,matCR)
        M[2*xi:2*xi+2,:]=np.dot(matT,MCC)
    return M

# def matML(dim, order,alpha,x,Npt,matCL,matCR):
#     M=np.zeros([2*2*Npt,2*(order+1)])
#     for xi in range(Npt):
#         matT=matTaylor(dim, order, alpha, x[xi])
#         MCC=matCInter(matCR,matCL)
#         M[2*xi:2*xi+2,:]=np.dot(matT,MCC)
#     for xi in range(Npt,2*Npt): 
#         matT=matTaylor(dim, order, alpha, x[xi])
#         M[2*xi:2*xi+2,:]=matT
#     return M

def matMm1(M):
    return invertMatrix(M)

def matR(dim,order,alpha,x,Npt,matCL,matCR):
    M=matMR(dim, order,alpha,x,Npt,matCL,matCR)
    Mm1=matMm1(M)
    return Mm1

# def matL(dim,order,alpha,x,Npt,matCL,matCR):
#     M=matML(dim, order,alpha,x,Npt,matCL,matCR)
#     Mm1=matMm1(M)
#     return Mm1


def modifiedValuesR(U,Uint,mR,j,Npt,dim,order, alpha, x):
    Ue=np.dot(mR,Uint)
    Ur=U.copy()
    for i in range(Npt,2*Npt):
        mtay=matTaylor(dim, order, alpha, x[i])
        Umod=np.dot(mtay, Ue)
        Ur[0,j-Npt+1+i]=Umod[0,0]
        Ur[1,j-Npt+1+i]=Umod[1,0]
    return Ur

def modifiedValuesL(U,Uint,mL,j,Npt,dim,order, alpha, x):
    Ue=np.dot(mL,Uint)
    Ur=U.copy()
    for i in range(Npt):
        mtay=matTaylor(dim, order, alpha, x[i])
        Umod=np.dot(mtay,Ue)
        Ur[0,j-Npt+1+i]=Umod[0,0]
        Ur[1,j-Npt+1+i]=Umod[1,0]
    return Ur


# def modifiedValuesR(U,xstep,alpha):
#     mT=matP(U,xstep,dim,order,alpha,x,Npt,rhoL,cL,rhoR,cR,dx)
#     Ur[:,j+Npt-2]=
#     Ur[:,j+Npt-1]=
#     Ur[:,j+Npt]=
#     U[:,j-Npt+1:j+Npt]=Ur
#     return U
    
    

# ESIM=3
# order=2*ESIM-1
# Npt=3
# Xp=np.linspace(-Npt+1, Npt,2*Npt)
# rhoL=rho[0]
# rhoR=rho[1]
# cL=cm[0]
# cR=cm[1]

# matCL=matC(rhoL,cL,order)
# matCR=matC(rhoR,cR,order)
# testcm0=matCInter(matCL,matCR)

# Testlor=matTaylor(2,order,Npt/2,Xp[3])

# test=matMR(2, order,0.5,Xp,Npt,matCL,matCR)
# #Mm1=invertMatrix(test)

# mT=matR(2,order,0.5,Xp,Npt,matCL,matCR)
        