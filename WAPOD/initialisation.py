#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 09:30:39 2025

@author: michael
"""

import numpy as np
import sympy as sp
import inputReading
import source

############ Geometrical configuration ############

###################################################


def timeStep(CFL, dx, cm):
    dt = CFL*dx/cm
    return dt

##### Material properties for Nmat materials #####
# % Construction of the operator A
####################################################


def initCauchyProblem(configuration, frontiere, mat, sourceparameters):
    X, xmin, xmax, Nx, dx, ESIM, Npt = inputReading.geometric(configuration)
    Nmat, rho, cm = inputReading.material(mat)
    fs, CFL, tshift = inputReading.source(sourceparameters)
    alphaF = inputReading.frontiere(frontiere)
    rhox, Celx = inputReading.affectMaterials(
        alphaF, rho, cm, xmin, xmax, Nx, dx)
    # dt=timeStep(CFL,dx,cm)
    coefsSource = np.zeros(2)
    if sourceparameters["Onde"] == "V":
        coefsSource[0] = 1/Celx[0]
        coefsSource[1] = -rhox[0]
    if sourceparameters["Onde"] == "S":
        coefsSource[0] = 1/Celx[0]
        coefsSource[1] = -rhox[0]
    V0 = np.zeros([2, Nx])
    for i in range(Nx):
        V0[0, i] = coefsSource[0] * \
            source.choice_timefct(sourceparameters, tshift-X[i]/Celx[0])
        V0[1, i] = coefsSource[1] * \
            source.choice_timefct(sourceparameters, tshift-X[i]/Celx[0])
    return V0

def initCauchyProblemH(configuration, frontiere, mat, sourceparameters,modulation):
    import homogenize
    X, xmin, xmax, Nx, dx, ESIM, Npt = inputReading.geometric(configuration)
    Nmat, rho, cm = inputReading.material(mat)
    fs, CFL, tshift = inputReading.source(sourceparameters)
    alphaF = inputReading.frontiere(frontiere)
    rhox, Celx = inputReading.affectMaterials(
        alphaF, rho, cm, xmin, xmax, Nx, dx)
    if frontiere["Contact_"+str(0)] == "masse-ressort":
        rhoeff, Ceff = homogenize.propHomo(
            0, rhox[0], Celx[0], frontiere, modulation)
    if frontiere["Contact_"+str(0)] == "masse-ressort module":
        rhoeff, Eeff, Ceff = homogenize.propHomoForm(
            rhox[0], Celx[0], frontiere, modulation)
        t = sp.symbols("t")
    reff=sp.lambdify(t, rhoeff, modules='numpy')
    Ceeff=sp.lambdify(t, Ceff, modules='numpy')
    # dt=timeStep(CFL,dx,cm)
    coefsSource = np.zeros(2)
    if sourceparameters["Onde"] == "V":
        coefsSource[0] = 1/Ceeff(0)
        coefsSource[1] = -reff(0)
    if sourceparameters["Onde"] == "S":
        coefsSource[0] = 1/Ceeff(0)
        coefsSource[1] = -reff(0)
    V0 = np.zeros([2, Nx])
    for i in range(Nx):
        V0[0, i] = coefsSource[0] * \
            source.choice_timefct(sourceparameters, tshift-X[i]/Ceeff(0))
        V0[1, i] = coefsSource[1] * \
            source.choice_timefct(sourceparameters, tshift-X[i]/Ceeff(0))
    return V0

def initCauchyProblemHTL(configuration, frontiere, mat, sourceparameters,modulation):
    import matMod
    X, xmin, xmax, Nx, dx, ESIM, Npt = inputReading.geometric(configuration)
    Nmat, rho, cm = inputReading.material(mat)
    fs, CFL, tshift = inputReading.source(sourceparameters)
    alphaF = inputReading.frontiere(frontiere)
    rhox, Celx = inputReading.affectMaterials(
        alphaF, rho, cm, xmin, xmax, Nx, dx)
    if mat["Nat_-1"] == "solides":
        rhoeff, Ceff = matMod.propHomo(
            0, rhox[0], Celx[0], frontiere, modulation)
    if mat["Nat_-1"] == "solides modules":
        frhoeff, fCeff = matMod.propMatHomoForm(mat, modulation)
        rhoeff=rhox[0]*frhoeff
        Ceff=Celx[0]*fCeff
        t = sp.symbols("t")
    reff=sp.lambdify(t, rhoeff, modules='numpy')
    Ceeff=sp.lambdify(t, Ceff, modules='numpy')
    # dt=timeStep(CFL,dx,cm)
    coefsSource = np.zeros(2)
    if sourceparameters["Onde"] == "V":
        coefsSource[0] = 1/Ceeff(0)
        coefsSource[1] = -reff(0)
    if sourceparameters["Onde"] == "S":
        coefsSource[0] = 1/Ceeff(0)
        coefsSource[1] = -reff(0)
    V0 = np.zeros([2, Nx])
    for i in range(Nx):
        V0[0, i] = coefsSource[0] * \
            source.choice_timefct(sourceparameters, tshift-X[i]/Ceeff(0))
        V0[1, i] = coefsSource[1] * \
            source.choice_timefct(sourceparameters, tshift-X[i]/Ceeff(0))
    return V0

def initPointSourceProblem(configuration, frontiere, mat, sourceparameters):
    X, xmin, xmax, Nx, dx, ESIM, Npt = inputReading.geometric(configuration)
    # Nmat,rho,cm=inputReading.material(mat)
    # fs,CFL,tshift=inputReading.source(sourceparameters)
    # dt=timeStep(CFL,dx,cm)
    # lambdaOnde=cm/fs
    # xsmin=cm*tshift-lambdaOnde
    # xsmax=cm*tshift
    V0 = np.zeros([2, Nx], dtype='complex')
    return V0


############### Simulation parameters ##############

####################################################

#####################################################
### Settings of the integration scheme ###


def coefScheme(K):
    if K == 2:
        # print("Lax-Wendroff")
        gamma = np.zeros([K+1, K])
        gamma[0, 0] = -1/2
        gamma[0, 1] = -1/2
        gamma[1, 0] = 0
        gamma[1, 1] = 1
        gamma[2, 0] = 1/2
        gamma[2, 1] = -1/2
        # print(gamma)
    elif K == 4:
        # print("ADER4")
        gamma = np.zeros([K+1, K])
        gamma[0, 0] = 1/12
        gamma[0, 1] = 1/24
        gamma[0, 2] = -1/12
        gamma[0, 3] = -1/24
        gamma[1, 0] = -2/3
        gamma[1, 1] = -2/3
        gamma[1, 2] = 1/6
        gamma[1, 3] = 1/6
        gamma[2, 0] = 0
        gamma[2, 1] = 5/4
        gamma[2, 2] = 0
        gamma[2, 3] = -1/4
        gamma[3, 0] = 2/3
        gamma[3, 1] = -2/3
        gamma[3, 2] = -1/6
        gamma[3, 3] = 1/6
        gamma[4, 0] = -1/12
        gamma[4, 1] = 1/24
        gamma[4, 2] = 1/12
        gamma[4, 3] = -1/24
    elif np.floor(K/2) != K/2:
        print("Odd number not accepted")
    elif K > 4:
        gamma = np.zeros([K+1, K])
        for i in range(K+1):
            for j in range(K):
                texte = 'coef '+str(i)+','+str(j)+":"
                gamma[i, j] = float(input(texte))
    return gamma

# % Computation of CS


def CSeig(K, s, A, dt, dx):
    eigenvalues, eigenvectors = np.linalg.eig(A)
    V = eigenvectors
    Λ = np.diag(eigenvalues)
    Vm1 = np.linalg.inv(V)
    gamma = coefScheme(K)
    csm = np.zeros(np.shape(A))
    for m in range(K):
        csym = gamma[s, m]*(dt/dx)**(m+1)*np.linalg.matrix_power(Λ, m+1)
        csm = csm+csym
    cs = np.dot(np.dot(V, csm), Vm1)
    return cs


def CS(K, s, A, dt, dx):
    gamma = coefScheme(K)
    csm = np.zeros(np.shape(A))
    for m in range(K):
        csym = gamma[s, m]*(dt/dx)**(m+1)*np.linalg.matrix_power(A, m+1)
        csm = csm+csym
    # cs=np.dot(np.dot(V,csm),Vm1)
    return csm

def preTMat(frontiere, rho, cm, order, modulation):
    Ninter = int(frontiere["nombre de frontieres "])
    matCtot = np.empty((4, Ninter), dtype=object)
    for Ni in range(Ninter):
        namorK = "beta_"+str(Ni)
        namorM = "zeta_"+str(Ni)
        if frontiere.get(namorK)==None:
            etaK=0.
        else:
            etaK= float(frontiere[namorK])
        if frontiere.get(namorM)==None:
            etaM=0.
        else:
            etaM= float(frontiere[namorM])
        nalpha = "Alpha_"+str(Ni)
        alpha = float(frontiere[nalpha])
        contacti = "Contact_"+str(Ni)
        if frontiere[contacti] == "parfait":
            matCtot[0, Ni] = immergedInterfaces.matC(rho[Ni], cm[Ni], order)
            matCtot[1, Ni] = immergedInterfaces.matC(rho[Ni+1], cm[Ni+1], order)
            matCtot[2, Ni] = immergedInterfaces.matCInter(
                matCtot[0, Ni], matCtot[1, Ni])
        elif frontiere[contacti] == "masse-ressort":
            nameKni = "Kn_"+str(Ni)
            nameMni = "Mn_"+str(Ni)
            if float(frontiere[nameKni]) == 0:
                Kni = 0
            else:
                Kni = float(frontiere[nameKni])
            F = immergedInterfaces.SMMatrix(Kni, float(frontiere[nameMni]))
            matCtot[0, Ni] = immergedInterfaces.matC_SM(
                rho[Ni], cm[Ni], order, F,-etaK,-etaM)
            matCtot[1, Ni] = immergedInterfaces.matC_SM(
                rho[Ni+1], cm[Ni+1], order, -F,etaK,etaM)
            matCtot[2, Ni] = immergedInterfaces.matCInter(
                matCtot[0, Ni], matCtot[1, Ni])
        elif frontiere[contacti] == "masse-ressort module" and modulation["Synchro"] == "Yes":
            fctModulM, t, x, OmegaM, phiM = immergedInterfaces.fctModM(
                modulation)
            fctModulF, t, x, OmegaF, phiF = immergedInterfaces.fctModF(
                modulation)
            fctModulZ, t, x, OmegaZ, phiZ = immergedInterfaces.fctModZ(
                modulation)
            fctModulB, t, x, OmegaB, phiB = immergedInterfaces.fctModB(
                modulation)
            if modulation["TypeMod"] == "Square NS":
                ratio = sp.Float(modulation["Ratio"])
            ampF = sp.Float(modulation["DeltaF"])
            ampM = sp.Float(modulation["DeltaM"])
            FreqM = sp.Float(modulation["Freq"])
            FreqF = sp.Float(modulation["Freq"])
            celM = sp.Float(modulation["ModC"])
            celF = sp.Float(modulation["ModC"])
            if modulation["Dephasage0"] == "Yes":
                phiMm = sp.Float(modulation["phi0"])
                phiFm = sp.Float(modulation["phi0"])
            else:
                phiFm = 0
                phiMm = 0
            DF, DM = sp.symbols('DF DM')
            nameKni = "Kn_"+str(Ni)
            nameMni = "Mn_"+str(Ni)
            if sp.Float(frontiere[nameKni]) == 0.0:
                Fmod = 0
            else:
                Fmod = 1*(1+DF*fctModulF)/sp.Float(frontiere[nameKni])
            Mmod = (1+DM*fctModulM)*sp.Float(frontiere[nameMni])
            nDeltaZ="DeltaZ"
            nDeltaB="DeltaB"
            if modulation.get(nDeltaZ)==None:
                ampZ = 0.
                FreqZ = 0
            else:
                ampZ = sp.Float(modulation["DeltaZ"])
                FreqZ = sp.Float(modulation["Freq"])
            if modulation.get(nDeltaB)==None:
                ampB=0.
                FreqB = 0.
            else:
                ampB= sp.Float(modulation["DeltaB"])
                FreqB = sp.Float(modulation["Freq"])
            # ampB = sp.Float(modulation["DeltaB"])
            # celZ = sp.Float(modulation["ModC"])
            # celB = sp.Float(modulation["ModC"])
            if modulation["Dephasage0"] == "Yes":
                phiZm = sp.Float(modulation["phi0"])
                phiBm = sp.Float(modulation["phi0"])
            else:
                phiZm = 0
                phiBm = 0
            DZ, DB = sp.symbols('DZ DB')
            zmod = (1+DZ*fctModulZ)*etaM#sp.Float(frontiere[nameZni])
            bmod = (1+DB*fctModulB)*etaK#sp.Float(frontiere[nameBni])
            formalB = immergedInterfaces.CMMatrix(Fmod, Mmod)
            formalE = immergedInterfaces.dampMatrix(bmod,zmod)
            matCtot[0, Ni] = immergedInterfaces.formatC_SM_modulated(
                rho[Ni], cm[Ni], order, formalB, x, t, OmegaM, phiM, OmegaF, phiF, DF, DM, modulation,formalE)
            matCtot[1, Ni] = immergedInterfaces.formatC_SM_modulated(
                rho[Ni+1], cm[Ni+1], order, -formalB, x, t, OmegaM, phiM, OmegaF, phiF, DF, DM, modulation,-formalE)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaM, FreqM*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaF, FreqF*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DF, ampF)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DM, ampM)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaB, FreqB*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaZ, FreqZ*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DZ, ampZ)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DB, ampB)
            if modulation["TypeMod"] == "Square NS":
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiF, phiFm+ratio)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiM, phiMm+ratio)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiZ, phiZm+ratio)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiB, phiBm+ratio)
            else:
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiF, phiFm)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiM, phiMm)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiZ, phiZm)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiB, phiBm)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaM, FreqM*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaF, FreqF*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DF, ampF)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DM, ampM)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaB, FreqB*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaZ, FreqZ*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DZ, ampZ)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DB, ampB)
            if modulation["TypeMod"] == "Square NS":
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiF, phiFm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiM, phiMm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiZ, phiZm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiB, phiBm+ratio)
            else:
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiF, phiFm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiM, phiMm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiZ, phiZm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiB, phiBm)
        elif frontiere[contacti] == "masse-ressort module" and modulation["Synchro"] == "No":
            fctModulM, t, x, OmegaM, phiM = immergedInterfaces.fctModM(
                modulation)
            fctModulF, t, x, OmegaF, phiF = immergedInterfaces.fctModF(
                modulation)
            fctModulZ, t, x, OmegaZ, phiZ = immergedInterfaces.fctModZ(
                modulation)
            fctModulB, t, x, OmegaB, phiB = immergedInterfaces.fctModB(
                modulation)
            ampF = sp.Float(modulation["DeltaF"])
            ampM = sp.Float(modulation["DeltaM"])
            if modulation["TypeMod"] == "Square NS":
                ratio = sp.Float(modulation["Ratio"])
            FreqM = sp.Float(modulation["FreqM"])
            FreqF = sp.Float(modulation["FreqF"])
            celM = sp.Float(modulation["ModM"])
            celF = sp.Float(modulation["ModF"])
            if modulation["Dephasage0"] == "Yes":
                phiFm = sp.Float(modulation["phiF0"])
                phiMm = sp.Float(modulation["phiM0"])
            else:
                phiFm = 0
                phiMm = 0
            DF, DM = sp.symbols('DF DM')
            nameKni = "Kn_"+str(Ni)
            nameMni = "Mn_"+str(Ni)
            if sp.Float(frontiere[nameKni]) == 0:
                Fmod = 0
            else:
                Fmod = (1+DF*fctModulF)/sp.Float(frontiere[nameKni])
            Mmod = (1+DM*fctModulM)*sp.Float(frontiere[nameMni])
            nDeltaZ="DeltaZ"
            nDeltaB="DeltaB"
            if modulation.get(nDeltaZ)==None:
                ampZ = 0.
                FreqZ=0
            else:
                ampZ = sp.Float(modulation["DeltaZ"])
                FreqZ = sp.Float(modulation["FreqZ"])
            if modulation.get(nDeltaB)==None:
                ampB=0.
                FreqB=0
            else:
                ampB=  sp.Float(modulation["DeltaB"])
                FreqB = sp.Float(modulation["FreqB"])
            # ampB = sp.Float(modulation["DeltaB"])
            # celZ = sp.Float(modulation["ModC"])
            # celB = sp.Float(modulation["ModC"])
            if modulation["Dephasage0"] == "Yes":
                phiZm = sp.Float(modulation["phi0"])
                phiBm = sp.Float(modulation["phi0"])
            else:
                phiZm = 0
                phiBm = 0
            DZ, DB = sp.symbols('DZ DB')
            zmod = (1+DZ*fctModulZ)*etaM#sp.Float(frontiere[nameZni])
            bmod = (1+DB*fctModulB)*etaK#sp.Float(frontiere[nameBni])
            formalB = immergedInterfaces.CMMatrix(Fmod, Mmod)
            formalE = immergedInterfaces.dampMatrix(bmod,zmod)
            matCtot[0, Ni] = immergedInterfaces.formatC_SM_modulated(
                rho[Ni], cm[Ni], order, formalB, x, t, OmegaM, phiM, OmegaF, phiF, DF, DM, modulation,formalE)
            matCtot[1, Ni] = immergedInterfaces.formatC_SM_modulated(
                rho[Ni+1], cm[Ni+1], order, -formalB, x, t, OmegaM, phiM, OmegaF, phiF, DF, DM, modulation,-formalE)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaF, FreqF*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaM, FreqM*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DF, ampF)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DM, ampM)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaB, FreqB*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaZ, FreqZ*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DZ, ampZ)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DB, ampB)
            if modulation["TypeMod"] == "Square NS":
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiF, phiFm+ratio)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiM, phiMm+ratio)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiZ, phiZm+ratio)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiB, phiBm+ratio)
            else:
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiF, phiFm)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiM, phiMm)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiZ, phiZm)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiB, phiBm)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaF, FreqF*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaM, FreqM*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DF, ampF)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DM, ampM)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaB, FreqB*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaZ, FreqZ*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DZ, ampZ)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DB, ampB)
            if modulation["TypeMod"] == "Square NS":
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiF, phiFm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiM, phiMm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiZ, phiZm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiB, phiBm+ratio)
            else:
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiF, phiFm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiM, phiMm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiZ, phiZm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiB, phiBm)
        # matCtot[2,Ni]=immergedInterfaces.matCInter(matCtot[0,Ni], matCtot[1,Ni])
        elif frontiere[contacti] == "masse-ressort module R" and modulation["Synchro"] == "Yes":
            nameKni = "Kn_"+str(Ni)
            nameMni = "Mn_"+str(Ni)
            if float(frontiere[nameKni]) == 0.0:
                Kni = 0
            else:
                Kni = float(frontiere[nameKni])
            fctModulM, t, x, OmegaM, phiM = immergedInterfaces.fctModM(modulation)
            fctModulF, t, x, OmegaF, phiF = immergedInterfaces.fctModF(modulation)
            fctModulZ, t, x, OmegaZ, phiZ = immergedInterfaces.fctModZ(
                modulation)
            fctModulB, t, x, OmegaB, phiB = immergedInterfaces.fctModB(
                modulation)
            ampF = sp.Float(modulation["DeltaF"])
            ampM = sp.Float(modulation["DeltaM"])
            if modulation["TypeMod"] == "Square NS":
                ratio = sp.Float(modulation["Ratio"])
            if modulation["Dephasage"] == "x-ct":
                FreqM = -sp.Float(modulation["Freq"])
                FreqF = -sp.Float(modulation["Freq"])
            else:
                FreqM = sp.Float(modulation["Freq"])
                FreqF = sp.Float(modulation["Freq"])
            if modulation["Dephasage0"] == "Yes":
                phiF0 = sp.Float(modulation["phi0"])
                phiM0 = sp.Float(modulation["phi0"])
            else:
                phiF0 = 0
                phiM0 = 0
            Deltax=float(frontiere["Alpha_1"])-float(frontiere["Alpha_0"])
            cmil=cm[1]
            tau=Deltax/cmil
            phiMm = +tau*FreqM*2*np.pi
            phiFm = +tau*FreqF*2*np.pi
            DF, DM = sp.symbols('DF DM')
            if sp.Float(frontiere[nameKni]) == 0.0:
                Fmod = 0
            else:
                Fmod = (1+DF*fctModulF)/sp.Float(frontiere[nameKni])
            Mmod = (1+DM*fctModulM)*sp.Float(frontiere[nameMni])
            nDeltaZ="DeltaZ"
            nDeltaB="DeltaB"
            if modulation.get(nDeltaZ)==None:
                ampZ = 0.
                FreqZ = 0
            else:
                ampZ = sp.Float(modulation["DeltaZ"])
                FreqZ = sp.Float(modulation["Freq"])
            if modulation.get(nDeltaB)==None:
                ampB=0.
                FreqB = 0.
            else:
                ampB= sp.Float(modulation["DeltaB"])
                FreqB = sp.Float(modulation["Freq"])
            # ampB = sp.Float(modulation["DeltaB"])
            # celZ = sp.Float(modulation["ModC"])
            # celB = sp.Float(modulation["ModC"])
            if modulation["Dephasage0"] == "Yes":
                phiZm = sp.Float(modulation["phi0"])
                phiBm = sp.Float(modulation["phi0"])
            else:
                phiZm = 0
                phiBm = 0
            DZ, DB = sp.symbols('DZ DB')
            zmod = (1+DZ*fctModulZ)*etaM#sp.Float(frontiere[nameZni])
            bmod = (1+DB*fctModulB)*etaK#sp.Float(frontiere[nameBni])
            formalB = immergedInterfaces.CMMatrix(Fmod, Mmod)
            formalE = immergedInterfaces.dampMatrix(bmod,zmod)
            matCtot[0, Ni] = immergedInterfaces.formatC_SM_modulated(
                rho[Ni], cm[Ni], order, formalB, x, t, OmegaM, phiM, OmegaF, phiF, DF, DM, modulation,formalE)
            matCtot[1, Ni] = immergedInterfaces.formatC_SM_modulated(
                rho[Ni+1], cm[Ni+1], order, -formalB, x, t, OmegaM, phiM, OmegaF, phiF, DF, DM, modulation,-formalE)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaM, FreqM*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaF, FreqF*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DF, ampF)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DM, ampM)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaB, FreqB*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaZ, FreqZ*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DZ, ampZ)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DB, ampB)
            if modulation["TypeMod"] == "Square NS":
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiF, phiFm+ratio)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiM, phiMm+ratio)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiZ, phiZm+ratio)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiB, phiBm+ratio)
            else:
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiF, phiFm)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiM, phiMm)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiZ, phiZm)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiB, phiBm)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaM, FreqM*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaF, FreqF*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DF, ampF)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DM, ampM)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaB, FreqB*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaZ, FreqZ*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DZ, ampZ)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DB, ampB)
            if modulation["TypeMod"] == "Square NS":
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiF, phiFm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiM, phiMm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiZ, phiZm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiB, phiBm+ratio)
            else:
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiF, phiFm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiM, phiMm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiZ, phiZm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiB, phiBm)
        elif frontiere[contacti] == "masse-ressort module key" and modulation["Synchro"] == "Yes":
            nameKni = "Kn_"+str(Ni)
            nameMni = "Mn_"+str(Ni)
            if float(frontiere[nameKni]) == 0.0:
                Kni = 0
            else:
                Kni = float(frontiere[nameKni])
            fctModulM, t, x, OmegaM, phiM = immergedInterfaces.fctModM(modulation)
            fctModulF, t, x, OmegaF, phiF = immergedInterfaces.fctModF(modulation)
            fctModulZ, t, x, OmegaZ, phiZ = immergedInterfaces.fctModZ(
                modulation)
            fctModulB, t, x, OmegaB, phiB = immergedInterfaces.fctModB(
                modulation)
            ampF = sp.Float(modulation["DeltaF"])
            ampM = sp.Float(modulation["DeltaM"])
            if modulation["TypeMod"] == "Square NS":
                ratio = sp.Float(modulation["Ratio"])
            if modulation["Dephasage"] == "x-ct":
                FreqM = -sp.Float(modulation["Freq"])
                FreqF = -sp.Float(modulation["Freq"])
            else:
                FreqM = sp.Float(modulation["Freq"])
                FreqF = sp.Float(modulation["Freq"])
            if modulation["Dephasage0"] == "Yes":
                phiF0 = sp.Float(modulation["phi0"])
                phiM0 = sp.Float(modulation["phi0"])
            else:
                phiF0 = 0
                phiM0 = 0
            Deltax=float(frontiere["Alpha_1"])-float(frontiere["Alpha_0"])
            cmil=cm[1]
            tau=Deltax/cmil
            phiMm = +tau*FreqM*2*np.pi
            phiFm = +tau*FreqF*2*np.pi
            DF, DM = sp.symbols('DF DM')
            if sp.Float(frontiere[nameKni]) == 0.0:
                Fmod = 0
            else:
                Fmod = (1+DF*fctModulF)/sp.Float(frontiere[nameKni])
            Mmod = (1+DM*fctModulM)*sp.Float(frontiere[nameMni])
            nDeltaZ="DeltaZ"
            nDeltaB="DeltaB"
            if modulation.get(nDeltaZ)==None:
                ampZ = 0.
                FreqZ = 0
            else:
                ampZ = sp.Float(modulation["DeltaZ"])
                FreqZ = sp.Float(modulation["Freq"])
            if modulation.get(nDeltaB)==None:
                ampB=0.
                FreqB = 0.
            else:
                ampB= sp.Float(modulation["DeltaB"])
                FreqB = sp.Float(modulation["Freq"])
            # ampB = sp.Float(modulation["DeltaB"])
            # celZ = sp.Float(modulation["ModC"])
            # celB = sp.Float(modulation["ModC"])
            if modulation["Dephasage0"] == "Yes":
                phiZm = sp.Float(modulation["phi0"])
                phiBm = sp.Float(modulation["phi0"])
            else:
                phiZm = 0
                phiBm = 0
            DZ, DB = sp.symbols('DZ DB')
            zmod = (1+DZ*fctModulZ)*etaM#sp.Float(frontiere[nameZni])
            bmod = (1+DB*fctModulB)*etaK#sp.Float(frontiere[nameBni])
            formalB = immergedInterfaces.CMMatrix(Fmod, Mmod)
            formalE = immergedInterfaces.dampMatrix(bmod,zmod)
            matCtot[0, Ni] = immergedInterfaces.formatC_SM_modulated(
                rho[Ni], cm[Ni], order, formalB, x, t, OmegaM, phiM, OmegaF, phiF, DF, DM, modulation,formalE)
            matCtot[1, Ni] = immergedInterfaces.formatC_SM_modulated(
                rho[Ni+1], cm[Ni+1], order, -formalB, x, t, OmegaM, phiM, OmegaF, phiF, DF, DM, modulation,-formalE)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaM, FreqM*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaF, FreqF*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DF, ampF)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DM, ampM)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaB, FreqB*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaZ, FreqZ*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DZ, ampZ)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DB, ampB)
            if modulation["TypeMod"] == "Square NS":
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiF, phiFm+ratio)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiM, phiMm+ratio)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiZ, phiZm+ratio)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiB, phiBm+ratio)
            else:
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiF, phiFm)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiM, phiMm)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiZ, phiZm)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiB, phiBm)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaM, FreqM*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaF, FreqF*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DF, ampF)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DM, ampM)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaB, FreqB*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaZ, FreqZ*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DZ, ampZ)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DB, ampB)
            if modulation["TypeMod"] == "Square NS":
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiF, phiFm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiM, phiMm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiZ, phiZm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiB, phiBm+ratio)
            else:
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiF, phiFm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiM, phiMm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiZ, phiZm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiB, phiBm)
        elif frontiere[contacti] == "masse-ressort module dephase" and modulation["Synchro"] == "Yes":
            nameKni = "Kn_"+str(Ni)
            nameMni = "Mn_"+str(Ni)
            if float(frontiere[nameKni]) == 0.0:
                Kni = 0
            else:
                Kni = float(frontiere[nameKni])
            fctModulM, t, x, OmegaM, phiM = immergedInterfaces.fctModM(modulation)
            fctModulF, t, x, OmegaF, phiF = immergedInterfaces.fctModF(modulation)
            fctModulZ, t, x, OmegaZ, phiZ = immergedInterfaces.fctModZ(modulation)
            fctModulB, t, x, OmegaB, phiB = immergedInterfaces.fctModB(modulation)
            ampF = sp.Float(modulation["DeltaF"])
            ampM = sp.Float(modulation["DeltaM"])
            if modulation["TypeMod"] == "Square NS":
                ratio = sp.Float(modulation["Ratio"])
            if modulation["Dephasage"] == "x-ct":
                FreqM = -sp.Float(modulation["Freq"])
                FreqF = -sp.Float(modulation["Freq"])
                FreqZ = -sp.Float(modulation["Freq"])
                FreqB = -sp.Float(modulation["Freq"])
            else:
                FreqM = sp.Float(modulation["Freq"])
                FreqF = sp.Float(modulation["Freq"])
                FreqZ = sp.Float(modulation["Freq"])
                FreqB = sp.Float(modulation["Freq"])
            if modulation["Dephasage0"] == "Yes":
                phiF0 = sp.Float(modulation["phi0"])
                phiM0 = sp.Float(modulation["phi0"])
                phiZ0 = sp.Float(modulation["phi0"])
                phiB0 = sp.Float(modulation["phi0"])
            else:
                phiF0 = 0
                phiM0 = 0
                phiZ0 = 0
                phiB0 = 0
            if modulation["Dephasage"] == "xct":
                celM = sp.Float(modulation["ModC"])
                celF = sp.Float(modulation["ModC"])
                celZ = sp.Float(modulation["ModC"])
                celB = sp.Float(modulation["ModC"])
                if celM == 0:
                    phiMm = phiM0
                else:
                    k = -FreqM*2*np.pi/celM
                    phiMm = k*alpha+phiM0
                if celF == 0:
                    phiFm = phiF0
                else:
                    k = -FreqF*2*np.pi/celF
                    phiFm = k*alpha+phiF0
                if celZ == 0:
                    phiZm = phiZ0
                else:
                    k = -FreqZ*2*np.pi/celZ
                    phiZm = k*alpha+phiZ0
                if celB == 0:
                    phiBm = phiB0
                else:
                    k = -FreqB*2*np.pi/celB
                    phiBm = k*alpha+phiB0
            nDeltaZ="DeltaZ"
            nDeltaB="DeltaB"
            if modulation.get(nDeltaZ)==None:
                ampZ = 0.
            else:
                ampZ = sp.Float(modulation["DeltaZ"])
            if modulation.get(nDeltaB)==None:
                ampB = 0.
            else:
                ampB = sp.Float(modulation["DeltaB"])
            DF, DM = sp.symbols('DF DM')
            nameKni = "Kn_"+str(Ni)
            nameMni = "Mn_"+str(Ni)
            if sp.Float(frontiere[nameKni]) == 0.0:
                Fmod = 0
            else:
                Fmod = 1*(1+DF*fctModulF)/sp.Float(frontiere[nameKni])
            Mmod = (1+DM*fctModulM)*sp.Float(frontiere[nameMni])
            DZ, DB = sp.symbols('DZ DB')
            zmod = (1+DZ*fctModulZ)*etaM#sp.Float(frontiere[nameZni])
            bmod = (1+DB*fctModulB)*etaK#sp.Float(frontiere[nameBni])
            formalB = immergedInterfaces.CMMatrix(Fmod, Mmod)
            formalE = immergedInterfaces.dampMatrix(bmod,zmod)
            matCtot[0, Ni] = immergedInterfaces.formatC_SM_modulated(
                rho[Ni], cm[Ni], order, formalB, x, t, OmegaM, phiM, OmegaF, phiF, DF, DM, modulation,formalE)
            matCtot[1, Ni] = immergedInterfaces.formatC_SM_modulated(
                rho[Ni+1], cm[Ni+1], order, -formalB, x, t, OmegaM, phiM, OmegaF, phiF, DF, DM, modulation,-formalE)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaM, FreqM*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaF, FreqF*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DF, ampF)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DM, ampM)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaB, FreqB*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaZ, FreqZ*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DZ, ampZ)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DB, ampB)
            if modulation["TypeMod"] == "Square NS":
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiF, phiFm+ratio)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiM, phiMm+ratio)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiZ, phiZm+ratio)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiB, phiBm+ratio)
            else:
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiF, phiFm)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiM, phiMm)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiZ, phiZm)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiB, phiBm)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaM, FreqM*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaF, FreqF*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DF, ampF)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DM, ampM)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaB, FreqB*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaZ, FreqZ*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DZ, ampZ)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DB, ampB)
            if modulation["TypeMod"] == "Square NS":
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiF, phiFm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiM, phiMm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiZ, phiZm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiB, phiBm+ratio)
            else:
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiF, phiFm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiM, phiMm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiZ, phiZm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiB, phiBm)
        elif frontiere[contacti] == "masse-ressort module dephase" and modulation["Synchro"] == "No":
            ## À finir .... modulation de l'amortissement
            fctModulM, t, x, OmegaM, phiM = immergedInterfaces.fctModM(
                modulation)
            fctModulF, t, x, OmegaF, phiF = immergedInterfaces.fctModF(
                modulation)
            fctModulZ, t, x, OmegaZ, phiZ = immergedInterfaces.fctModZ(
                modulation)
            fctModulB, t, x, OmegaB, phiB = immergedInterfaces.fctModB(
                modulation)
            ampF = sp.Float(modulation["DeltaF"])
            ampM = sp.Float(modulation["DeltaM"])
            if modulation["TypeMod"] == "Square NS":
                ratio = sp.Float(modulation["Ratio"])
            if modulation["Dephasage"] == "x-ct":
                FreqM = -sp.Float(modulation["FreqM"])
                FreqF = -sp.Float(modulation["FreqF"])
            else:
                FreqM = sp.Float(modulation["FreqM"])
                FreqF = sp.Float(modulation["FreqF"])
            if modulation["Dephasage0"] == "Yes":
                phiF0 = sp.Float(modulation["phiF0"])
                phiM0 = sp.Float(modulation["phiM0"])
            else:
                phiF0 = 0
                phiM0 = 0
            if modulation["Dephasage"] == "xct":
                celM = sp.Float(modulation["ModM"])
                celF = sp.Float(modulation["ModF"])
                if celM == 0:
                    phiMm = phiM0
                else:
                    kM = -FreqM*2*np.pi/celM
                    phiMm = kM*alpha+phiM0
                if celF == 0:
                    phiFm = phiF0
                else:
                    kF = -FreqF*2*np.pi/celF
                    phiFm = kF*alpha+phiF0
            elif modulation["Dephasage"] == "t-x/c":
                celM = sp.Float(modulation["ModM"])
                celF = sp.Float(modulation["ModF"])
                if celM == 0:
                    phiMm = phiM0
                else:
                    phiMm = -FreqM*2*np.pi*alpha/celM+phiM0
                if celF == 0:
                    phiFm = phiF0
                else:
                    phiFm = -FreqF*2*np.pi*alpha/celF+phiF0
            DF, DM = sp.symbols('DF DM')
            nameKni = "Kn_"+str(Ni)
            nameMni = "Mn_"+str(Ni)
            if sp.Float(frontiere[nameKni]) == 0.0:
                Fmod = 0
            else:
                Fmod = (1+DF*fctModulF)/sp.Float(frontiere[nameKni])
            Mmod = (1+DM*fctModulM)*sp.Float(frontiere[nameMni])
            formalB = immergedInterfaces.CMMatrix(Fmod, Mmod)
            matCtot[0, Ni] = immergedInterfaces.formatC_SM_modulated(
                rho[Ni], cm[Ni], order, formalB, x, t, OmegaM, phiM, OmegaF, phiF, DF, DM, modulation,-etaK,-etaM)
            matCtot[1, Ni] = immergedInterfaces.formatC_SM_modulated(
                rho[Ni+1], cm[Ni+1], order, -formalB, x, t, OmegaM, phiM, OmegaF, phiF, DF, DM, modulation,etaK,etaM)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaM, FreqM*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaF, FreqF*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DF, ampF)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DM, ampM)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaB, FreqB*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaZ, FreqZ*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DZ, ampZ)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DB, ampB)
            if modulation["TypeMod"] == "Square NS":
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiF, phiFm+ratio)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiM, phiMm+ratio)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiZ, phiZm+ratio)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiB, phiBm+ratio)
            else:
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiF, phiFm)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiM, phiMm)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiZ, phiZm)
                matCtot[0, Ni] = matCtot[0, Ni].subs(phiB, phiBm)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaM, FreqM*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaF, FreqF*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DF, ampF)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DM, ampM)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaB, FreqB*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaZ, FreqZ*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DZ, ampZ)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DB, ampB)
            if modulation["TypeMod"] == "Square NS":
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiF, phiFm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiM, phiMm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiZ, phiZm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiB, phiBm+ratio)
            else:
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiF, phiFm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiM, phiMm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiZ, phiZm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiB, phiBm)
        # matCtot[2,Ni]=immergedInterfaces.matCInter(matCtot[0,Ni], matCtot[1,Ni])
    return matCtot
