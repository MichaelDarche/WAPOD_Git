
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

########### Integration scheme ############
# Finite differences without source term


def FD_cauchy(U0, Nt, Nx, K, rhox, cmx, dt, dx, film):
    U = U0.copy()
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    csTot = constructCS(K, Nx, rhox, cmx, dt, dx)
    for ts in range(Nt):
        # print(ts)
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Us = np.zeros([2, 1])
            # rhos=rhox[xstep]
            # cms=cmx[xstep]
            # A=inputReading.comA(rhos,cms)
            for s in range(K+1):
                cs = csTot[s, xstep]
                # cs=CS(K,s,A,dt,dx)
                sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                Us[0, 0] = Us[0, 0]+sumUs[0]
                Us[1, 0] = Us[1, 0]+sumUs[1]
            Un[0, xstep] = U[0, xstep]-Us[0, 0]
            Un[1, xstep] = U[1, xstep]-Us[1, 0]
        U = Un.copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "yes":
        return U, V, Sig
    else:
        return U


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
        elif frontiere[contacti] == "masse-ressort module independant":
            fctModulM, t, x, OmegaM, phiM = immergedInterfaces.fctModMidpt(
                frontiere,Ni)
            fctModulF, t, x, OmegaF, phiF = immergedInterfaces.fctModFidpt(
                frontiere,Ni)
            fctModulZ, t, x, OmegaZ, phiZ = immergedInterfaces.fctModZidpt(
                frontiere,Ni)
            fctModulB, t, x, OmegaB, phiB = immergedInterfaces.fctModBidpt(
                frontiere,Ni)
            DeltaFi="DeltaF_"+str(Ni)
            ampF = sp.Float(frontiere[DeltaFi])
            DeltaMi="DeltaM_"+str(Ni)
            ampM = sp.Float(frontiere[DeltaMi])
            TypeModi="TypeMod_"+str(Ni)
            if frontiere[TypeModi] == "Square NS":
                nRatio="Ratio_"+str(Ni)
                ratio = sp.Float(frontiere[nRatio])
            if frontiere["Synchro_"+str(Ni)]=="Yes":
                FreqMi="Freq_"+str(Ni)
                FreqFi="Freq_"+str(Ni)
                celM = sp.Float(frontiere["ModC_"+str(Ni)])
                celF = sp.Float(frontiere["ModC_"+str(Ni)])
            else:
                FreqMi="FreqM_"+str(Ni)
                FreqFi="FreqF_"+str(Ni)
                celM = sp.Float(frontiere["ModM_"+str(Ni)])
                celF = sp.Float(frontiere["ModF_"+str(Ni)])
            FreqM = sp.Float(frontiere[FreqMi])
            FreqF = sp.Float(frontiere[FreqFi])
            if frontiere["Dephasage0_"+str(Ni)] == "Yes":
                phiFm = sp.Float(frontiere["phiF0_"+str(Ni)])
                phiMm = sp.Float(frontiere["phiM0_"+str(Ni)])
            elif frontiere["Dephasage0_"+str(Ni)] == "x2Pi":
                phiFm = 2*sp.pi*sp.Float(frontiere["phiF0_"+str(Ni)])
                phiMm = 2*sp.pi*sp.Float(frontiere["phiM0_"+str(Ni)])
            else:
                phiFm = 0
                phiMm = 0
            DF, DM = sp.symbols('DF DM')
            nameKni = "Kn_"+str(Ni)
            nameMni = "Mn_"+str(Ni)
            if sp.Float(frontiere[nameKni]) == 0.0:
                Fmod = 0
            else:
                Fmod = (1+DF*fctModulF)/sp.Float(frontiere[nameKni])
            Mmod = (1+DM*fctModulM)*sp.Float(frontiere[nameMni])
            nDeltaZ="DeltaZ_"+str(Ni)
            nDeltaB="DeltaB_"+str(Ni)
            if frontiere.get(nDeltaZ)==None:
                ampZ = 0.
                FreqZ=0
            else:
                ampZ = sp.Float(frontiere["DeltaZ_"+str(Ni)])
                if frontiere["Synchro_"+str(Ni)]=="Yes":
                    FreqZ = sp.Float(frontiere["Freq_"+str(Ni)])
                else:
                    FreqZ = sp.Float(frontiere["FreqZ_"+str(Ni)])
            if frontiere.get(nDeltaB)==None:
                ampB=0.
                FreqB=0
            else:
                ampB=  sp.Float(frontiere["DeltaB_"+str(Ni)])
                if frontiere["Synchro_"+str(Ni)]=="Yes":
                    FreqB = sp.Float(frontiere["Freq_"+str(Ni)])
                else:
                    FreqB = sp.Float(frontiere["FreqB_"+str(Ni)])
            # ampB = sp.Float(modulation["DeltaB"])
            # celZ = sp.Float(modulation["ModC"])
            # celB = sp.Float(modulation["ModC"])
            if frontiere["Dephasage0_"+str(Ni)] == "Yes":
                phiZm = sp.Float(frontiere["phi0_"+str(Ni)])
                phiBm = sp.Float(frontiere["phi0_"+str(Ni)])
            elif frontiere["Dephasage0_"+str(Ni)] == "x2Pi":
                phiZm = 2*sp.pi*sp.Float(frontiere["phi0_"+str(Ni)])
                phiBm = 2*sp.pi*sp.Float(frontiere["phi0_"+str(Ni)])
            else:
                phiZm = 0
                phiBm = 0
            DZ, DB = sp.symbols('DZ DB')
            zmod = (1+DZ*fctModulZ)*etaM#sp.Float(frontiere[nameZni])
            bmod = (1+DB*fctModulB)*etaK#sp.Float(frontiere[nameBni])
            formalB = immergedInterfaces.CMMatrix(Fmod, Mmod)
            formalE = immergedInterfaces.dampMatrix(bmod,zmod)
            matCtot[0, Ni] = immergedInterfaces.formatC_SM_modulated_indep(
                rho[Ni], cm[Ni], order, formalB, x, t, OmegaM, phiM, OmegaF, phiF, DF, DM, frontiere,Ni,formalE)
            matCtot[1, Ni] = immergedInterfaces.formatC_SM_modulated_indep(
                rho[Ni+1], cm[Ni+1], order, -formalB, x, t, OmegaM, phiM, OmegaF, phiF, DF, DM, frontiere,Ni,-formalE)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaF, FreqF*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaM, FreqM*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DF, ampF)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DM, ampM)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaB, FreqB*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaZ, FreqZ*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DZ, ampZ)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DB, ampB)
            if frontiere["TypeMod_"+str(Ni)] == "Square NS":
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
            if frontiere["TypeMod_"+str(Ni)] == "Square NS":
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiF, phiFm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiM, phiMm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiZ, phiZm+ratio)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiB, phiBm+ratio)
            else:
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiF, phiFm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiM, phiMm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiZ, phiZm)
                matCtot[1, Ni] = matCtot[1, Ni].subs(phiB, phiBm)
        elif frontiere[contacti] == "masse-ressort module 2" and modulation["Synchro"] == "Yes":
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
            ampF = -sp.Float(modulation["DeltaF"])
            ampM = -sp.Float(modulation["DeltaM"])
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
                ampZ =- sp.Float(modulation["DeltaZ"])
                FreqZ = sp.Float(modulation["Freq"])
            if modulation.get(nDeltaB)==None:
                ampB=0.
                FreqB = 0.
            else:
                ampB= -sp.Float(modulation["DeltaB"])
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
        elif frontiere[contacti] == "masse-ressort module 2" and modulation["Synchro"] == "No":
            fctModulM, t, x, OmegaM, phiM = immergedInterfaces.fctModM(
                modulation)
            fctModulF, t, x, OmegaF, phiF = immergedInterfaces.fctModF(
                modulation)
            fctModulZ, t, x, OmegaZ, phiZ = immergedInterfaces.fctModZ(
                modulation)
            fctModulB, t, x, OmegaB, phiB = immergedInterfaces.fctModB(
                modulation)
            ampF = -sp.Float(modulation["DeltaF"])
            ampM = -sp.Float(modulation["DeltaM"])
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
                ampB=  -sp.Float(modulation["DeltaB"])
                FreqB = -sp.Float(modulation["FreqB"])
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


def preTMatK(frontiere, rho, cm, order, modulation):
    Ninter = int(frontiere["nombre de frontieres "])
    matCtot = np.empty((4, Ninter), dtype=object)
    for Ni in range(Ninter):
        nalpha = "Alpha_"+str(Ni)
        alpha = float(frontiere[nalpha])
        contacti = "Contact_"+str(Ni)
        if frontiere[contacti] == "parfait":
            matCtot[0, Ni] = immergedInterfaces.matC(rho[Ni], cm[Ni], order)
            matCtot[1, Ni] = immergedInterfaces.matC(
                rho[Ni+1], cm[Ni+1], order)
            matCtot[2, Ni] = immergedInterfaces.matCInter(
                matCtot[0, Ni], matCtot[1, Ni])
        if frontiere[contacti] == "masse-ressort":
            nameKni = "Kn_"+str(Ni)
            nameMni = "Mn_"+str(Ni)
            F = immergedInterfaces.SMMatrix(
                float(frontiere[nameKni]), float(frontiere[nameMni]))
            matCtot[0, Ni] = immergedInterfaces.matC_SM(
                rho[Ni], cm[Ni], order, F)
            matCtot[1, Ni] = immergedInterfaces.matC_SM(
                rho[Ni+1], cm[Ni+1], order, -F)
            matCtot[2, Ni] = immergedInterfaces.matCInter(
                matCtot[0, Ni], matCtot[1, Ni])
        if frontiere[contacti] == "masse-ressort module":
            fctModulM, t, x, OmegaM, phiM = immergedInterfaces.fctModM(
                modulation)
            fctModulF, t, x, OmegaF, phiF = immergedInterfaces.fctModF(
                modulation)
            ampF = sp.Float(modulation["DeltaF"])
            ampM = sp.Float(modulation["DeltaM"])
            FreqM = sp.Float(modulation["FreqM"])
            FreqF = sp.Float(modulation["FreqF"])
            celM = sp.Float(modulation["ModM"])
            celF = sp.Float(modulation["ModF"])
            if celM == 0:
                phiMm = 0
            else:
                phiMm = FreqM*2*np.pi*alpha/celM
            if celF == 0:
                phiFm = 0
            else:
                phiFm = FreqF*2*np.pi*alpha/celF
            DF, DM = sp.symbols('DF DM')
            nameKni = "Kn_"+str(Ni)
            nameMni = "Mn_"+str(Ni)
            Fmod = (1+DF*fctModulF)/sp.Float(frontiere[nameKni])
            Mmod = (1+DM*fctModulM)*sp.Float(frontiere[nameMni])
            formalB = immergedInterfaces.CMMatrix(Fmod, Mmod)
            matCtot[0, Ni] = immergedInterfaces.formatC_SM_modulated(
                rho[Ni], cm[Ni], order, formalB, x, t, OmegaM, phiM, OmegaF, phiF, DF, DM, modulation)
            matCtot[1, Ni] = immergedInterfaces.formatC_SM_modulated(
                rho[Ni+1], cm[Ni+1], order, -formalB, x, t, OmegaM, phiM, OmegaF, phiF, DF, DM, modulation)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaM, FreqM*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(OmegaF, FreqF*2*np.pi)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DF, ampF)
            matCtot[0, Ni] = matCtot[0, Ni].subs(DM, ampM)
            matCtot[0, Ni] = matCtot[0, Ni].subs(phiF, phiFm)
            matCtot[0, Ni] = matCtot[0, Ni].subs(phiM, phiMm)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaM, FreqM*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(OmegaF, FreqF*2*np.pi)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DF, ampF)
            matCtot[1, Ni] = matCtot[1, Ni].subs(DM, ampM)
            matCtot[1, Ni] = matCtot[1, Ni].subs(phiF, phiFm)
            matCtot[1, Ni] = matCtot[1, Ni].subs(phiM, phiMm)
            # matCtot[2,Ni]=immergedInterfaces.matCInter(matCtot[0,Ni], matCtot[1,Ni])
    return matCtot, Fmod, Mmod


def constructCS(K, Nx, rho, c, dt, dx):
    csTot = np.empty((K+1, Nx), dtype=object)
    for x in range(int(K/2), Nx-int(K/2)):
        for s in range(K+1):
            A = inputReading.comA(rho[x], c[x])
            csTot[s, x] = CS(K, s, A, dt, dx)
    return csTot

def constructCSWillis(K, Nx, rho, c,S1, dt, dx):
    csTot = np.empty((K+1, Nx), dtype=object)
    for x in range(int(K/2), Nx-int(K/2)):
        for s in range(K+1):
            A = inputReading.comAWillis(rho[x], c[x],S1[x])
            csTot[s, x] = CSeig(K, s, A, dt, dx)
    return csTot


def FD_cauchyII(U0, Nt, Nx, K, rhox, cmx, dt, dx, config, mat, frontiere, modulation, film,rCSVM="no",rtraces="no"):
    X, xmin, xmax, Nx, dx, ESIM, Npt = inputReading.geometric(config)
    order = 2*ESIM-1
    Nmat, rho, cm = inputReading.material(mat)
    Ninter = int(frontiere["nombre de frontieres "])
    alpha = inputReading.frontiere(frontiere)
    alpha[-1] = xmax+0.0001*xmax
    dim = 2
    U = U0.copy()
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    csTot = constructCS(K, Nx, rhox, cmx, dt, dx)
    formalCMat = preTMat(frontiere, rho, cm, order, modulation)
    t = sp.symbols('t')
    FMatrix = np.empty((4, Ninter), dtype=object)
    CSVM=np.zeros([Nt,4*Ninter])
    traces=np.zeros([Nt,4*Ninter])
    for Ni in range(Ninter):
        contacti = "Contact_"+str(Ni)
        if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
            FMatrix[:, Ni] = formalCMat[:, Ni]
        else:
            FMatrix[0, Ni] = sp.lambdify(t, formalCMat[0, Ni], modules="numpy")
            FMatrix[1, Ni] = sp.lambdify(t, formalCMat[1, Ni], modules="numpy")
    matCtot = np.empty((4, Ninter), dtype=object)
    for ts in range(Nt):
        print(ts)
        interS = 0
        # Un = np.zeros([2,Nx])
        for Ni in range(Ninter):
            contacti = "Contact_"+str(Ni)
            if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
                matCtot[:, Ni] = formalCMat[:, Ni]
            else:
                matCtot[0, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[0, Ni], t, ts*dt)
                matCtot[1, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[1, Ni], t, ts*dt)
                matCtot[2, Ni] = immergedInterfaces.matCInter(
                    matCtot[0, Ni], matCtot[1, Ni])
        # print(matCtot)
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            if (xstep+xmin/dx) <= alpha[interS]/dx-K/2:
                Us = np.zeros([2, 1])
                # rhos=rhox[xstep]
                # cms=cmx[xstep]
                # A=inputReading.comA(rhos,cms)
                for s in range(K+1):
                    cs = csTot[s, xstep]
                    # cs=CS(K,s,A,dt,dx)
                    sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                    Us[0, 0] = Us[0, 0]+sumUs[0]
                    Us[1, 0] = Us[1, 0]+sumUs[1]
                Un[0, xstep] = U[0, xstep]-Us[0, 0]
                Un[1, xstep] = U[1, xstep]-Us[1, 0]
            else:
                xj = int(alpha[interS]/dx)
                if xstep == xj-int(K/2)+1:
                    matCL = matCtot[0, interS]
                    matCR = matCtot[1, interS]
                    x = np.linspace(xj-Npt+1, xj+Npt, 2*Npt)*dx
                    mR = immergedInterfaces.matR(
                        dim, order, alpha[interS], x, Npt, matCL, matCR)
                    mL = np.dot(matCtot[2, interS], mR)
                    Uint = np.reshape(U[:, xj-Npt+1:xj+Npt+1],
                                      (dim*2*Npt, 1), order='F')
                    UmR = immergedInterfaces.modifiedValuesR(
                        U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                    UmL = immergedInterfaces.modifiedValuesL(
                        U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                    Ujum=-np.dot(mR,Uint)+np.dot(mL,Uint)
                    Umoy=np.dot(mR,Uint)+np.dot(mL,Uint)
                    CSVM[ts,0+interS*4]=Ujum[0,0]
                    CSVM[ts,1+interS*4]=Ujum[1,0]
                    CSVM[ts,2+interS*4]=1/2*Umoy[0,0]
                    CSVM[ts,3+interS*4]=1/2*Umoy[1,0]
                    tracesL=np.dot(mL,Uint)
                    tracesR=np.dot(mR,Uint)
                    traces[ts,0+interS*4]=tracesL[0,0]
                    traces[ts,1+interS*4]=tracesL[1,0]
                    traces[ts,2+interS*4]=tracesR[0,0]
                    traces[ts,3+interS*4]=tracesR[1,0]
                if K/2*dx > alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin <= alpha[interS]:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                elif xstep*dx+xmin > alpha[interS] and K/2*dx+alpha[interS] >= (xstep*dx+xmin):
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                else:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    if xstep == xj+int(K/2)+1:
                        interS = interS+1
        U = Un.copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    elif rCSVM=="no" and rtraces=="no":
        return U, V[::film], Sig[::film]
    elif rCSVM=="yes" and rtraces=="no":
        return U, V[::film], Sig[::film],CSVM
    elif rCSVM=="no" and rtraces=="yes":
        return U, V[::film], Sig[::film],traces
    else:
        return U, V[::film], Sig[::film],CSVM,traces

def FD_cauchyIIKey(U0, Nt, Nx, K, rhox, cmx, dt, dx, config, mat, frontiere, modulation, film,rCSVM="no",rtraces="no"):
    X, xmin, xmax, Nx, dx, ESIM, Npt = inputReading.geometric(config)
    order = 2*ESIM-1
    Nmat, rho, cm = inputReading.material(mat)
    Ninter = int(frontiere["nombre de frontieres "])
    alpha = inputReading.frontiere(frontiere)
    alpha[-1] = xmax+0.0001*xmax
    dim = 2
    U = U0.copy()
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    csTot = constructCS(K, Nx, rhox, cmx, dt, dx)
    formalCMat = preTMat(frontiere, rho, cm, order, modulation)
    t = sp.symbols('t')
    FMatrix = np.empty((4, Ninter), dtype=object)
    CSVM=np.zeros([Nt,4*Ninter])
    traces=np.zeros([Nt,4*Ninter])
    for Ni in range(Ninter):
        contacti = "Contact_"+str(Ni)
        if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
            FMatrix[:, Ni] = formalCMat[:, Ni]
        else:
            FMatrix[0, Ni] = sp.lambdify(t, formalCMat[0, Ni], modules="numpy")
            FMatrix[1, Ni] = sp.lambdify(t, formalCMat[1, Ni], modules="numpy")
    matCtot = np.empty((4, Ninter), dtype=object)
    for ts in range(Nt):
        print(Nt-ts)
        interS = 0
        # Un = np.zeros([2,Nx])
        for Ni in range(Ninter):
            contacti = "Contact_"+str(Ni)
            if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
                matCtot[:, Ni] = formalCMat[:, Ni]
            else:
                matCtot[0, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[0, Ni], t, -(Nt-ts)*dt)
                matCtot[1, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[1, Ni], t, -(Nt-ts)*dt)
                matCtot[2, Ni] = immergedInterfaces.matCInter(
                    matCtot[0, Ni], matCtot[1, Ni])
        # print(matCtot)
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            if (xstep+xmin/dx) <= alpha[interS]/dx-K/2:
                Us = np.zeros([2, 1])
                # rhos=rhox[xstep]
                # cms=cmx[xstep]
                # A=inputReading.comA(rhos,cms)
                for s in range(K+1):
                    cs = csTot[s, xstep]
                    # cs=CS(K,s,A,dt,dx)
                    sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                    Us[0, 0] = Us[0, 0]+sumUs[0]
                    Us[1, 0] = Us[1, 0]+sumUs[1]
                Un[0, xstep] = U[0, xstep]-Us[0, 0]
                Un[1, xstep] = U[1, xstep]-Us[1, 0]
            else:
                xj = int(alpha[interS]/dx)
                if xstep == xj-int(K/2)+1:
                    matCL = matCtot[0, interS]
                    matCR = matCtot[1, interS]
                    x = np.linspace(xj-Npt+1, xj+Npt, 2*Npt)*dx
                    mR = immergedInterfaces.matR(
                        dim, order, alpha[interS], x, Npt, matCL, matCR)
                    mL = np.dot(matCtot[2, interS], mR)
                    Uint = np.reshape(U[:, xj-Npt+1:xj+Npt+1],
                                      (dim*2*Npt, 1), order='F')
                    UmR = immergedInterfaces.modifiedValuesR(
                        U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                    UmL = immergedInterfaces.modifiedValuesL(
                        U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                    Ujum=-np.dot(mR,Uint)+np.dot(mL,Uint)
                    Umoy=np.dot(mR,Uint)+np.dot(mL,Uint)
                    CSVM[ts,0+interS*4]=Ujum[0,0]
                    CSVM[ts,1+interS*4]=Ujum[1,0]
                    CSVM[ts,2+interS*4]=1/2*Umoy[0,0]
                    CSVM[ts,3+interS*4]=1/2*Umoy[1,0]
                    tracesL=np.dot(mL,Uint)
                    tracesR=np.dot(mR,Uint)
                    traces[ts,0+interS*4]=tracesL[0,0]
                    traces[ts,1+interS*4]=tracesL[1,0]
                    traces[ts,2+interS*4]=tracesR[0,0]
                    traces[ts,3+interS*4]=tracesR[1,0]
                if K/2*dx > alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin <= alpha[interS]:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                elif xstep*dx+xmin > alpha[interS] and K/2*dx+alpha[interS] >= (xstep*dx+xmin):
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                else:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    if xstep == xj+int(K/2)+1:
                        interS = interS+1
        U = Un.copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    elif rCSVM=="no" and rtraces=="no":
        return U, V[::film], Sig[::film]
    elif rCSVM=="yes" and rtraces=="no":
        return U, V[::film], Sig[::film],CSVM
    elif rCSVM=="no" and rtraces=="yes":
        return U, V[::film], Sig[::film],traces
    else:
        return U, V[::film], Sig[::film],CSVM,traces
    
    
def FD_cauchyIIK(U0, Nt, Nx, K, rhox, cmx, dt, dx, config, mat, frontiere, modulation, film):
    X, xmin, xmax, Nx, dx, ESIM, Npt = inputReading.geometric(config)
    order = 2*ESIM-1
    Nmat, rho, cm = inputReading.material(mat)
    Ninter = int(frontiere["nombre de frontieres "])
    alpha = inputReading.frontiere(frontiere)
    alpha[-1] = xmax+0.0001*xmax
    dim = 2
    U = U0.copy()
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    csTot = constructCS(K, Nx, rhox, cmx, dt, dx)
    formalCMat = preTMatK(frontiere, rho, cm, order, modulation)
    t = sp.symbols('t')
    FMatrix = np.empty((4, Ninter), dtype=object)
    for Ni in range(Ninter):
        contacti = "Contact_"+str(Ni)
        if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
            FMatrix[:, Ni] = formalCMat[:, Ni]
        else:
            FMatrix[0, Ni] = sp.lambdify(t, formalCMat[0, Ni], modules="numpy")
            FMatrix[1, Ni] = sp.lambdify(t, formalCMat[1, Ni], modules="numpy")
    matCtot = np.empty((4, Ninter), dtype=object)
    for ts in range(Nt):
        print(ts)
        interS = 0
        # Un = np.zeros([2,Nx])
        for Ni in range(Ninter):
            contacti = "Contact_"+str(Ni)
            if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
                matCtot[:, Ni] = formalCMat[:, Ni]
            else:
                matCtot[0, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[0, Ni], t, ts*dt)
                matCtot[1, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[1, Ni], t, ts*dt)
                matCtot[2, Ni] = immergedInterfaces.matCInter(
                    matCtot[0, Ni], matCtot[1, Ni])
        # print(matCtot)
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            if (xstep+xmin/dx) <= alpha[interS]/dx-K/2:
                Us = np.zeros([2, 1])
                # rhos=rhox[xstep]
                # cms=cmx[xstep]
                # A=inputReading.comA(rhos,cms)
                for s in range(K+1):
                    cs = csTot[s, xstep]
                    # cs=CS(K,s,A,dt,dx)
                    sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                    Us[0, 0] = Us[0, 0]+sumUs[0]
                    Us[1, 0] = Us[1, 0]+sumUs[1]
                Un[0, xstep] = U[0, xstep]-Us[0, 0]
                Un[1, xstep] = U[1, xstep]-Us[1, 0]
            else:
                xj = int(alpha[interS]/dx)
                if xstep == xj-int(K/2)+1:
                    matCL = matCtot[0, interS]
                    matCR = matCtot[1, interS]
                    x = np.linspace(xj-Npt+1, xj+Npt, 2*Npt)*dx
                    mR = immergedInterfaces.matR(
                        dim, order, alpha[interS], x, Npt, matCL, matCR)
                    mL = np.dot(matCtot[2, interS], mR)
                    Uint = np.reshape(U[:, xj-Npt+1:xj+Npt+1],
                                      (dim*2*Npt, 1), order='F')
                    UmR = immergedInterfaces.modifiedValuesR(
                        U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                    UmL = immergedInterfaces.modifiedValuesL(
                        U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                if K/2*dx > alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin <= alpha[interS]:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                elif xstep*dx+xmin > alpha[interS] and K/2*dx+alpha[interS] >= (xstep*dx+xmin):
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)-periodic medium with homogeneous cells connected by imperfect
                        sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                else:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    if xstep == xj+int(K/2)+1:
                        interS = interS+1
        U = Un.copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]
# Finite differences with source term


def dirac_source(x, x_0, dx):
    if abs(x-x_0) < dx/2:
        return 1.0/dx
    elif abs(x-x_0) == dx/2:
        return 1.0/dx/2
    else:
        return 0.0

# def sce_dirac(sce_t,x,x_0,dx):
#     def sced(t,x):
#         sce=dirac_source(x, x_0, dx)*sce_t
#         return sce
#     return sced


def FD_sourceV(U0, Nt, Nx, K, rho, cm, dt, dx, x_0, sce_t, film):
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    csTot = constructCS(K, Nx, rho, cm, dt, dx)
    for ts in range(Nt):
        # print(ts)
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Us = np.zeros([2, 1])
            for s in range(K+1):
                cs = csTot[s, xstep]
                sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                Us[0, 0] = Us[0, 0]+sumUs[0]
                Us[1, 0] = Us[1, 0]+sumUs[1]
            Un[0, xstep] = U[0, xstep]-Us[0, 0]
            Un[1, xstep] = U[1, xstep]-Us[1, 0]
            # if xstep==int(x_0/dx):
            #     print(Un[:,xstep])
            Utn[0, xstep] = Un[0, xstep]+dt * dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
            
            Utn[1, xstep] = Un[1, xstep]
        print(sce_t[ts+1])
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]

def FD_sourceVWillis(U0, Nt, Nx, K, rho, cm,S1, dt, dx, x_0, sce_t, film):
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    csTot = constructCSWillis(K, Nx, rho, cm,S1, dt, dx)
    for ts in range(Nt):
        print(ts)
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Us = np.zeros([2, 1])
            for s in range(K+1):
                cs = csTot[s, xstep]
                sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                Us[0, 0] = Us[0, 0]+sumUs[0]
                Us[1, 0] = Us[1, 0]+sumUs[1]
            Un[0, xstep] = U[0, xstep]-Us[0, 0]
            Un[1, xstep] = U[1, xstep]-Us[1, 0]
            # if xstep==int(x_0/dx):
            #     print(Un[:,xstep])
            Utn[0, xstep] = Un[0, xstep]+dt * dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
            Utn[1, xstep] = Un[1, xstep]+dt * S1[xstep]* dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]


def FD_sourceVIIHomo0(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film):
    import homogenize
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    if frontiere["Contact_"+str(0)] == "masse-ressort":
        rhoeff, Ceff = homogenize.propHomo(
            0, rhox[0], cmx[0], frontiere, modulation)
    for ts in range(Nt):
        print(ts)
        if frontiere["Contact_"+str(0)] == "masse-ressort module":
            rhoeff, Ceff = homogenize.propHomo(
                ts*dt, rhox[0], cmx[0], frontiere, modulation)
            Ref = np.linspace(rhoeff, rhoeff, Nx)
            Cef = np.linspace(Ceff, Ceff, Nx)
        csTot = constructCS(K, Nx, Ref, Cef, dt, dx)
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Us = np.zeros([2, 1])
            for s in range(K+1):
                cs = csTot[s, xstep]
                sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                Us[0, 0] = Us[0, 0]+sumUs[0]
                Us[1, 0] = Us[1, 0]+sumUs[1]
            Un[0, xstep] = U[0, xstep]-Us[0, 0]
            Un[1, xstep] = U[1, xstep]-Us[1, 0]
            Utn[0, xstep] = Un[0, xstep]+dt * \
                dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[0]/rhoeff
            Utn[1, xstep] = Un[1, xstep]
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]


def comEP(rho, c):
    E = c**2*rho
    return E


def comP(rho, c):
    E = comEP(rho, c)
    P = np.array([[0, -rho], [-1/E, 0]])
    Pm1 = np.linalg.inv(P)
    return P, Pm1


def comPC(rho, c):
    E = comEP(rho, c)
    P = np.array([[rho, 0], [-0, 1/E]])
    Pm1 = np.linalg.inv(P)
    return P, Pm1


def comPF(rho, E):
    t = sp.symbols("t")
    P = sp.Matrix([[0, -1/E], [-rho, 0]])
    Pm1 = P.inv()
    dP = sp.diff(P, t)
    return Pm1, dP

def comPFa(rho, E,Qc,Qm):
    t = sp.symbols("t")
    P = sp.Matrix([[0, -1/E], [-rho, 0]])
    Pm1 = P.inv()
    dP = sp.diff(P, t)
    Q=sp.Matrix([[0, -Qc], [-Qm, 0]])
    return Pm1, dP,Q

def comPFadA(rho, E,Qc,Qm):
    t = sp.symbols("t")
    P = sp.Matrix([[0, -1/E], [-rho, 0]])
    Pm1 = P.inv()
    dP = sp.diff(P, t)
    dPm1 = sp.diff(Pm1, t)
    ddPm1 = sp.diff(dPm1, t)
    dddPm1 = sp.diff(ddPm1, t)
    Q=sp.Matrix([[0, -Qc], [-Qm, 0]])
    return Pm1, dP,Q,dPm1,ddPm1,dddPm1

def coefADERdA(dt,dx,Pm1, dP,Q,dPm1,ddPm1,dddPm1):
    alpha1=-dt/dx*(Pm1+dPm1*dt/2+ddPm1*dt**2/6+dddPm1*dt**3/24)
    alpha2=(dt/dx)**2*(Pm1*Pm1*1/2+(2*dPm1*Pm1+Pm1*dPm1)*dt/6+(3*ddPm1*Pm1+3*dPm1*dPm1+Pm1*ddPm1)*dt**2/24)
    alpha3=-(dt/dx)**3*(Pm1*Pm1*Pm1/6+(3*dPm1*Pm1*Pm1+2*Pm1*dPm1*Pm1+Pm1*Pm1*dPm1)*dt/24)
    alpha4=(dt/dx)**4*Pm1*Pm1*Pm1*Pm1/24
    coefsdA=[alpha1,alpha2,alpha3,alpha4]
    return coefsdA

def CSdA(K, s, A, dt, dx):
    gamma = coefScheme(K)
    csm = np.zeros(np.shape(A))
    for m in range(K):
        csym = gamma[s, m]*(dt/dx)**(m+1)*np.linalg.matrix_power(A, m+1)
        csm = csm+csym
    # cs=np.dot(np.dot(V,csm),Vm1)
    return csm

def constructCSP(K, Nx, rho, c, dt, dx):
    csTot = np.empty((K+1, Nx), dtype=object)
    for x in range(int(K/2), Nx-int(K/2)):
        for s in range(K+1):
            A, Pm1 = comP(rho[x], c[x])
            csTot[s, x] = CS(K, s, Pm1, dt, dx)
    return csTot


def constructCSGS(K, Nx, Ph, dt, dx):
    csTot = np.empty((K+1, Nx), dtype=object)
    for x in range(int(K/2), Nx-int(K/2)):
        for s in range(K+1):
            csTot[s, x] = CS(K, s, Ph, dt, dx)
    return csTot

def constructCSGSP(K, Nx, Ph, dt, dx):
    csTot = np.empty((K+1, Nx), dtype=object)
    for x in range(int(0), Nx):
        for s in range(K+1):
            csTot[s, x] = CS(K, s, Ph, dt, dx)
    return csTot

def constructCSC(K, Nx, rho, c, dt, dx):
    csTot = np.empty((K+1, Nx), dtype=object)
    for x in range(int(K/2), Nx-int(K/2)):
        for s in range(K+1):
            R = inputReading.comR(rho[x], c[x])
            csTot[s, x] = CS(K, s, R, dt, dx)
    return csTot


def FD_sourceVIICara(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film):
    X, xmin, xmax, Nx, dx, ESIM, Npt = inputReading.geometric(config)
    order = 2*ESIM-1
    Nmat, rho, cm = inputReading.material(mat)
    Ninter = int(frontiere["nombre de frontieres "])
    alpha = inputReading.frontiere(frontiere)
    alpha[-1] = xmax+0.0001*xmax
    dim = 2
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    csTot = constructCSC(K, Nx, rhox, cmx, dt, dx)
    formalCMat = preTMat(frontiere, 1/(rho*cm**2), cm, order, modulation)
    t = sp.symbols('t')
    FMatrix = np.empty((4, Ninter), dtype=object)
    for Ni in range(Ninter):
        contacti = "Contact_"+str(Ni)
        if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
            FMatrix[:, Ni] = formalCMat[:, Ni]
        else:
            FMatrix[0, Ni] = sp.lambdify(t, formalCMat[0, Ni], modules="numpy")
            FMatrix[1, Ni] = sp.lambdify(t, formalCMat[1, Ni], modules="numpy")
    matCtot = np.empty((4, Ninter), dtype=object)
    # formalCMat=immergedInterfaces.formatC_SM_modulated(rho,cm,order,formalB,x,t,Lambda,Omega,DF,DM)#,alpha,ts,LambdaM,OmegaM,ampF,ampM)
    P, Pm1 = comPC(rhox[0], cmx[0])
    for ts in range(Nt):
        print(ts)
        interS = 0
        # Un = np.zeros([2,Nx])
        for Ni in range(Ninter):
            contacti = "Contact_"+str(Ni)
            if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
                matCtot[:, Ni] = formalCMat[:, Ni]
            else:
                matCtot[0, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[0, Ni], t, ts*dt)
                matCtot[1, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[1, Ni], t, ts*dt)
                matCtot[2, Ni] = immergedInterfaces.matCInter(
                    matCtot[0, Ni], matCtot[1, Ni])
        # print(matCtot)
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            if xstep <= x_0:
                if (xstep+xmin/dx) <= alpha[interS]/dx-K/2:
                    Vs = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        Vm = np.dot(P, U[:, xstep+s-int(K/2)])
                        sumVs = np.dot(cs, Vm)
                        Vs[0, 0] = Vs[0, 0]+sumVs[0]
                        Vs[1, 0] = Vs[1, 0]+sumVs[1]
                    RU = np.dot(Pm1, Vs)
                    Un[0, xstep] = U[0, xstep]-RU[0, 0]
                    Un[1, xstep] = U[1, xstep]-RU[1, 0]
                    Utn[0, xstep] = Un[0, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = Un[1, xstep]
                else:
                    xj = int(alpha[interS]/dx)
                    if xstep == xj-int(K/2)+1:
                        matCL = matCtot[0, interS]
                        matCR = matCtot[1, interS]
                        x = np.linspace(xj-Npt+1, xj+Npt, 2*Npt)*dx
                        mR = immergedInterfaces.matR(
                            dim, order, alpha[interS], x, Npt, matCL, matCR)
                        mL = np.dot(matCtot[2, interS], mR)
                        # print(mL)
                        Vv = np.dot(P, U)
                        Uint = np.reshape(
                            Vv[:, xj-Npt+1:xj+Npt+1], (dim*2*Npt, 1), order='F')
                        UmR = immergedInterfaces.modifiedValuesR(
                            Vv, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                        UmL = immergedInterfaces.modifiedValuesL(
                            Vv, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                    if K/2*dx > alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin <= alpha[interS]:
                        Vs = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumVs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                            Vs[0, 0] = Vs[0, 0]+sumVs[0]
                            Vs[1, 0] = Vs[1, 0]+sumVs[1]
                        RU = np.dot(Pm1, Vs)
                        Un[0, xstep] = U[0, xstep]-RU[0, 0]
                        Un[1, xstep] = U[1, xstep]-RU[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]
                    elif xstep*dx+xmin > alpha[interS] and K/2*dx+alpha[interS] >= (xstep*dx+xmin):
                        Vs = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumVs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                            Vs[0, 0] = Vs[0, 0]+sumVs[0]
                            Vs[1, 0] = Vs[1, 0]+sumVs[1]
                        RU = np.dot(Pm1, Vs)
                        Un[0, xstep] = U[0, xstep]-RU[0, 0]
                        Un[1, xstep] = U[1, xstep]-RU[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]
                    else:
                        Vs = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            Vm = np.dot(P, U[:, xstep+s-int(K/2)])
                            sumVs = np.dot(cs, Vm)
                            Vs[0, 0] = Vs[0, 0]+sumVs[0]
                            Vs[1, 0] = Vs[1, 0]+sumVs[1]
                        RU = np.dot(Pm1, Vs)
                        Un[0, xstep] = U[0, xstep]-RU[0, 0]
                        Un[1, xstep] = U[1, xstep]-RU[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]
                        if xstep == xj+int(K/2)+1:
                            interS = interS+1
            if xstep > x_0:
                if (xstep+xmin/dx) < alpha[interS]/dx-K/2:
                    Vs = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        Vm = np.dot(P, U[:, xstep+s-int(K/2)])
                        sumVs = np.dot(cs, Vm)
                        Vs[0, 0] = Vs[0, 0]+sumVs[0]
                        Vs[1, 0] = Vs[1, 0]+sumVs[1]
                    RU = np.dot(Pm1, Vs)
                    Un[0, xstep] = U[0, xstep]-RU[0, 0]
                    Un[1, xstep] = U[1, xstep]-RU[1, 0]
                    Utn[0, xstep] = Un[0, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = Un[1, xstep]
                else:
                    xj = int(alpha[interS]/dx)
                    if xstep == xj-int(K/2):
                        matCL = matCtot[0, interS]
                        matCR = matCtot[1, interS]
                        x = np.linspace(xj-Npt, xj+Npt-1, 2*Npt)*dx
                        mR = immergedInterfaces.matR(
                            dim, order, alpha[interS], x, Npt, matCL, matCR)
                        mL = np.dot(matCtot[2, interS], mR)
                        # print(mL)
                        Vv = np.dot(P, U)
                        Uint = np.reshape(
                            Vv[:, xj-Npt:xj+Npt], (dim*2*Npt, 1), order='F')
                        UmR = immergedInterfaces.modifiedValuesRbis(
                            Vv, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                        UmL = immergedInterfaces.modifiedValuesLbis(
                            Vv, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                    if K/2*dx >= alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin < alpha[interS]:
                        Vs = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumVs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                            Vs[0, 0] = Vs[0, 0]+sumVs[0]
                            Vs[1, 0] = Vs[1, 0]+sumVs[1]
                        RU = np.dot(Pm1, Vs)
                        Un[0, xstep] = U[0, xstep]-RU[0, 0]
                        Un[1, xstep] = U[1, xstep]-RU[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]
                    elif xstep*dx+xmin >= alpha[interS] and K/2*dx+alpha[interS] > (xstep*dx+xmin):
                        Vs = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumVs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                            Vs[0, 0] = Vs[0, 0]+sumVs[0]
                            Vs[1, 0] = Vs[1, 0]+sumVs[1]
                        RU = np.dot(Pm1, Vs)
                        Un[0, xstep] = U[0, xstep]-RU[0, 0]
                        Un[1, xstep] = U[1, xstep]-RU[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]
                    else:
                        Vs = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            Vm = np.dot(P, U[:, xstep+s-int(K/2)])
                            sumVs = np.dot(cs, Vm)
                            Vs[0, 0] = V[0, 0]+sumVs[0]
                            Vs[1, 0] = Vs[1, 0]+sumVs[1]
                        RU = np.dot(Pm1, Vs)
                        Un[0, xstep] = U[0, xstep]-RU[0, 0]
                        Un[1, xstep] = U[1, xstep]-RU[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]
                        if xstep == xj+int(K/2)+1:
                            interS = interS+1
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]


def FD_sourceVIIHomoCara(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film):
    import homogenize
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    if frontiere["Contact_"+str(0)] == "masse-ressort":
        rhoeff, Ceff = homogenize.propHomo(
            0, rhox[0], cmx[0], frontiere, modulation)
    for ts in range(Nt):
        print(ts)
        if frontiere["Contact_"+str(0)] == "masse-ressort module":
            rhoeff, Ceff = homogenize.propHomo(
                ts*dt, rhox[0], cmx[0], frontiere, modulation)
            P, Pm1 = comP(rhoeff, Ceff)
            Ref = np.linspace(rhoeff, rhoeff, Nx)
            Cef = np.linspace(Ceff, Ceff, Nx)
        csTot = constructCSP(K, Nx, Ref, Cef, dt, dx)
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Vs = np.zeros([2, 1])
            for s in range(K+1):
                cs = csTot[s, xstep]
                Vm = np.dot(P, U[:, xstep+s-int(K/2)])
                sumVs = np.dot(cs, Vm)
                Vs[0, 0] = Vs[0, 0]+sumVs[0]
                Vs[1, 0] = Vs[1, 0]+sumVs[1]
            RU = np.dot(Pm1, Vs)
            # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Un[0, xstep] = U[0, xstep]-RU[0, 0]
            Un[1, xstep] = U[1, xstep]-RU[1, 0]
            Utn[0, xstep] = Un[0, xstep]+dt * \
                dirac_source(xstep*dx, x_0, dx) * \
                sce_t[ts+1]*rhox[int(x_0/dx)]/rhoeff
            Utn[1, xstep] = Un[1, xstep]
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]


def FD_sourceSIIHomoCara(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film):
    import homogenize
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    if frontiere["Contact_"+str(0)] == "masse-ressort":
        rhoeff, Ceff = homogenize.propHomo(
            0, rhox[0], cmx[0], frontiere, modulation)
    for ts in range(Nt):
        print(ts)
        if frontiere["Contact_"+str(0)] == "masse-ressort module":
            rhoeff, Ceff = homogenize.propHomo(
                ts*dt, rhox[0], cmx[0], frontiere, modulation)
            P, Pm1 = comP(rhoeff, Ceff)
            Ref = np.linspace(rhoeff, rhoeff, Nx)
            Cef = np.linspace(Ceff, Ceff, Nx)
        csTot = constructCSP(K, Nx, Ref, Cef, dt, dx)
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Vs = np.zeros([2, 1])
            for s in range(K+1):
                cs = csTot[s, xstep]
                Vm = np.dot(P, U[:, xstep+s-int(K/2)])
                sumVs = np.dot(cs, Vm)
                Vs[0, 0] = Vs[0, 0]+sumVs[0]
                Vs[1, 0] = Vs[1, 0]+sumVs[1]
            RU = np.dot(Pm1, Vs)
            # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Un[0, xstep] = U[0, xstep]-RU[0, 0]
            Un[1, xstep] = U[1, xstep]-RU[1, 0]
            Utn[0, xstep] = Un[0, xstep]
            Utn[1, xstep] = Un[1, xstep]+dt*dirac_source(
                xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/rhoeff*(cmx[int(x_0/dx)]/Ceff)**2
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]
#

#Godunov-Splitting for time interfaces
def FD_sourceVII_GS_HTL(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film):
    import matMod
    t = sp.symbols("t")
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    if mat["Nat_-1"] == "solides":
        S= np.zeros([2,2])
    if mat["Nat_-1"] == "solides modules":
        R0=rhox[0]
        facR,facC=matMod.propMatHomoForm(mat,modulation)
        fE=facR*facC*facC
        E0=rhox[0]*cmx[0]*cmx[0]
        Pm1, dP = comPF(R0*facR,E0*fE)
        if modulation["TypeMod"]=="Square":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="SquareSym":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Square NS":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Sinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus QP":
            S = Pm1*dP
        Sfunc = sp.lambdify(t, S, modules='numpy')
        Pm1func = sp.lambdify(t, Pm1, modules='numpy')
    for ts in range(int(Nt)):
        print(ts)
        tts=ts*dt
        Ph = Pm1func(tts)
        Ref = R0*np.linspace(float(facR.subs(t, tts)),float(facR.subs(t, tts)), Nx)
        csTot = constructCSGS(K, Nx, Ph, dt, dx)
        Ue = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Us = np.zeros([2, 1])
            for s in range(K+1):
                cs = csTot[s, xstep]
                sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                Us[0, 0] = Us[0, 0]+sumUs[0]
                Us[1, 0] = Us[1, 0]+sumUs[1]
                # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Ue[0, xstep] = U[0, xstep]-Us[0, 0]
            Ue[1, xstep] = U[1, xstep]-Us[1, 0]
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Un = np.zeros([2, Nx])
        # Ref=np.linspace(float(rhoeff.subs(t,ts*dt)),float(rhoeff.subs(t,ts*dt)),Nx)
        for xstep in range(int(K/2), Nx-int(K/2)):
            Sn = Sfunc(tts)
            expSn = [[np.exp(-dt*Sn[0, 0]), 0], [0, np.exp(-dt*Sn[1, 1])]]
            Un[:, xstep] = np.dot(Ue[:, xstep], expSn)
            # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Utn[0, xstep] = Un[0, xstep]+dt * \
                dirac_source(xstep*dx, x_0, dx) * \
                sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            Utn[1, xstep] = Un[1, xstep]
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]
    
    
#Godunov-Splitting for time interfaces   
def FD_sourceVIIHomoGS(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film):
    import homogenize
    t = sp.symbols("t")
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    if frontiere["Contact_"+str(0)] == "masse-ressort":
        rhoeff, Ceff = homogenize.propHomo(
            0, rhox[0], cmx[0], frontiere, modulation)
    if frontiere["Contact_"+str(0)] == "masse-ressort module":
        rhoeff, Eeff, Ceff,QCeff,QMeff = homogenize.propHomoForm(
            rhox[0], cmx[0], frontiere, modulation)
        Pm1, dP = comPF(rhoeff, Eeff)
        if modulation["TypeMod"]=="Square":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="SquareSym":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Square NS":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Sinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus QP":
            S = Pm1*dP
        Sfunc = sp.lambdify(t, S, modules='numpy')
        Pm1func = sp.lambdify(t, Pm1, modules='numpy')
    if frontiere["Contact_"+str(0)] == "masse-ressort module dephase":
        rhoeff, Eeff, Ceff = homogenize.propHomoForm(
            rhox[0], cmx[0], frontiere, modulation)
        Pm1, dP = comPF(rhoeff, Eeff)
        if modulation["TypeMod"]=="Square":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="SquareSym":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Square NS":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Sinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus QP":
            S = Pm1*dP
        Sfunc = sp.lambdify(t, S, modules='numpy')
        Pm1func = sp.lambdify(t, Pm1, modules='numpy')
    for ts in range(int(Nt)):
        print(ts)
        Ph = Pm1func(ts*dt)
        Ref = np.linspace(float(rhoeff.subs(t, ts*dt)),
                          float(rhoeff.subs(t, ts*dt)), Nx)
        csTot = constructCSGS(K, Nx, Ph, dt, dx)
        Ue = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Us = np.zeros([2, 1])
            for s in range(K+1):
                cs = csTot[s, xstep]
                sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                Us[0, 0] = Us[0, 0]+sumUs[0]
                Us[1, 0] = Us[1, 0]+sumUs[1]
                # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Ue[0, xstep] = U[0, xstep]-Us[0, 0]
            Ue[1, xstep] = U[1, xstep]-Us[1, 0]
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Un = np.zeros([2, Nx])
        # Ref=np.linspace(float(rhoeff.subs(t,ts*dt)),float(rhoeff.subs(t,ts*dt)),Nx)
        for xstep in range(int(K/2), Nx-int(K/2)):
            Sn = Sfunc(dt*ts)
            expSn = [[np.exp(-dt*Sn[0, 0]), 0], [0, np.exp(-dt*Sn[1, 1])]]
            Un[:, xstep] = np.dot(Ue[:, xstep], expSn)
            # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Utn[0, xstep] = Un[0, xstep]+dt * \
                dirac_source(xstep*dx, x_0, dx) * \
                sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            Utn[1, xstep] = Un[1, xstep]
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]
    

    
#Cauchy problem
def FD_CauchyVIIHomoGS(U0, Nt, Nx, K, rhox, cmx, dt, dx, config, mat, frontiere, modulation, film):
    import homogenize
    t = sp.symbols("t")
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    if frontiere["Contact_"+str(0)] == "masse-ressort":
        rhoeff, Ceff = homogenize.propHomo(
            0, rhox[0], cmx[0], frontiere, modulation)
    if frontiere["Contact_"+str(0)] == "masse-ressort module":
        rhoeff, Eeff, Ceff = homogenize.propHomoForm(
            rhox[0], cmx[0], frontiere, modulation)
        Pm1, dP = comPF(rhoeff, Eeff)
        if modulation["TypeMod"]=="Square":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="SquareSym":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Square NS":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Sinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus QP":
            S = Pm1*dP
        Sfunc = sp.lambdify(t, S, modules='numpy')
        Pm1func = sp.lambdify(t, Pm1, modules='numpy')
    if frontiere["Contact_"+str(0)] == "masse-ressort module dephase":
        rhoeff, Eeff, Ceff = homogenize.propHomoForm(
            rhox[0], cmx[0], frontiere, modulation)
        Pm1, dP = comPF(rhoeff, Eeff)
        if modulation["TypeMod"]=="Square":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="SquareSym":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Square NS":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Sinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus QP":
            S = Pm1*dP
        Sfunc = sp.lambdify(t, S, modules='numpy')
        Pm1func = sp.lambdify(t, Pm1, modules='numpy')
    for ts in range(int(Nt)):
        print(ts)
        Ph = Pm1func(ts*dt)
        #Ref = np.linspace(float(rhoeff.subs(t, ts*dt)),
        #                  float(rhoeff.subs(t, ts*dt)), Nx)
        csTot = constructCSGS(K, Nx, Ph, dt, dx)
        Ue = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Us = np.zeros([2, 1])
            for s in range(K+1):
                cs = csTot[s, xstep]
                sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                Us[0, 0] = Us[0, 0]+sumUs[0]
                Us[1, 0] = Us[1, 0]+sumUs[1]
                # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Ue[0, xstep] = U[0, xstep]-Us[0, 0]
            Ue[1, xstep] = U[1, xstep]-Us[1, 0]
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Un = np.zeros([2, Nx])
        # Ref=np.linspace(float(rhoeff.subs(t,ts*dt)),float(rhoeff.subs(t,ts*dt)),Nx)
        for xstep in range(int(K/2), Nx-int(K/2)):
            Sn = Sfunc(dt*ts)
            expSn = [[np.exp(-dt*Sn[0, 0]), 0], [0, np.exp(-dt*Sn[1, 1])]]
            Un[:, xstep] = np.dot(Ue[:, xstep], expSn)
            # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Utn[0, xstep] = Un[0, xstep]#+dt * \
               # dirac_source(xstep*dx, x_0, dx) * \
                #sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            Utn[1, xstep] = Un[1, xstep]
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]



def FD_CauchyVIIHomoGS_Period(U0, Nt, Nx, K, rhox, cmx, dt, dx, config, mat, frontiere, modulation, film):
    import homogenize
    t = sp.symbols("t")
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    if frontiere["Contact_"+str(0)] == "masse-ressort":
        rhoeff, Ceff = homogenize.propHomo(
            0, rhox[0], cmx[0], frontiere, modulation)
    if frontiere["Contact_"+str(0)] == "masse-ressort module":
        rhoeff, Eeff, Ceff = homogenize.propHomoForm(
            rhox[0], cmx[0], frontiere, modulation)
        Pm1, dP = comPF(rhoeff, Eeff)
        if modulation["TypeMod"]=="Square":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="SquareSym":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Square NS":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Sinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus QP":
            S = Pm1*dP
        Sfunc = sp.lambdify(t, S, modules='numpy')
        Pm1func = sp.lambdify(t, Pm1, modules='numpy')
    for ts in range(int(Nt)):
        print(ts)
        Ph = Pm1func(ts*dt)
        #Ref = np.linspace(float(rhoeff.subs(t, ts*dt)),
        #                  float(rhoeff.subs(t, ts*dt)), Nx)
        csTot = constructCSGSP(K, Nx, Ph, dt, dx)
        Ue = np.zeros([2, Nx])
        for xstep in range(int(0), Nx):
            Us = np.zeros([2, 1])
            for s in range(K+1):
                cs = csTot[s, xstep]
                if xstep+s-int(K/2)>=Nx:
                    sumUs = np.dot(cs, U[:, xstep-Nx+s-int(K/2)])
                elif xstep+s-int(K/2)<0:
                    sumUs = np.dot(cs, U[:, Nx+xstep+s-int(K/2)])
                else:
                    sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                Us[0, 0] = Us[0, 0]+sumUs[0]
                Us[1, 0] = Us[1, 0]+sumUs[1]
                # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Ue[0, xstep] = U[0, xstep]-Us[0, 0]
            Ue[1, xstep] = U[1, xstep]-Us[1, 0]
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Un = np.zeros([2, Nx])
        # Ref=np.linspace(float(rhoeff.subs(t,ts*dt)),float(rhoeff.subs(t,ts*dt)),Nx)
        for xstep in range(int(0), Nx):
            Sn = Sfunc(dt*ts)
            expSn = [[np.exp(-dt*Sn[0, 0]), 0], [0, np.exp(-dt*Sn[1, 1])]]
            Un[:, xstep] = np.dot(Ue[:, xstep], expSn)
            # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Utn[0, xstep] = Un[0, xstep]#+dt * \
               # dirac_source(xstep*dx, x_0, dx) * \
                #sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            Utn[1, xstep] = Un[1, xstep]
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]
    
    
def FD_CauchyVIIHomoGS_Period_HTL(U0, Nt, Nx, K, rhox, cmx, dt, dx, config, mat, modulation, film):
    import matMod
    t = sp.symbols("t")
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    if mat["Nat_-1"] == "solides":
        S= np.zeros([2,2])
    if mat["Nat_-1"] == "solides modules":
        R0=rhox[0]
        facR,facC=matMod.propMatHomoForm(mat,modulation)
        fE=facR*facC*facC
        E0=rhox[0]*cmx[0]*cmx[0]
        Pm1, dP = comPF(R0*facR,E0*fE)
        if modulation["TypeMod"]=="Square":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="SquareSym":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Square NS":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Sinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus QP":
            S = Pm1*dP
        Sfunc = sp.lambdify(t, S, modules='numpy')
        Pm1func = sp.lambdify(t, Pm1, modules='numpy')
    for ts in range(int(Nt)):
        print(ts)
        Ph = Pm1func(ts*dt)
        #Ref = np.linspace(float(rhoeff.subs(t, ts*dt)),
        #                  float(rhoeff.subs(t, ts*dt)), Nx)
        csTot = constructCSGSP(K, Nx, Ph, dt, dx)
        Ue = np.zeros([2, Nx])
        for xstep in range(int(0), Nx):
            Us = np.zeros([2, 1])
            for s in range(K+1):
                cs = csTot[s, xstep]
                if xstep+s-int(K/2)>=Nx:
                    sumUs = np.dot(cs, U[:, xstep-Nx+s-int(K/2)])
                elif xstep+s-int(K/2)<0:
                    sumUs = np.dot(cs, U[:, Nx+xstep+s-int(K/2)])
                else:
                    sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                Us[0, 0] = Us[0, 0]+sumUs[0]
                Us[1, 0] = Us[1, 0]+sumUs[1]
                # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Ue[0, xstep] = U[0, xstep]-Us[0, 0]
            Ue[1, xstep] = U[1, xstep]-Us[1, 0]
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Un = np.zeros([2, Nx])
        # Ref=np.linspace(float(rhoeff.subs(t,ts*dt)),float(rhoeff.subs(t,ts*dt)),Nx)
        for xstep in range(int(0), Nx):
            Sn = Sfunc(dt*ts)
            expSn = [[np.exp(-dt*Sn[0, 0]), 0], [0, np.exp(-dt*Sn[1, 1])]]
            Un[:, xstep] = np.dot(Ue[:, xstep], expSn)
            # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Utn[0, xstep] = Un[0, xstep]#+dt * \
               # dirac_source(xstep*dx, x_0, dx) * \
                #sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            Utn[1, xstep] = Un[1, xstep]
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]

def FD_sourceSIIHomoGS(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film):
    import homogenize
    t = sp.symbols("t")
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    if frontiere["Contact_"+str(0)] == "masse-ressort":
        rhoeff, Ceff = homogenize.propHomo(
            0, rhox[0], cmx[0], frontiere, modulation)
    if frontiere["Contact_"+str(0)] == "masse-ressort module":
        rhoeff, Eeff, Ceff = homogenize.propHomoForm(
            rhox[0], cmx[0], frontiere, modulation)
        Pm1, dP = comPF(rhoeff, Eeff)
        if modulation["TypeMod"]=="Square":
            S= Pm1*dP*0 
        if modulation["TypeMod"]=="SquareSym":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Square NS":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Sinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus QP":
            S = Pm1*dP
        Sfunc = sp.lambdify(t, S, modules='numpy')
        Pm1func = sp.lambdify(t, Pm1, modules='numpy')
    for ts in range(int(Nt)):
        print(ts)
        Ph = Pm1func(ts*dt)
        Ref = np.linspace(float(rhoeff.subs(t, ts*dt)),
                          float(rhoeff.subs(t, ts*dt)), Nx)
        Cef = np.linspace(float(Ceff.subs(t, (2*ts+1)*dt)),
                          float(Ceff.subs(t, (2*ts+1)*dt)), Nx)
        csTot = constructCSGS(K, Nx, Ph, dt, dx)
        Ue = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Us = np.zeros([2, 1])
            for s in range(K+1):
                cs = csTot[s, xstep]
                sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                Us[0, 0] = Us[0, 0]+sumUs[0]
                Us[1, 0] = Us[1, 0]+sumUs[1]
                # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Ue[0, xstep] = U[0, xstep]-Us[0, 0]
            Ue[1, xstep] = U[1, xstep]-Us[1, 0]
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Un = np.zeros([2, Nx])
        # Ref=np.linspace(float(rhoeff.subs(t,ts*dt)),float(rhoeff.subs(t,ts*dt)),Nx)
        for xstep in range(int(K/2), Nx-int(K/2)):
            Sn = Sfunc(dt*ts)
            expSn = [[np.exp(-dt*Sn[0, 0]), 0], [0, np.exp(-dt*Sn[1, 1])]]
            Un[:, xstep] = np.dot(Ue[:, xstep], expSn)
            # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Utn[0, xstep] = Un[0, xstep]
            Utn[1, xstep] = Un[1, xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(
                x_0/dx)]/Ref[xstep]*(cmx[int(x_0/dx)]/Cef[xstep])**2
            U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]


def FD_sourceVIIHomoSS(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film):
    import homogenize
    t = sp.symbols("t")
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    if frontiere["Contact_"+str(0)] == "masse-ressort":
        rhoeff, Ceff = homogenize.propHomo(
            0, rhox[0], cmx[0], frontiere, modulation)
    if frontiere["Contact_"+str(0)] == "masse-ressort module":
        rhoeff, Eeff, Ceff,QCeff,QMeff = homogenize.propHomoForm(rhox[0], cmx[0], frontiere, modulation)
        Pm1, dP = comPF(rhoeff, Eeff)
        if modulation["TypeMod"]=="Square":
            S= Pm1*dP*0 
        if modulation["TypeMod"]=="SquareSym":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Square NS":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Sinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus QP":
            S = Pm1*dP
        Sfunc = sp.lambdify(t, S, modules='numpy')
        Pm1func = sp.lambdify(t, Pm1, modules='numpy')
        rhofunc=sp.lambdify(t,rhoeff, modules='numpy')
        #Cfunc=sp.lambdify(t,Ceff, modules='numpy')
    if frontiere["Contact_"+str(0)] == "masse-ressort module dephase":
        rhoeff, Eeff, Ceff = homogenize.propHomoForm(rhox[0], cmx[0], frontiere, modulation)
        Pm1, dP = comPF(rhoeff, Eeff)
        if modulation["TypeMod"]=="Square":
            S= Pm1*dP*0 
        if modulation["TypeMod"]=="SquareSym":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Square NS":
            S= Pm1*dP*0 
        if modulation["TypeMod"]=="Sinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus QP":
            S = Pm1*dP
        Sfunc = sp.lambdify(t, S, modules='numpy')
        Pm1func = sp.lambdify(t,Pm1, modules='numpy')
        rhofunc=sp.lambdify(t,rhoeff, modules='numpy')
        #Cfunc=sp.lambdify(t,Ceff, modules='numpy')
    for ts in range(int(Nt)):
        print(ts)
        Ph = Pm1func(ts*dt)
        Ref = np.linspace(rhofunc(ts*dt),rhofunc(ts*dt), Nx)
        # Cef=np.linspace(float(Ceff.subs(t,(2*ts+1)*dt)),float(Ceff.subs(t,(2*ts+1)*dt)),Nx)
        csTot = constructCSGS(K, Nx, Ph, dt, dx)
        Ui = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Sn = Sfunc(dt*(ts-1/2))
            expSn = [[np.exp(-dt/2*Sn[0, 0]), 0], [0, np.exp(-dt/2*Sn[1, 1])]]
            Ui[:, xstep] = np.dot(U[:, xstep], expSn)
            # # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Ue = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Us = np.zeros([2, 1])
            for s in range(K+1):
                cs = csTot[s, xstep]
                sumUs = np.dot(cs, Ui[:, xstep+s-int(K/2)])
                Us[0, 0] = Us[0, 0]+sumUs[0]
                Us[1, 0] = Us[1, 0]+sumUs[1]
                # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Ue[0, xstep] = Ui[0, xstep]-Us[0, 0]
            Ue[1, xstep] = Ui[1, xstep]-Us[1, 0]
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Un = np.zeros([2, Nx])
        # Ref=np.linspace(float(rhoeff.subs(t,ts*dt)),float(rhoeff.subs(t,ts*dt)),Nx)
        for xstep in range(int(K/2), Nx-int(K/2)):
            Sn = Sfunc(dt*(ts+1/2))
            expSn = [[np.exp(-dt/2*Sn[0, 0]), 0], [0, np.exp(-dt/2*Sn[1, 1])]]
            Un[:, xstep] = np.dot(Ue[:, xstep], expSn)
            # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Utn[0, xstep] = Un[0, xstep]+dt * \
                dirac_source(xstep*dx, x_0, dx) * \
                sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            Utn[1, xstep] = Un[1, xstep]
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]

def FD_sourceVIIHomoSS2(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film):
    import homogenize
    t = sp.symbols("t")
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    if frontiere["Contact_"+str(0)] == "parfait":
        rhoeff, Ceff = homogenize.propHomo(
            0, rhox[0], cmx[0], frontiere, modulation)
        Pm1, dP = comPF(rhoeff, rhoeff*Ceff*Ceff)
        Sfunc = sp.lambdify(t, Pm1*dP*0, modules='numpy')
        Pm1func = sp.lambdify(t, Pm1, modules='numpy')
        rhofunc=sp.lambdify(t,rhoeff, modules='numpy')
    if frontiere["Contact_"+str(0)] == "masse-ressort":
        rhoeff, Ceff = homogenize.propHomo(
            0, rhox[0], cmx[0], frontiere, modulation)
    if frontiere["Contact_"+str(0)] == "masse-ressort module":
        rhoeff, Eeff, Ceff,QCeff,QMeff = homogenize.propHomoForm(rhox[0], cmx[0], frontiere, modulation)
        Pm1, dP = comPF(rhoeff, Eeff)
        if modulation["TypeMod"]=="Square":
            S= Pm1*dP*0 
        if modulation["TypeMod"]=="SquareSym":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Square NS":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Sinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus QP":
            S = Pm1*dP
        Sfunc = sp.lambdify(t, S, modules='numpy')
        Pm1func = sp.lambdify(t, Pm1, modules='numpy')
        rhofunc=sp.lambdify(t,rhoeff, modules='numpy')
    if frontiere["Contact_"+str(0)] == "masse-ressort module independant":
        rhoeff, Eeff, Ceff,QCeff,QMeff = homogenize.propHomoFormI(rhox[0], cmx[0], frontiere, nbpCell)
        Pm1, dP = comPF(rhoeff, Eeff)
        if modulation["TypeMod"]=="Square":
            S= Pm1*dP*0 
        if modulation["TypeMod"]=="SquareSym":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Square NS":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Sinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus QP":
            S = Pm1*dP
        Sfunc = sp.lambdify(t, S, modules='numpy')
        Pm1func = sp.lambdify(t, Pm1, modules='numpy')
        rhofunc=sp.lambdify(t,rhoeff, modules='numpy')
        #Cfunc=sp.lambdify(t,Ceff, modules='numpy')
    if frontiere["Contact_"+str(0)] == "masse-ressort module dephase":
        rhoeff, Eeff, Ceff = homogenize.propHomoForm(rhox[0], cmx[0], frontiere, modulation)
        Pm1, dP = comPF(rhoeff, Eeff)
        if modulation["TypeMod"]=="Square":
            S= Pm1*dP*0 
        if modulation["TypeMod"]=="SquareSym":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Square NS":
            S= Pm1*dP*0 
        if modulation["TypeMod"]=="Sinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus QP":
            S = Pm1*dP
        Sfunc = sp.lambdify(t, S, modules='numpy')
        Pm1func = sp.lambdify(t,Pm1, modules='numpy')
        rhofunc=sp.lambdify(t,rhoeff, modules='numpy')
        #Cfunc=sp.lambdify(t,Ceff, modules='numpy')
    for ts in range(int(Nt)):
        print(ts)
        Ph = Pm1func(ts*dt)
        Ref = rhofunc(ts*dt)
        # Cef=np.linspace(float(Ceff.subs(t,(2*ts+1)*dt)),float(Ceff.subs(t,(2*ts+1)*dt)),Nx)
        csTot = constructCSGS(K, Nx, Ph, dt, dx)
        Ui = np.zeros([2, Nx])
        Sn = Sfunc(dt*(ts))
        expSn = [[np.exp(-dt/2*Sn[0, 0]), 0], [0, np.exp(-dt/2*Sn[1, 1])]]
        for xstep in range(int(K/2), Nx-int(K/2)):
            Ui[:, xstep] = np.dot(U[:, xstep], expSn)
            # # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Ue = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Us = np.zeros([2, 1])
            for s in range(K+1):
                cs = csTot[s, xstep]
                sumUs = np.dot(cs, Ui[:, xstep+s-int(K/2)])
                Us[0, 0] = Us[0, 0]+sumUs[0]
                Us[1, 0] = Us[1, 0]+sumUs[1]
                # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Ue[0, xstep] = Ui[0, xstep]-Us[0, 0]
            Ue[1, xstep] = Ui[1, xstep]-Us[1, 0]
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Un = np.zeros([2, Nx])
        Sn = Sfunc(dt*(ts+1))
        expSn = [[np.exp(-dt/2*Sn[0, 0]), 0], [0, np.exp(-dt/2*Sn[1, 1])]]
        # Ref=np.linspace(float(rhoeff.subs(t,ts*dt)),float(rhoeff.subs(t,ts*dt)),Nx)
        for xstep in range(int(K/2), Nx-int(K/2)):
            Un[:, xstep] = np.dot(Ue[:, xstep], expSn)
            # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Utn[0, xstep] = Un[0, xstep]+dt * \
                dirac_source(xstep*dx, x_0, dx) * \
                sce_t[ts+1]*rhox[int(x_0/dx)]/Ref
            Utn[1, xstep] = Un[1, xstep]
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]

def FD_sourceVIIHomoSSP(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, Xmin,nbpCell, film):
    import homogenize
    t = sp.symbols("t")
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    rhoeff, Eeff, Ceff,QCeff,QMeff = homogenize.propHomoFormP(Xmin,mat,frontiere,nbpCell)
    print(rhoeff,Eeff,Ceff)
    Pm1, dP,Q = comPFa(rhoeff, Eeff,QCeff,QMeff)
    if frontiere["TypeMod_"+str(0)]=="Square":
        S= Pm1*(dP*0 +Q)
    if frontiere["TypeMod_"+str(0)]=="SquareSym":
        S= Pm1*(dP*0+Q)
    if frontiere["TypeMod_"+str(0)]=="Square NS":
        S= Pm1*(dP*0+Q)
    if frontiere["TypeMod_"+str(0)]=="Sinus":
        S = Pm1*(dP+Q)
    if frontiere["TypeMod_"+str(0)]=="Cosinus":
        S = Pm1*(dP+Q)
    if frontiere["TypeMod_"+str(0)]=="Cosinus QP":
        S = Pm1*(dP+Q)
    Sfunc = sp.lambdify(t, S, modules='numpy')
    Pm1func = sp.lambdify(t, Pm1, modules='numpy')
    rhofunc=sp.lambdify(t,rhoeff, modules='numpy')
    for ts in range(int(Nt)):
        print(ts)
        Ph = Pm1func(ts*dt)
        Ref = rhofunc(ts*dt)
        # Cef=np.linspace(float(Ceff.subs(t,(2*ts+1)*dt)),float(Ceff.subs(t,(2*ts+1)*dt)),Nx)
        csTot = constructCSGS(K, Nx, Ph, dt, dx)
        Ui = np.zeros([2, Nx])
        Sn = Sfunc(dt*(ts))
        expSn = [[np.exp(-dt/2*Sn[0, 0]), 0], [0, np.exp(-dt/2*Sn[1, 1])]]
        for xstep in range(int(K/2), Nx-int(K/2)):
            Ui[:, xstep] = np.dot(U[:, xstep], expSn)
            # # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Ue = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Us = np.zeros([2, 1])
            for s in range(K+1):
                cs = csTot[s, xstep]
                sumUs = np.dot(cs, Ui[:, xstep+s-int(K/2)])
                Us[0, 0] = Us[0, 0]+sumUs[0]
                Us[1, 0] = Us[1, 0]+sumUs[1]
                # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Ue[0, xstep] = Ui[0, xstep]-Us[0, 0]
            Ue[1, xstep] = Ui[1, xstep]-Us[1, 0]
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Un = np.zeros([2, Nx])
        Sn = Sfunc(dt*(ts+1))
        expSn = [[np.exp(-dt/2*Sn[0, 0]), 0], [0, np.exp(-dt/2*Sn[1, 1])]]
        # Ref=np.linspace(float(rhoeff.subs(t,ts*dt)),float(rhoeff.subs(t,ts*dt)),Nx)
        for xstep in range(int(K/2), Nx-int(K/2)):
            Un[:, xstep] = np.dot(Ue[:, xstep], expSn)
            # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Utn[0, xstep] = Un[0, xstep]+dt * \
                dirac_source(xstep*dx, x_0, dx) * \
                sce_t[ts+1]*rhox[int(x_0/dx)]/Ref
            Utn[1, xstep] = Un[1, xstep]
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]

def FD_sourceVIIHomoSSHdA(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, modulation, film):
    import matMod
    gamma = coefScheme(K)
    t = sp.symbols("t")
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    if mat["Nat_-1"] == "solides":
        S= np.zeros([2,2])
    if mat["Nat_-1"] == "solides modules ER" or mat["Nat_-1"] == "solides modules CR":
        R0=rhox[0]
        facR,facC=matMod.propMatHomoForm(mat,modulation)
        fE=facR*facC*facC
        E0=rhox[0]*cmx[0]*cmx[0]
        Pm1, dP1,Q,dP,ddP,dddP = comPFadA(R0*facR,E0*fE,0,0)
        alphasi=coefADERdA(dt,dx,Pm1, dP1,Q,dP,ddP,dddP)
    # Pm1, dP1,Q,dP,ddP,dddP = comPFadA(rhoeff, Eeff,QCeff,QMeff) 
    # alphasi=coefADERdA(dt,dx,rhoeff, Eeff,QCeff,QMeff)
    if modulation["TypeMod"]=="Square":
        S= Pm1*(dP1*0 +Q)
        alphasi=coefADERdA(dt,dx,Pm1, 0*dP1,Q,0*dP,0*ddP,0*dddP)
    elif modulation["TypeMod"]=="SquareSym":
        S= Pm1*(dP1*0+Q)
        alphasi=coefADERdA(dt,dx,Pm1, 0*dP1,Q,0*dP,0*ddP,0*dddP)
    elif modulation["TypeMod"]=="Square NS":
        S= Pm1*(dP1*0+Q)
        alphasi=coefADERdA(dt,dx,Pm1, 0*dP1,Q,0*dP,0*ddP,0*dddP)
    elif modulation["TypeMod"]=="Sinus":
        S = Pm1*(dP1+Q)
    elif modulation["TypeMod"]=="Cosinus":
        S = Pm1*(dP1+Q)
    elif modulation["TypeMod"]=="Cosinus QP":
        S = Pm1*(dP1+Q)
    Sfunc = sp.lambdify(t, S, modules='numpy')
    #Pm1func = sp.lambdify(t, Pm1, modules='numpy')
    rhofunc=sp.lambdify(t,R0*facR, modules='numpy')
    # Afunc0=alphasi[0]*gamma[0,0]+alphasi[1]*gamma[0,1]+alphasi[2]*gamma[0,2]+alphasi[3]*gamma[0,3]
    # Afunc1 =alphasi[0]*gamma[1,0]+alphasi[1]*gamma[1,1]+alphasi[2]*gamma[1,2]+alphasi[3]*gamma[1,3]
    # Afunc2=alphasi[0]*gamma[2,0]+alphasi[1]*gamma[2,1]+alphasi[2]*gamma[2,2]+alphasi[3]*gamma[2,1]
    # Afunc3=alphasi[0]*gamma[3,0]+alphasi[1]*gamma[3,1]+alphasi[2]*gamma[3,2]+alphasi[3]*gamma[3,3]
    # Afunc4=alphasi[0]*gamma[4,0]+alphasi[1]*gamma[4,1]+alphasi[2]*gamma[4,2]+alphasi[3]*gamma[4,3]
    Afunc0=alphasi[0]
    Afunc1 =alphasi[1]
    Afunc2=alphasi[2]
    Afunc3=alphasi[3]
    print(Afunc0)
    Am2 = sp.lambdify(t,Afunc0 , modules='numpy')
    Am1= sp.lambdify(t, Afunc1, modules='numpy')
    A0 = sp.lambdify(t,Afunc2, modules='numpy')
    A1 = sp.lambdify(t, Afunc3, modules='numpy')
    # A2 = sp.lambdify(t, Afunc4, modules='numpy')
    for ts in range(int(Nt)):
        print(ts)
        Atm2 = Am2(ts*dt)
        Atm1 = Am1(ts*dt)
        At0 = A0(ts*dt)
        At1 = A1(ts*dt)
        # At2 = A2(ts*dt)
        Ref = rhofunc(ts*dt)
        Ui = np.zeros([2, Nx])
        Sn = Sfunc(dt*(ts))
        expSn = [[np.exp(-dt/2*Sn[0, 0]), 0], [0, np.exp(-dt/2*Sn[1, 1])]]
        for xstep in range(int(K/2), Nx-int(K/2)):
            Ui[:, xstep] = np.dot(U[:, xstep], expSn)
            # # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Ue = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            # Us = np.zeros([2, 1])
            #Ue[:, xstep]=np.dot(Atm2,Ui[:,xstep-int(K/2)])+np.dot(Atm1,Ui[:,xstep-int(K/2+1)])+np.dot(At0,Ui[:,xstep])+np.dot(At1,Ui[:,xstep+1])+np.dot(At0,Ui[:,xstep+int(K/2)])
            dxU=(gamma[0,0]*Ui[:,xstep-int(K/2)]+gamma[1,0]*Ui[:,xstep-int(K/2)+1]+gamma[2,0]*Ui[:,xstep-int(K/2)+2]+gamma[3,0]*Ui[:,xstep-int(K/2)+3]+gamma[4,0]*Ui[:,xstep-int(K/2)+4])
            dxxU=-2*(gamma[0,1]*Ui[:,xstep-int(K/2)]+gamma[1,1]*Ui[:,xstep-int(K/2)+1]+gamma[2,1]*Ui[:,xstep-int(K/2)+2]+gamma[3,1]*Ui[:,xstep-int(K/2)+3]+gamma[4,1]*Ui[:,xstep-int(K/2)+4])
            dxxxU=6*(gamma[0,2]*Ui[:,xstep-int(K/2)]+gamma[1,2]*Ui[:,xstep-int(K/2)+1]+gamma[2,2]*Ui[:,xstep-int(K/2)+2]+gamma[3,2]*Ui[:,xstep-int(K/2)+3]+gamma[4,2]*Ui[:,xstep-int(K/2)+4])
            dxxxxU=-24*(gamma[0,3]*Ui[:,xstep-int(K/2)]+gamma[1,3]*Ui[:,xstep-int(K/2)+1]+gamma[2,3]*Ui[:,xstep-int(K/2)+2]+gamma[3,3]*Ui[:,xstep-int(K/2)+3]+gamma[4,3]*Ui[:,xstep-int(K/2)+4])
            Ue[:, xstep]=Ui[:,xstep]+np.dot(Atm2,dxU)+np.dot(Atm1,dxxU)+np.dot(At0,dxxxU)+np.dot(At1,dxxxxU)
            #     gamma[s,]
            #     sumUs = np.dot(cs, Ui[:, xstep+s-int(K/2)])
            #     Us[0, 0] = Us[0, 0]+sumUs[0]
            #     Us[1, 0] = Us[1, 0]+sumUs[1]
            #     # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            # Ue[0, xstep] = Ui[0, xstep]-Us[0, 0]
            # Ue[1, xstep] = Ui[1, xstep]-Us[1, 0]
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Un = np.zeros([2, Nx])
        Sn = Sfunc(dt*(ts+1))
        expSn = [[np.exp(-dt/2*Sn[0, 0]), 0], [0, np.exp(-dt/2*Sn[1, 1])]]
        # Ref=np.linspace(float(rhoeff.subs(t,ts*dt)),float(rhoeff.subs(t,ts*dt)),Nx)
        for xstep in range(int(K/2), Nx-int(K/2)):
            Un[:, xstep] = np.dot(Ue[:, xstep], expSn)
            # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Utn[0, xstep] = Un[0, xstep]+dt * \
                dirac_source(xstep*dx, x_0, dx) * \
                sce_t[ts+1]*rhox[int(x_0/dx)]/Ref
            Utn[1, xstep] = Un[1, xstep]
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]

def FD_sourceVIIHomoSSPdA(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, Xmin,nbpCell, film):
    import homogenize
    gamma = coefScheme(K)
    t = sp.symbols("t")
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    rhoeff, Eeff, Ceff,QCeff,QMeff = homogenize.propHomoFormP(Xmin,mat,frontiere,nbpCell)
    print(rhoeff,Eeff,Ceff)
    Pm1, dP1,Q,dP,ddP,dddP = comPFadA(rhoeff, Eeff,QCeff,QMeff) 
    alphasi=coefADERdA(dt,dx,Pm1, dP1,Q,dP,ddP,dddP)
    if frontiere["TypeMod_"+str(0)]=="Square":
        S= Pm1*(dP1*0 +Q)
    if frontiere["TypeMod_"+str(0)]=="SquareSym":
        S= Pm1*(dP1*0+Q)
    if frontiere["TypeMod_"+str(0)]=="Square NS":
        S= Pm1*(dP1*0+Q)
    if frontiere["TypeMod_"+str(0)]=="Sinus":
        S = Pm1*(dP1+Q)
    if frontiere["TypeMod_"+str(0)]=="Cosinus":
        S = Pm1*(dP1+Q)
    if frontiere["TypeMod_"+str(0)]=="Cosinus QP":
        S = Pm1*(dP1+Q)
    Sfunc = sp.lambdify(t, S, modules='numpy')
    #Pm1func = sp.lambdify(t, Pm1, modules='numpy')
    rhofunc=sp.lambdify(t,rhoeff, modules='numpy')
    # Afunc0=alphasi[0]*gamma[0,0]+alphasi[1]*gamma[0,1]+alphasi[2]*gamma[0,2]+alphasi[3]*gamma[0,3]
    # Afunc1 =alphasi[0]*gamma[1,0]+alphasi[1]*gamma[1,1]+alphasi[2]*gamma[1,2]+alphasi[3]*gamma[1,3]
    # Afunc2=alphasi[0]*gamma[2,0]+alphasi[1]*gamma[2,1]+alphasi[2]*gamma[2,2]+alphasi[3]*gamma[2,1]
    # Afunc3=alphasi[0]*gamma[3,0]+alphasi[1]*gamma[3,1]+alphasi[2]*gamma[3,2]+alphasi[3]*gamma[3,3]
    # Afunc4=alphasi[0]*gamma[4,0]+alphasi[1]*gamma[4,1]+alphasi[2]*gamma[4,2]+alphasi[3]*gamma[4,3]
    Afunc0=alphasi[0]
    Afunc1 =alphasi[1]
    Afunc2=alphasi[2]
    Afunc3=alphasi[3]
    Am2 = sp.lambdify(t,Afunc0 , modules='numpy')
    Am1= sp.lambdify(t, Afunc1, modules='numpy')
    A0 = sp.lambdify(t,Afunc2, modules='numpy')
    A1 = sp.lambdify(t, Afunc3, modules='numpy')
    # A2 = sp.lambdify(t, Afunc4, modules='numpy')
    for ts in range(int(Nt)):
        print(ts)
        Atm2 = Am2(ts*dt)
        Atm1 = Am1(ts*dt)
        At0 = A0(ts*dt)
        At1 = A1(ts*dt)
        # At2 = A2(ts*dt)
        Ref = rhofunc(ts*dt)
        Ui = np.zeros([2, Nx])
        Sn = Sfunc(dt*(ts))
        expSn = [[np.exp(-dt/2*Sn[0, 0]), 0], [0, np.exp(-dt/2*Sn[1, 1])]]
        for xstep in range(int(K/2), Nx-int(K/2)):
            Ui[:, xstep] = np.dot(U[:, xstep], expSn)
            # # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Ue = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            # Us = np.zeros([2, 1])
            #Ue[:, xstep]=np.dot(Atm2,Ui[:,xstep-int(K/2)])+np.dot(Atm1,Ui[:,xstep-int(K/2+1)])+np.dot(At0,Ui[:,xstep])+np.dot(At1,Ui[:,xstep+1])+np.dot(At0,Ui[:,xstep+int(K/2)])
            dxU=(gamma[0,0]*Ui[:,xstep-int(K/2)]+gamma[1,0]*Ui[:,xstep-int(K/2)+1]+gamma[2,0]*Ui[:,xstep-int(K/2)+2]+gamma[3,0]*Ui[:,xstep-int(K/2)+3]+gamma[4,0]*Ui[:,xstep-int(K/2)+4])
            dxxU=-2*(gamma[0,1]*Ui[:,xstep-int(K/2)]+gamma[1,1]*Ui[:,xstep-int(K/2)+1]+gamma[2,1]*Ui[:,xstep-int(K/2)+2]+gamma[3,1]*Ui[:,xstep-int(K/2)+3]+gamma[4,1]*Ui[:,xstep-int(K/2)+4])
            dxxxU=6*(gamma[0,2]*Ui[:,xstep-int(K/2)]+gamma[1,2]*Ui[:,xstep-int(K/2)+1]+gamma[2,2]*Ui[:,xstep-int(K/2)+2]+gamma[3,2]*Ui[:,xstep-int(K/2)+3]+gamma[4,2]*Ui[:,xstep-int(K/2)+4])
            dxxxxU=-24*(gamma[0,3]*Ui[:,xstep-int(K/2)]+gamma[1,3]*Ui[:,xstep-int(K/2)+1]+gamma[2,3]*Ui[:,xstep-int(K/2)+2]+gamma[3,3]*Ui[:,xstep-int(K/2)+3]+gamma[4,3]*Ui[:,xstep-int(K/2)+4])
            Ue[:, xstep]=Ui[:,xstep]+np.dot(Atm2,dxU)+np.dot(Atm1,dxxU)+np.dot(At0,dxxxU)+np.dot(At1,dxxxxU)
            #     gamma[s,]
            #     sumUs = np.dot(cs, Ui[:, xstep+s-int(K/2)])
            #     Us[0, 0] = Us[0, 0]+sumUs[0]
            #     Us[1, 0] = Us[1, 0]+sumUs[1]
            #     # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            # Ue[0, xstep] = Ui[0, xstep]-Us[0, 0]
            # Ue[1, xstep] = Ui[1, xstep]-Us[1, 0]
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Un = np.zeros([2, Nx])
        Sn = Sfunc(dt*(ts+1))
        expSn = [[np.exp(-dt/2*Sn[0, 0]), 0], [0, np.exp(-dt/2*Sn[1, 1])]]
        # Ref=np.linspace(float(rhoeff.subs(t,ts*dt)),float(rhoeff.subs(t,ts*dt)),Nx)
        for xstep in range(int(K/2), Nx-int(K/2)):
            Un[:, xstep] = np.dot(Ue[:, xstep], expSn)
            # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Utn[0, xstep] = Un[0, xstep]+dt * \
                dirac_source(xstep*dx, x_0, dx) * \
                sce_t[ts+1]*rhox[int(x_0/dx)]/Ref
            Utn[1, xstep] = Un[1, xstep]
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]
    
def FD_CauchyVIIHomoSS_Period(U0, Nt, Nx, K, rhox, cmx, dt, dx, config, mat, frontiere, modulation, film):
    import homogenize
    t = sp.symbols("t")
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    if frontiere["Contact_"+str(0)] == "masse-ressort":
        rhoeff, Ceff = homogenize.propHomo(
            0, rhox[0], cmx[0], frontiere, modulation)
    if frontiere["Contact_"+str(0)] == "masse-ressort module":
        rhoeff, Eeff, Ceff = homogenize.propHomoForm(rhox[0], cmx[0], frontiere, modulation)
        Pm1, dP = comPF(rhoeff, Eeff)
        if modulation["TypeMod"]=="Square":
            S= Pm1*dP*0 
        if modulation["TypeMod"]=="SquareSym":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Square NS":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Sinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus QP":
            S = Pm1*dP
        Sfunc = sp.lambdify(t, S, modules='numpy')
        Pm1func = sp.lambdify(t, Pm1, modules='numpy')
        rhofunc=sp.lambdify(t,rhoeff, modules='numpy')
        #Cfunc=sp.lambdify(t,Ceff, modules='numpy')
    if frontiere["Contact_"+str(0)] == "masse-ressort module dephase":
        rhoeff, Eeff, Ceff = homogenize.propHomoForm(rhox[0], cmx[0], frontiere, modulation)
        Pm1, dP = comPF(rhoeff, Eeff)
        if modulation["TypeMod"]=="Square":
            S= Pm1*dP*0 
        if modulation["TypeMod"]=="SquareSym":
            S= Pm1*dP*0
        if modulation["TypeMod"]=="Square NS":
            S= Pm1*dP*0 
        if modulation["TypeMod"]=="Sinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus":
            S = Pm1*dP
        if modulation["TypeMod"]=="Cosinus QP":
            S = Pm1*dP
        Sfunc = sp.lambdify(t, S, modules='numpy')
        Pm1func = sp.lambdify(t, Pm1, modules='numpy')
        rhofunc=sp.lambdify(t,rhoeff, modules='numpy')
        #Cfunc=sp.lambdify(t,Ceff, modules='numpy')
    for ts in range(int(Nt)):
        print(ts)
        Ph = Pm1func(ts*dt)
        Ref = np.linspace(rhofunc(ts*dt),
                          rhofunc(ts*dt), Nx)
        # Cef=np.linspace(float(Ceff.subs(t,(2*ts+1)*dt)),float(Ceff.subs(t,(2*ts+1)*dt)),Nx)
        csTot = constructCSGS(K, Nx, Ph, dt, dx)
        Ui = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Sn = Sfunc(dt*(ts-1/2))
            expSn = [[np.exp(-dt/2*Sn[0, 0]), 0], [0, np.exp(-dt/2*Sn[1, 1])]]
            Ui[:, xstep] = np.dot(U[:, xstep], expSn)
            # # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Ue = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Us = np.zeros([2, 1])
            for s in range(K+1):
                cs = csTot[s, xstep]
                if xstep+s-int(K/2)>=Nx:
                    sumUs = np.dot(cs, U[:, xstep-Nx+s-int(K/2)])
                elif xstep+s-int(K/2)<0:
                    sumUs = np.dot(cs, U[:, Nx+xstep+s-int(K/2)])
                else:
                    sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                Us[0, 0] = Us[0, 0]+sumUs[0]
                Us[1, 0] = Us[1, 0]+sumUs[1]
                # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Ue[0, xstep] = Ui[0, xstep]-Us[0, 0]
            Ue[1, xstep] = Ui[1, xstep]-Us[1, 0]
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Un = np.zeros([2, Nx])
        # Ref=np.linspace(float(rhoeff.subs(t,ts*dt)),float(rhoeff.subs(t,ts*dt)),Nx)
        for xstep in range(int(K/2), Nx-int(K/2)):
            Sn = Sfunc(dt*(ts+1/2))
            expSn = [[np.exp(-dt/2*Sn[0, 0]), 0], [0, np.exp(-dt/2*Sn[1, 1])]]
            Un[:, xstep] = np.dot(Ue[:, xstep], expSn)
            # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Utn[0, xstep] = Un[0, xstep]#+dt * \
                #dirac_source(xstep*dx, x_0, dx) * \
                #sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            Utn[1, xstep] = Un[1, xstep]
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]

def FD_sourceSIIHomoSS(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film):
    import homogenize
    t = sp.symbols("t")
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    if frontiere["Contact_"+str(0)] == "masse-ressort":
        rhoeff, Ceff = homogenize.propHomo(
            0, rhox[0], cmx[0], frontiere, modulation)
    if frontiere["Contact_"+str(0)] == "masse-ressort module":
        rhoeff, Eeff, Ceff = homogenize.propHomoForm(
            rhox[0], cmx[0], frontiere, modulation)
        Pm1, dP = comPF(rhoeff, Eeff)
        S = Pm1*dP
        Sfunc = sp.lambdify(t, S, modules='numpy')
        Pm1func = sp.lambdify(t, Pm1, modules='numpy')
    for ts in range(int(Nt)):
        print(ts)
        Ph = Pm1func(ts*dt)
        Ref = np.linspace(float(rhoeff.subs(t, ts*dt)),
                          float(rhoeff.subs(t, ts*dt)), Nx)
        Cef = np.linspace(float(Ceff.subs(t, (2*ts+1)*dt)),
                          float(Ceff.subs(t, (2*ts+1)*dt)), Nx)
        csTot = constructCSGS(K, Nx, Ph, dt, dx)
        Ui = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Sn = Sfunc(dt*(ts-1/2))
            expSn = [[np.exp(-dt/2*Sn[0, 0]), 0], [0, np.exp(-dt/2*Sn[1, 1])]]
            Ui[:, xstep] = np.dot(U[:, xstep], expSn)
            # # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Ue = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Us = np.zeros([2, 1])
            for s in range(K+1):
                cs = csTot[s, xstep]
                sumUs = np.dot(cs, Ui[:, xstep+s-int(K/2)])
                Us[0, 0] = Us[0, 0]+sumUs[0]
                Us[1, 0] = Us[1, 0]+sumUs[1]
                # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Ue[0, xstep] = Ui[0, xstep]-Us[0, 0]
            Ue[1, xstep] = Ui[1, xstep]-Us[1, 0]
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Un = np.zeros([2, Nx])
        # Ref=np.linspace(float(rhoeff.subs(t,ts*dt)),float(rhoeff.subs(t,ts*dt)),Nx)
        for xstep in range(int(K/2), Nx-int(K/2)):
            Sn = Sfunc(dt*(ts+1/2))
            expSn = [[np.exp(-dt/2*Sn[0, 0]), 0], [0, np.exp(-dt/2*Sn[1, 1])]]
            Un[:, xstep] = np.dot(Ue[:, xstep], expSn)
            Utn[0, xstep] = Un[0, xstep]
            Utn[1, xstep] = Un[1, xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(
                x_0/dx)]/Ref[xstep]*(cmx[int(x_0/dx)]/Cef[xstep])
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]


def FD_sourceVIIHomoSSC(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film):
    import homogenize
    t = sp.symbols("t")
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    if frontiere["Contact_"+str(0)] == "masse-ressort":
        rhoeff, Ceff = homogenize.propHomo(
            0, rhox[0], cmx[0], frontiere, modulation)
    if frontiere["Contact_"+str(0)] == "masse-ressort module":
        rhoeff, Eeff, Ceff = homogenize.propHomoForm(
            rhox[0], cmx[0], frontiere, modulation)
        Pm1, dP = comPF(rhoeff, Eeff)
        S = Pm1*dP
        Sfunc = sp.lambdify(t, S, modules='numpy')
        Pm1func = sp.lambdify(t, Pm1, modules='numpy')
    for ts in range(int(Nt)):
        print(ts)
        Ph = Pm1func(ts*dt)
        Ref = np.linspace(float(rhoeff.subs(t, ts*dt)),
                          float(rhoeff.subs(t, ts*dt)), Nx)
        # Cef=np.linspace(float(Ceff.subs(t,(2*ts+1)*dt)),float(Ceff.subs(t,(2*ts+1)*dt)),Nx)
        csTot = constructCSGS(K, Nx, Ph, dt, dx)
        Ui = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Sn1 = Sfunc(dt*(ts-1))
            Sn2 = Sfunc(dt*(ts))
            Sn = (Sn1+Sn2)/2
            expSn = [[np.exp(-dt/2*Sn[0, 0]), 0], [0, np.exp(-dt/2*Sn[1, 1])]]
            Ui[:, xstep] = np.dot(U[:, xstep], expSn)
            # # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Ue = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Us = np.zeros([2, 1])
            for s in range(K+1):
                cs = csTot[s, xstep]
                sumUs = np.dot(cs, Ui[:, xstep+s-int(K/2)])
                Us[0, 0] = Us[0, 0]+sumUs[0]
                Us[1, 0] = Us[1, 0]+sumUs[1]
                # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Ue[0, xstep] = Ui[0, xstep]-Us[0, 0]
            Ue[1, xstep] = Ui[1, xstep]-Us[1, 0]
            # Utn[0,xstep]=Un[0,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            # Utn[1,xstep]=Un[1,xstep]
        Un = np.zeros([2, Nx])
        # Ref=np.linspace(float(rhoeff.subs(t,ts*dt)),float(rhoeff.subs(t,ts*dt)),Nx)
        for xstep in range(int(K/2), Nx-int(K/2)):
            Sn1 = Sfunc(dt*(ts))
            Sn2 = Sfunc(dt*(ts+1))
            Sn = (Sn1+Sn2)/2
            expSn = [[np.exp(-dt/2*Sn[0, 0]), 0], [0, np.exp(-dt/2*Sn[1, 1])]]
            Un[:, xstep] = np.dot(Ue[:, xstep], expSn)
            # Un[:,xstep]=U[:,xstep]-np.dot(Pm1,Vs)
            Utn[0, xstep] = Un[0, xstep]+dt * \
                dirac_source(xstep*dx, x_0, dx) * \
                sce_t[ts+1]*rhox[int(x_0/dx)]/Ref[xstep]
            Utn[1, xstep] = Un[1, xstep]
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]


def FD_sourceS(U0, Nt, Nx, K, rho, cm, dt, dx, x_0, sce_t, film):
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    csTot = constructCS(K, Nx, rho, cm, dt, dx)
    for ts in range(Nt):
        # print(ts)
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            Us = np.zeros([2, 1])
            for s in range(K+1):
                cs = csTot[s, xstep]
                sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                Us[0, 0] = Us[0, 0]+sumUs[0]
                Us[1, 0] = Us[1, 0]+sumUs[1]
            Un[0, xstep] = U[0, xstep]-Us[0, 0]
            Un[1, xstep] = U[1, xstep]-Us[1, 0]
            Utn[0, xstep] = Un[0, xstep]
            Utn[1, xstep] = Un[1, xstep]+dt * \
                dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]


def FD_sourceVII(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film,rCSVM="no",rtraces="no",printT=10):
    X, xmin, xmax, Nx, dx, ESIM, Npt = inputReading.geometric(config)
    order = 2*ESIM-1
    Nmat, rho, cm = inputReading.material(mat)
    Ninter = int(frontiere["nombre de frontieres "])
    alpha = inputReading.frontiere(frontiere)
    alpha[-1] = xmax+0.0001*xmax
    dim = 2
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    csTot = constructCS(K, Nx, rhox, cmx, dt, dx)
    formalCMat = preTMat(frontiere, rho, cm, order, modulation)
    t = sp.symbols('t')
    FMatrix = np.empty((4, Ninter), dtype=object)
    CSVM=np.zeros([Nt,4*Ninter])
    traces=np.zeros([Nt,4*Ninter])
    for Ni in range(Ninter):
        contacti = "Contact_"+str(Ni)
        if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
            FMatrix[:, Ni] = formalCMat[:, Ni]
        else:
            FMatrix[0, Ni] = sp.lambdify(t, formalCMat[0, Ni], modules="numpy")
            FMatrix[1, Ni] = sp.lambdify(t, formalCMat[1, Ni], modules="numpy")
    matCtot = np.empty((4, Ninter), dtype=object)
    # formalCMat=immergedInterfaces.formatC_SM_modulated(rho,cm,order,formalB,x,t,Lambda,Omega,DF,DM)#,alpha,ts,LambdaM,OmegaM,ampF,ampM)
    for ts in range(Nt):
        if (ts)//printT==(ts)/printT:
            print(ts)
        interS = 0
        # Un = np.zeros([2,Nx])
        for Ni in range(Ninter):
            contacti = "Contact_"+str(Ni)
            if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
                matCtot[:, Ni] = formalCMat[:, Ni]
            else:
                matCtot[0, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[0, Ni], t, ts*dt)
                matCtot[1, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[1, Ni], t, ts*dt)
                matCtot[2, Ni] = immergedInterfaces.matCInter(
                    matCtot[0, Ni], matCtot[1, Ni])
        # print(matCtot)
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            if xstep <= x_0:
                if (xstep+xmin/dx) <= alpha[interS]/dx-K/2:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = Un[1, xstep]
                else:
                    xj = int(alpha[interS]/dx)
                    if xstep == xj-int(K/2)+1:
                        matCL = matCtot[0, interS]
                        matCR = matCtot[1, interS]
                        x = np.linspace(xj-Npt+1, xj+Npt, 2*Npt)*dx
                        mR = immergedInterfaces.matR(
                            dim, order, alpha[interS], x, Npt, matCL, matCR)
                        mL = np.dot(matCtot[2, interS], mR)
                        # print(mL)
                        Uint = np.reshape(
                            U[:, xj-Npt+1:xj+Npt+1], (dim*2*Npt, 1), order='F')
                        UmR = immergedInterfaces.modifiedValuesR(
                            U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                        UmL = immergedInterfaces.modifiedValuesL(
                            U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                        Ujum=-np.dot(mR,Uint)+np.dot(mL,Uint)
                        Umoy=np.dot(mR,Uint)+np.dot(mL,Uint)
                        CSVM[ts,0+interS*4]=Ujum[0,0]
                        CSVM[ts,1+interS*4]=Ujum[1,0]
                        CSVM[ts,2+interS*4]=1/2*Umoy[0,0]
                        CSVM[ts,3+interS*4]=1/2*Umoy[1,0]
                        tracesL=np.dot(mL,Uint)
                        tracesR=np.dot(mR,Uint)
                        traces[ts,0+interS*4]=tracesL[0,0]
                        traces[ts,1+interS*4]=tracesL[1,0]
                        traces[ts,2+interS*4]=tracesR[0,0]
                        traces[ts,3+interS*4]=tracesR[1,0]
                    if K/2*dx > alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin <= alpha[interS]:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]
                    elif xstep*dx+xmin > alpha[interS] and K/2*dx+alpha[interS] >= (xstep*dx+xmin):
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]
                    else:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]
                        if xstep == xj+int(K/2)+1:
                            interS = interS+1
            if xstep > x_0:
                if (xstep+xmin/dx) < alpha[interS]/dx-K/2:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = Un[1, xstep]
                else:
                    xj = int(alpha[interS]/dx)
                    if xstep == xj-int(K/2):
                        matCL = matCtot[0, interS]
                        matCR = matCtot[1, interS]
                        x = np.linspace(xj-Npt, xj+Npt-1, 2*Npt)*dx
                        mR = immergedInterfaces.matR(
                            dim, order, alpha[interS], x, Npt, matCL, matCR)
                        mL = np.dot(matCtot[2, interS], mR)
                        #print(Npt)
                        Uint = np.reshape(
                            U[:, xj-Npt:xj+Npt], (dim*2*Npt, 1), order='F')
                        UmR = immergedInterfaces.modifiedValuesRbis(
                            U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                        UmL = immergedInterfaces.modifiedValuesLbis(
                            U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                        Ujum=-np.dot(mR,Uint)+np.dot(mL,Uint)
                        Umoy=np.dot(mR,Uint)+np.dot(mL,Uint)
                        CSVM[ts,0+interS*4]=Ujum[0,0]
                        CSVM[ts,1+interS*4]=Ujum[1,0]
                        CSVM[ts,2+interS*4]=1/2*Umoy[0,0]
                        CSVM[ts,3+interS*4]=1/2*Umoy[1,0]
                        tracesL=np.dot(mL,Uint)
                        tracesR=np.dot(mR,Uint)
                        traces[ts,0+interS*4]=tracesL[0,0]
                        traces[ts,1+interS*4]=tracesL[1,0]
                        traces[ts,2+interS*4]=tracesR[0,0]
                        traces[ts,3+interS*4]=tracesR[1,0]
                    if K/2*dx >= alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin < alpha[interS]:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]
                    elif xstep*dx+xmin >= alpha[interS] and K/2*dx+alpha[interS] > (xstep*dx+xmin):
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]
                    else:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]
                        if xstep == xj+int(K/2)+1:
                            interS = interS+1
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    elif rCSVM=="no" and rtraces=="no":
        return U, V[::film], Sig[::film]
    elif rCSVM=="no"and rtraces=="yes":
        return U, V[::film], Sig[::film],traces
    elif rCSVM=="yes"and rtraces=="no":
        return U, V[::film], Sig[::film],CSVM
    else:
        return U, V[::film], Sig[::film],CSVM,traces
    
def FD_sourceVIItn(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film,rCSVM="no",rtraces="no",printT=10):
    X, xmin, xmax, Nx, dx, ESIM, Npt = inputReading.geometric(config)
    order = 2*ESIM-1
    Nmat, rho, cm = inputReading.material(mat)
    Ninter = int(frontiere["nombre de frontieres "])
    alpha = inputReading.frontiere(frontiere)
    alpha[-1] = xmax+0.0001*xmax
    dim = 2
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    csTot = constructCS(K, Nx, rhox, cmx, dt, dx)
    formalCMat = preTMat(frontiere, rho, cm, order, modulation)
    t = sp.symbols('t')
    FMatrix = np.empty((4, Ninter), dtype=object)
    CSVM=np.zeros([Nt,4*Ninter])
    traces=np.zeros([Nt,4*Ninter])
    for Ni in range(Ninter):
        contacti = "Contact_"+str(Ni)
        if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
            FMatrix[:, Ni] = formalCMat[:, Ni]
        else:
            FMatrix[0, Ni] = sp.lambdify(t, formalCMat[0, Ni], modules="numpy")
            FMatrix[1, Ni] = sp.lambdify(t, formalCMat[1, Ni], modules="numpy")
    matCtot = np.empty((4, Ninter), dtype=object)
    # formalCMat=immergedInterfaces.formatC_SM_modulated(rho,cm,order,formalB,x,t,Lambda,Omega,DF,DM)#,alpha,ts,LambdaM,OmegaM,ampF,ampM)
    for ts in range(Nt):
        if (ts)//printT==(ts)/printT:
            print(ts)
        interS = 0
        # Un = np.zeros([2,Nx])
        for Ni in range(Ninter):
            contacti = "Contact_"+str(Ni)
            if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
                matCtot[:, Ni] = formalCMat[:, Ni]
            else:
                matCtot[0, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[0, Ni], t, ts*dt)
                matCtot[1, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[1, Ni], t, ts*dt)
                matCtot[2, Ni] = immergedInterfaces.matCInter(
                    matCtot[0, Ni], matCtot[1, Ni])
        # print(matCtot)
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            if xstep <= x_0:
                if (xstep+xmin/dx) <= alpha[interS]/dx-K/2:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts]
                    Utn[1, xstep] = Un[1, xstep]
                else:
                    xj = int(alpha[interS]/dx)
                    if xstep == xj-int(K/2)+1:
                        matCL = matCtot[0, interS]
                        matCR = matCtot[1, interS]
                        x = np.linspace(xj-Npt+1, xj+Npt, 2*Npt)*dx
                        mR = immergedInterfaces.matR(
                            dim, order, alpha[interS], x, Npt, matCL, matCR)
                        mL = np.dot(matCtot[2, interS], mR)
                        # print(mL)
                        Uint = np.reshape(
                            U[:, xj-Npt+1:xj+Npt+1], (dim*2*Npt, 1), order='F')
                        UmR = immergedInterfaces.modifiedValuesR(
                            U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                        UmL = immergedInterfaces.modifiedValuesL(
                            U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                        Ujum=-np.dot(mR,Uint)+np.dot(mL,Uint)
                        Umoy=np.dot(mR,Uint)+np.dot(mL,Uint)
                        CSVM[ts,0+interS*4]=Ujum[0,0]
                        CSVM[ts,1+interS*4]=Ujum[1,0]
                        CSVM[ts,2+interS*4]=1/2*Umoy[0,0]
                        CSVM[ts,3+interS*4]=1/2*Umoy[1,0]
                        tracesL=np.dot(mL,Uint)
                        tracesR=np.dot(mR,Uint)
                        traces[ts,0+interS*4]=tracesL[0,0]
                        traces[ts,1+interS*4]=tracesL[1,0]
                        traces[ts,2+interS*4]=tracesR[0,0]
                        traces[ts,3+interS*4]=tracesR[1,0]
                    if K/2*dx > alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin <= alpha[interS]:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts]
                        Utn[1, xstep] = Un[1, xstep]
                    elif xstep*dx+xmin > alpha[interS] and K/2*dx+alpha[interS] >= (xstep*dx+xmin):
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts]
                        Utn[1, xstep] = Un[1, xstep]
                    else:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts]
                        Utn[1, xstep] = Un[1, xstep]
                        if xstep == xj+int(K/2)+1:
                            interS = interS+1
            if xstep > x_0:
                if (xstep+xmin/dx) < alpha[interS]/dx-K/2:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts]
                    Utn[1, xstep] = Un[1, xstep]
                else:
                    xj = int(alpha[interS]/dx)
                    if xstep == xj-int(K/2):
                        matCL = matCtot[0, interS]
                        matCR = matCtot[1, interS]
                        x = np.linspace(xj-Npt, xj+Npt-1, 2*Npt)*dx
                        mR = immergedInterfaces.matR(
                            dim, order, alpha[interS], x, Npt, matCL, matCR)
                        mL = np.dot(matCtot[2, interS], mR)
                        #print(Npt)
                        Uint = np.reshape(
                            U[:, xj-Npt:xj+Npt], (dim*2*Npt, 1), order='F')
                        UmR = immergedInterfaces.modifiedValuesRbis(
                            U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                        UmL = immergedInterfaces.modifiedValuesLbis(
                            U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                        Ujum=-np.dot(mR,Uint)+np.dot(mL,Uint)
                        Umoy=np.dot(mR,Uint)+np.dot(mL,Uint)
                        CSVM[ts,0+interS*4]=Ujum[0,0]
                        CSVM[ts,1+interS*4]=Ujum[1,0]
                        CSVM[ts,2+interS*4]=1/2*Umoy[0,0]
                        CSVM[ts,3+interS*4]=1/2*Umoy[1,0]
                        tracesL=np.dot(mL,Uint)
                        tracesR=np.dot(mR,Uint)
                        traces[ts,0+interS*4]=tracesL[0,0]
                        traces[ts,1+interS*4]=tracesL[1,0]
                        traces[ts,2+interS*4]=tracesR[0,0]
                        traces[ts,3+interS*4]=tracesR[1,0]
                    if K/2*dx >= alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin < alpha[interS]:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts]
                        Utn[1, xstep] = Un[1, xstep]
                    elif xstep*dx+xmin >= alpha[interS] and K/2*dx+alpha[interS] > (xstep*dx+xmin):
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts]
                        Utn[1, xstep] = Un[1, xstep]
                    else:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts]
                        Utn[1, xstep] = Un[1, xstep]
                        if xstep == xj+int(K/2)+1:
                            interS = interS+1
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    elif rCSVM=="no" and rtraces=="no":
        return U, V[::film], Sig[::film]
    elif rCSVM=="no"and rtraces=="yes":
        return U, V[::film], Sig[::film],traces
    elif rCSVM=="yes"and rtraces=="no":
        return U, V[::film], Sig[::film],CSVM
    else:
        return U, V[::film], Sig[::film],CSVM,traces
    
def FD_sourceVIIOneWay(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film,rCSVM="no",rtraces="no",printT=10):
    X, xmin, xmax, Nx, dx, ESIM, Npt = inputReading.geometric(config)
    order = 2*ESIM-1
    Nmat, rho, cm = inputReading.material(mat)
    Ninter = int(frontiere["nombre de frontieres "])
    alpha = inputReading.frontiere(frontiere)
    alpha[-1] = xmax+0.0001*xmax
    dim = 2
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    csTot = constructCS(K, Nx, rhox, cmx, dt, dx)
    formalCMat = preTMat(frontiere, rho, cm, order, modulation)
    t = sp.symbols('t')
    FMatrix = np.empty((4, Ninter), dtype=object)
    CSVM=np.zeros([Nt,4*Ninter])
    traces=np.zeros([Nt,4*Ninter])
    for Ni in range(Ninter):
        contacti = "Contact_"+str(Ni)
        if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
            FMatrix[:, Ni] = formalCMat[:, Ni]
        else:
            FMatrix[0, Ni] = sp.lambdify(t, formalCMat[0, Ni], modules="numpy")
            FMatrix[1, Ni] = sp.lambdify(t, formalCMat[1, Ni], modules="numpy")
    matCtot = np.empty((4, Ninter), dtype=object)
    # formalCMat=immergedInterfaces.formatC_SM_modulated(rho,cm,order,formalB,x,t,Lambda,Omega,DF,DM)#,alpha,ts,LambdaM,OmegaM,ampF,ampM)
    for ts in range(Nt):
        if (ts)//printT==(ts)/printT:
            print(ts)
        interS = 0
        # Un = np.zeros([2,Nx])
        for Ni in range(Ninter):
            contacti = "Contact_"+str(Ni)
            if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
                matCtot[:, Ni] = formalCMat[:, Ni]
            else:
                matCtot[0, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[0, Ni], t, ts*dt)
                matCtot[1, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[1, Ni], t, ts*dt)
                matCtot[2, Ni] = immergedInterfaces.matCInter(
                    matCtot[0, Ni], matCtot[1, Ni])
        # print(matCtot)
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            if xstep <= x_0:
                if (xstep+xmin/dx) <= alpha[interS]/dx-K/2:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = Un[1, xstep]-rhox[xstep]*cmx[xstep]*dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                else:
                    xj = int(alpha[interS]/dx)
                    if xstep == xj-int(K/2)+1:
                        matCL = matCtot[0, interS]
                        matCR = matCtot[1, interS]
                        x = np.linspace(xj-Npt+1, xj+Npt, 2*Npt)*dx
                        mR = immergedInterfaces.matR(
                            dim, order, alpha[interS], x, Npt, matCL, matCR)
                        mL = np.dot(matCtot[2, interS], mR)
                        # print(mL)
                        Uint = np.reshape(
                            U[:, xj-Npt+1:xj+Npt+1], (dim*2*Npt, 1), order='F')
                        UmR = immergedInterfaces.modifiedValuesR(
                            U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                        UmL = immergedInterfaces.modifiedValuesL(
                            U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                        Ujum=-np.dot(mR,Uint)+np.dot(mL,Uint)
                        Umoy=np.dot(mR,Uint)+np.dot(mL,Uint)
                        CSVM[ts,0+interS*4]=Ujum[0,0]
                        CSVM[ts,1+interS*4]=Ujum[1,0]
                        CSVM[ts,2+interS*4]=1/2*Umoy[0,0]
                        CSVM[ts,3+interS*4]=1/2*Umoy[1,0]
                        tracesL=np.dot(mL,Uint)
                        tracesR=np.dot(mR,Uint)
                        traces[ts,0+interS*4]=tracesL[0,0]
                        traces[ts,1+interS*4]=tracesL[1,0]
                        traces[ts,2+interS*4]=tracesR[0,0]
                        traces[ts,3+interS*4]=tracesR[1,0]
                    if K/2*dx > alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin <= alpha[interS]:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]-rhox[xstep]*cmx[xstep]*dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    elif xstep*dx+xmin > alpha[interS] and K/2*dx+alpha[interS] >= (xstep*dx+xmin):
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]-rhox[xstep]*cmx[xstep]*dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    else:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]-rhox[xstep]*cmx[xstep]*dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        if xstep == xj+int(K/2)+1:
                            interS = interS+1
            if xstep > x_0:
                if (xstep+xmin/dx) < alpha[interS]/dx-K/2:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = Un[1, xstep]-rhox[xstep]*cmx[xstep]*dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                else:
                    xj = int(alpha[interS]/dx)
                    if xstep == xj-int(K/2):
                        matCL = matCtot[0, interS]
                        matCR = matCtot[1, interS]
                        x = np.linspace(xj-Npt, xj+Npt-1, 2*Npt)*dx
                        mR = immergedInterfaces.matR(
                            dim, order, alpha[interS], x, Npt, matCL, matCR)
                        mL = np.dot(matCtot[2, interS], mR)
                        #print(Npt)
                        Uint = np.reshape(
                            U[:, xj-Npt:xj+Npt], (dim*2*Npt, 1), order='F')
                        UmR = immergedInterfaces.modifiedValuesRbis(
                            U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                        UmL = immergedInterfaces.modifiedValuesLbis(
                            U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                        Ujum=-np.dot(mR,Uint)+np.dot(mL,Uint)
                        Umoy=np.dot(mR,Uint)+np.dot(mL,Uint)
                        CSVM[ts,0+interS*4]=Ujum[0,0]
                        CSVM[ts,1+interS*4]=Ujum[1,0]
                        CSVM[ts,2+interS*4]=1/2*Umoy[0,0]
                        CSVM[ts,3+interS*4]=1/2*Umoy[1,0]
                        tracesL=np.dot(mL,Uint)
                        tracesR=np.dot(mR,Uint)
                        traces[ts,0+interS*4]=tracesL[0,0]
                        traces[ts,1+interS*4]=tracesL[1,0]
                        traces[ts,2+interS*4]=tracesR[0,0]
                        traces[ts,3+interS*4]=tracesR[1,0]
                    if K/2*dx >= alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin < alpha[interS]:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]-rhox[xstep]*cmx[xstep]*dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    elif xstep*dx+xmin >= alpha[interS] and K/2*dx+alpha[interS] > (xstep*dx+xmin):
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]-rhox[xstep]*cmx[xstep]*dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    else:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]-rhox[xstep]*cmx[xstep]*dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        if xstep == xj+int(K/2)+1:
                            interS = interS+1
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    elif rCSVM=="no" and rtraces=="no":
        return U, V[::film], Sig[::film]
    elif rCSVM=="no"and rtraces=="yes":
        return U, V[::film], Sig[::film],traces
    elif rCSVM=="yes"and rtraces=="no":
        return U, V[::film], Sig[::film],CSVM
    else:
        return U, V[::film], Sig[::film],CSVM,traces
    
def FD_sourceVIIKey(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film,rCSVM="no",rtraces="no",printT=10):
    X, xmin, xmax, Nx, dx, ESIM, Npt = inputReading.geometric(config)
    order = 2*ESIM-1
    Nmat, rho, cm = inputReading.material(mat)
    Ninter = int(frontiere["nombre de frontieres "])
    alpha = inputReading.frontiere(frontiere)
    alpha[-1] = xmax+0.0001*xmax
    dim = 2
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    csTot = constructCS(K, Nx, rhox, cmx, dt, dx)
    formalCMat = preTMat(frontiere, rho, cm, order, modulation)
    t = sp.symbols('t')
    FMatrix = np.empty((4, Ninter), dtype=object)
    CSVM=np.zeros([Nt,4*Ninter])
    traces=np.zeros([Nt,4*Ninter])
    for Ni in range(Ninter):
        contacti = "Contact_"+str(Ni)
        if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
            FMatrix[:, Ni] = formalCMat[:, Ni]
        else:
            FMatrix[0, Ni] = sp.lambdify(t, formalCMat[0, Ni], modules="numpy")
            FMatrix[1, Ni] = sp.lambdify(t, formalCMat[1, Ni], modules="numpy")
    matCtot = np.empty((4, Ninter), dtype=object)
    # formalCMat=immergedInterfaces.formatC_SM_modulated(rho,cm,order,formalB,x,t,Lambda,Omega,DF,DM)#,alpha,ts,LambdaM,OmegaM,ampF,ampM)
    for ts in range(Nt):
        if (Nt-ts)//printT==(Nt-ts)/printT:
            print(Nt-ts)
        interS = 0
        # Un = np.zeros([2,Nx])
        for Ni in range(Ninter):
            contacti = "Contact_"+str(Ni)
            if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
                matCtot[:, Ni] = formalCMat[:, Ni]
            else:
                matCtot[0, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[0, Ni], t, -(Nt+1-ts)*dt)
                matCtot[1, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[1, Ni], t, -(Nt+1-ts)*dt)
                matCtot[2, Ni] = immergedInterfaces.matCInter(
                    matCtot[0, Ni], matCtot[1, Ni])
        # print(matCtot)
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            if xstep <= x_0:
                if (xstep+xmin/dx) <= alpha[interS]/dx-K/2:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = Un[1, xstep]
                else:
                    xj = int(alpha[interS]/dx)
                    if xstep == xj-int(K/2)+1:
                        matCL = matCtot[0, interS]
                        matCR = matCtot[1, interS]
                        x = np.linspace(xj-Npt+1, xj+Npt, 2*Npt)*dx
                        mR = immergedInterfaces.matR(
                            dim, order, alpha[interS], x, Npt, matCL, matCR)
                        mL = np.dot(matCtot[2, interS], mR)
                        # print(mL)
                        Uint = np.reshape(
                            U[:, xj-Npt+1:xj+Npt+1], (dim*2*Npt, 1), order='F')
                        UmR = immergedInterfaces.modifiedValuesR(
                            U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                        UmL = immergedInterfaces.modifiedValuesL(
                            U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                        Ujum=-np.dot(mR,Uint)+np.dot(mL,Uint)
                        Umoy=np.dot(mR,Uint)+np.dot(mL,Uint)
                        CSVM[ts,0+interS*4]=Ujum[0,0]
                        CSVM[ts,1+interS*4]=Ujum[1,0]
                        CSVM[ts,2+interS*4]=1/2*Umoy[0,0]
                        CSVM[ts,3+interS*4]=1/2*Umoy[1,0]
                        tracesL=np.dot(mL,Uint)
                        tracesR=np.dot(mR,Uint)
                        traces[ts,0+interS*4]=tracesL[0,0]
                        traces[ts,1+interS*4]=tracesL[1,0]
                        traces[ts,2+interS*4]=tracesR[0,0]
                        traces[ts,3+interS*4]=tracesR[1,0]
                    if K/2*dx > alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin <= alpha[interS]:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]
                    elif xstep*dx+xmin > alpha[interS] and K/2*dx+alpha[interS] >= (xstep*dx+xmin):
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]
                    else:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]
                        if xstep == xj+int(K/2)+1:
                            interS = interS+1
            if xstep > x_0:
                if (xstep+xmin/dx) < alpha[interS]/dx-K/2:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = Un[1, xstep]
                else:
                    xj = int(alpha[interS]/dx)
                    if xstep == xj-int(K/2):
                        matCL = matCtot[0, interS]
                        matCR = matCtot[1, interS]
                        x = np.linspace(xj-Npt, xj+Npt-1, 2*Npt)*dx
                        mR = immergedInterfaces.matR(
                            dim, order, alpha[interS], x, Npt, matCL, matCR)
                        mL = np.dot(matCtot[2, interS], mR)
                        #print(Npt)
                        Uint = np.reshape(
                            U[:, xj-Npt:xj+Npt], (dim*2*Npt, 1), order='F')
                        UmR = immergedInterfaces.modifiedValuesRbis(
                            U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                        UmL = immergedInterfaces.modifiedValuesLbis(
                            U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                        Ujum=-np.dot(mR,Uint)+np.dot(mL,Uint)
                        Umoy=np.dot(mR,Uint)+np.dot(mL,Uint)
                        CSVM[ts,0+interS*4]=Ujum[0,0]
                        CSVM[ts,1+interS*4]=Ujum[1,0]
                        CSVM[ts,2+interS*4]=1/2*Umoy[0,0]
                        CSVM[ts,3+interS*4]=1/2*Umoy[1,0]
                        tracesL=np.dot(mL,Uint)
                        tracesR=np.dot(mR,Uint)
                        traces[ts,0+interS*4]=tracesL[0,0]
                        traces[ts,1+interS*4]=tracesL[1,0]
                        traces[ts,2+interS*4]=tracesR[0,0]
                        traces[ts,3+interS*4]=tracesR[1,0]
                    if K/2*dx >= alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin < alpha[interS]:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]
                    elif xstep*dx+xmin >= alpha[interS] and K/2*dx+alpha[interS] > (xstep*dx+xmin):
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]
                    else:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]
                        if xstep == xj+int(K/2)+1:
                            interS = interS+1
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    elif rCSVM=="no" and rtraces=="no":
        return U, V[::film], Sig[::film]
    elif rCSVM=="no"and rtraces=="yes":
        return U, V[::film], Sig[::film],traces
    elif rCSVM=="yes"and rtraces=="no":
        return U, V[::film], Sig[::film],CSVM
    else:
        return U, V[::film], Sig[::film],CSVM,traces

def FD_sourceVIIKeyOW(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film,rCSVM="no",rtraces="no",printT=10):
    X, xmin, xmax, Nx, dx, ESIM, Npt = inputReading.geometric(config)
    order = 2*ESIM-1
    Nmat, rho, cm = inputReading.material(mat)
    Ninter = int(frontiere["nombre de frontieres "])
    alpha = inputReading.frontiere(frontiere)
    alpha[-1] = xmax+0.0001*xmax
    dim = 2
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    csTot = constructCS(K, Nx, rhox, cmx, dt, dx)
    formalCMat = preTMat(frontiere, rho, cm, order, modulation)
    t = sp.symbols('t')
    FMatrix = np.empty((4, Ninter), dtype=object)
    CSVM=np.zeros([Nt,4*Ninter])
    traces=np.zeros([Nt,4*Ninter])
    for Ni in range(Ninter):
        contacti = "Contact_"+str(Ni)
        if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
            FMatrix[:, Ni] = formalCMat[:, Ni]
        else:
            FMatrix[0, Ni] = sp.lambdify(t, formalCMat[0, Ni], modules="numpy")
            FMatrix[1, Ni] = sp.lambdify(t, formalCMat[1, Ni], modules="numpy")
    matCtot = np.empty((4, Ninter), dtype=object)
    # formalCMat=immergedInterfaces.formatC_SM_modulated(rho,cm,order,formalB,x,t,Lambda,Omega,DF,DM)#,alpha,ts,LambdaM,OmegaM,ampF,ampM)
    for ts in range(Nt):
        if (Nt-ts)//printT==(Nt-ts)/printT:
            print(Nt-ts)
        interS = 0
        # Un = np.zeros([2,Nx])
        for Ni in range(Ninter):
            contacti = "Contact_"+str(Ni)
            if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
                matCtot[:, Ni] = formalCMat[:, Ni]
            else:
                matCtot[0, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[0, Ni], t, -(Nt+1-ts)*dt)
                matCtot[1, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[1, Ni], t, -(Nt+1-ts)*dt)
                matCtot[2, Ni] = immergedInterfaces.matCInter(
                    matCtot[0, Ni], matCtot[1, Ni])
        # print(matCtot)
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            if xstep <= x_0:
                if (xstep+xmin/dx) <= alpha[interS]/dx-K/2:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = Un[1, xstep]+rhox[xstep]*cmx[xstep]*dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                else:
                    xj = int(alpha[interS]/dx)
                    if xstep == xj-int(K/2)+1:
                        matCL = matCtot[0, interS]
                        matCR = matCtot[1, interS]
                        x = np.linspace(xj-Npt+1, xj+Npt, 2*Npt)*dx
                        mR = immergedInterfaces.matR(
                            dim, order, alpha[interS], x, Npt, matCL, matCR)
                        mL = np.dot(matCtot[2, interS], mR)
                        # print(mL)
                        Uint = np.reshape(
                            U[:, xj-Npt+1:xj+Npt+1], (dim*2*Npt, 1), order='F')
                        UmR = immergedInterfaces.modifiedValuesR(
                            U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                        UmL = immergedInterfaces.modifiedValuesL(
                            U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                        Ujum=-np.dot(mR,Uint)+np.dot(mL,Uint)
                        Umoy=np.dot(mR,Uint)+np.dot(mL,Uint)
                        CSVM[ts,0+interS*4]=Ujum[0,0]
                        CSVM[ts,1+interS*4]=Ujum[1,0]
                        CSVM[ts,2+interS*4]=1/2*Umoy[0,0]
                        CSVM[ts,3+interS*4]=1/2*Umoy[1,0]
                        tracesL=np.dot(mL,Uint)
                        tracesR=np.dot(mR,Uint)
                        traces[ts,0+interS*4]=tracesL[0,0]
                        traces[ts,1+interS*4]=tracesL[1,0]
                        traces[ts,2+interS*4]=tracesR[0,0]
                        traces[ts,3+interS*4]=tracesR[1,0]
                    if K/2*dx > alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin <= alpha[interS]:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]+rhox[xstep]*cmx[xstep]*dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    elif xstep*dx+xmin > alpha[interS] and K/2*dx+alpha[interS] >= (xstep*dx+xmin):
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]+rhox[xstep]*cmx[xstep]*dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    else:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]+rhox[xstep]*cmx[xstep]*dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        if xstep == xj+int(K/2)+1:
                            interS = interS+1
            if xstep > x_0:
                if (xstep+xmin/dx) < alpha[interS]/dx-K/2:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = Un[1, xstep]+rhox[xstep]*cmx[xstep]*dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                else:
                    xj = int(alpha[interS]/dx)
                    if xstep == xj-int(K/2):
                        matCL = matCtot[0, interS]
                        matCR = matCtot[1, interS]
                        x = np.linspace(xj-Npt, xj+Npt-1, 2*Npt)*dx
                        mR = immergedInterfaces.matR(
                            dim, order, alpha[interS], x, Npt, matCL, matCR)
                        mL = np.dot(matCtot[2, interS], mR)
                        #print(Npt)
                        Uint = np.reshape(
                            U[:, xj-Npt:xj+Npt], (dim*2*Npt, 1), order='F')
                        UmR = immergedInterfaces.modifiedValuesRbis(
                            U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                        UmL = immergedInterfaces.modifiedValuesLbis(
                            U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                        Ujum=-np.dot(mR,Uint)+np.dot(mL,Uint)
                        Umoy=np.dot(mR,Uint)+np.dot(mL,Uint)
                        CSVM[ts,0+interS*4]=Ujum[0,0]
                        CSVM[ts,1+interS*4]=Ujum[1,0]
                        CSVM[ts,2+interS*4]=1/2*Umoy[0,0]
                        CSVM[ts,3+interS*4]=1/2*Umoy[1,0]
                        tracesL=np.dot(mL,Uint)
                        tracesR=np.dot(mR,Uint)
                        traces[ts,0+interS*4]=tracesL[0,0]
                        traces[ts,1+interS*4]=tracesL[1,0]
                        traces[ts,2+interS*4]=tracesR[0,0]
                        traces[ts,3+interS*4]=tracesR[1,0]
                    if K/2*dx >= alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin < alpha[interS]:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]+rhox[xstep]*cmx[xstep]*dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    elif xstep*dx+xmin >= alpha[interS] and K/2*dx+alpha[interS] > (xstep*dx+xmin):
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]+rhox[xstep]*cmx[xstep]*dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    else:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        Utn[1, xstep] = Un[1, xstep]+rhox[xstep]*cmx[xstep]*dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        if xstep == xj+int(K/2)+1:
                            interS = interS+1
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    elif rCSVM=="no" and rtraces=="no":
        return U, V[::film], Sig[::film]
    elif rCSVM=="no"and rtraces=="yes":
        return U, V[::film], Sig[::film],traces
    elif rCSVM=="yes"and rtraces=="no":
        return U, V[::film], Sig[::film],CSVM
    else:
        return U, V[::film], Sig[::film],CSVM,traces

def FD_sourceVIIComplexSF(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film):
    X, xmin, xmax, Nx, dx, ESIM, Npt = inputReading.geometric(config)
    order = 2*ESIM-1
    Nmat, rho, cm = inputReading.material(mat)
    Ninter = int(frontiere["nombre de frontieres "])
    alpha = inputReading.frontiere(frontiere)
    alpha[-1] = xmax+0.0001*xmax
    dim = 2
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)], dtype='complex_')
    V = np.zeros([int(Nt), int(Nx)], dtype='complex_')
    Sig = np.zeros([int(Nt), int(Nx)], dtype='complex_')
    csTot = constructCS(K, Nx, rhox, cmx, dt, dx)
    formalCMat, Fform, Kform = preTMatK(frontiere, rho, cm, order, modulation)
    t = sp.symbols('t')
    FMatrix = np.empty((4, Ninter), dtype=object)

    for Ni in range(Ninter):
        contacti = "Contact_"+str(Ni)
        if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
            FMatrix[:, Ni] = formalCMat[:, Ni]
        else:
            FMatrix[0, Ni] = sp.lambdify(t, formalCMat[0, Ni], modules="numpy")
            FMatrix[1, Ni] = sp.lambdify(t, formalCMat[1, Ni], modules="numpy")
    matCtot = np.empty((4, Ninter), dtype=object)
    # formalCMat=immergedInterfaces.formatC_SM_modulated(rho,cm,order,formalB,x,t,Lambda,Omega,DF,DM)#,alpha,ts,LambdaM,OmegaM,ampF,ampM)
    for ts in range(Nt):
        print(ts)
        interS = 0
        # Un = np.zeros([2,Nx])
        for Ni in range(Ninter):
            contacti = "Contact_"+str(Ni)
            if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
                matCtot[:, Ni] = formalCMat[:, Ni]
            else:
                matCtot[0, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[0, Ni], t, ts*dt)
                matCtot[1, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[1, Ni], t, ts*dt)
                matCtot[2, Ni] = immergedInterfaces.matCInter(
                    matCtot[0, Ni], matCtot[1, Ni])
        # print(matCtot)
        Un = np.zeros([2, Nx], dtype='complex_')
        for xstep in range(int(K/2), Nx-int(K/2)):
            if xstep <= x_0:
                if (xstep+xmin/dx) <= alpha[interS]/dx-K/2:
                    Us = np.zeros([2, 1], dtype='complex_')
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]
                    Utn[1, xstep] = Un[1, xstep]
                else:
                    xj = int(alpha[interS]/dx)
                    if xstep == xj-int(K/2)+1:
                        matCL = matCtot[0, interS]
                        matCR = matCtot[1, interS]
                        x = np.linspace(xj-Npt+1, xj+Npt, 2*Npt)*dx
                        mR = immergedInterfaces.matR(
                            dim, order, alpha[interS], x, Npt, matCL, matCR)
                        mL = np.dot(matCtot[2, interS], mR)
                        # print(mL)
                        Uint = np.reshape(
                            U[:, xj-Npt+1:xj+Npt+1], (dim*2*Npt, 1), order='F')
                        UmR = immergedInterfaces.modifiedValuesR(
                            U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                        UmL = immergedInterfaces.modifiedValuesL(
                            U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                    if K/2*dx > alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin <= alpha[interS]:
                        Us = np.zeros([2, 1], dtype='complex_')
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]
                        Utn[1, xstep] = Un[1, xstep]
                    elif xstep*dx+xmin > alpha[interS] and K/2*dx+alpha[interS] >= (xstep*dx+xmin):
                        Us = np.zeros([2, 1], dtype='complex_')
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]
                        Utn[1, xstep] = Un[1, xstep]
                    else:
                        Us = np.zeros([2, 1], dtype='complex_')
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]
                        Utn[1, xstep] = Un[1, xstep]
                        if xstep == xj+int(K/2)+1:
                            interS = interS+1
            if xstep > x_0:
                if (xstep+xmin/dx) < alpha[interS]/dx-K/2:
                    Us = np.zeros([2, 1], dtype='complex_')
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]
                    Utn[1, xstep] = Un[1, xstep]
                else:
                    xj = int(alpha[interS]/dx)
                    if xstep == xj-int(K/2):
                        matCL = matCtot[0, interS]
                        matCR = matCtot[1, interS]
                        x = np.linspace(xj-Npt, xj+Npt-1, 2*Npt)*dx
                        mR = immergedInterfaces.matR(
                            dim, order, alpha[interS], x, Npt, matCL, matCR)
                        mL = np.dot(matCtot[2, interS], mR)
                        # print(mL)
                        Uint = np.reshape(
                            U[:, xj-Npt:xj+Npt], (dim*2*Npt, 1), order='F')
                        UmR = immergedInterfaces.modifiedValuesRbis(
                            U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                        UmL = immergedInterfaces.modifiedValuesLbis(
                            U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                    if K/2*dx >= alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin < alpha[interS]:
                        Us = np.zeros([2, 1], dtype='complex_')
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]
                        Utn[1, xstep] = Un[1, xstep]
                    elif xstep*dx+xmin >= alpha[interS] and K/2*dx+alpha[interS] > (xstep*dx+xmin):
                        Us = np.zeros([2, 1], dtype='complex_')
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]
                        Utn[1, xstep] = Un[1, xstep]
                    else:
                        Us = np.zeros([2, 1], dtype='complex_')
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]
                        Utn[1, xstep] = Un[1, xstep]
                        if xstep == xj+int(K/2)+1:
                            interS = interS+1
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]


def FD_sourceSII(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film):
    X, xmin, xmax, Nx, dx, ESIM, Npt = inputReading.geometric(config)
    order = 2*ESIM-1
    Nmat, rho, cm = inputReading.material(mat)
    Ninter = int(frontiere["nombre de frontieres "])
    alpha = inputReading.frontiere(frontiere)
    alpha[-1] = xmax+0.0001*xmax
    dim = 2
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    csTot = constructCS(K, Nx, rhox, cmx, dt, dx)
    formalCMat = preTMat(frontiere, rho, cm, order, modulation)
    t = sp.symbols('t')
    FMatrix = np.empty((4, Ninter), dtype=object)
    for Ni in range(Ninter):
        contacti = "Contact_"+str(Ni)
        if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
            FMatrix[:, Ni] = formalCMat[:, Ni]
        else:
            FMatrix[0, Ni] = sp.lambdify(t, formalCMat[0, Ni], modules="numpy")
            FMatrix[1, Ni] = sp.lambdify(t, formalCMat[1, Ni], modules="numpy")
    matCtot = np.empty((4, Ninter), dtype=object)
    # formalCMat=immergedInterfaces.formatC_SM_modulated(rho,cm,order,formalB,x,t,Lambda,Omega,DF,DM)#,alpha,ts,LambdaM,OmegaM,ampF,ampM)
    for ts in range(Nt):
        print(ts)
        interS = 0
        # Un = np.zeros([2,Nx])
        for Ni in range(Ninter):
            contacti = "Contact_"+str(Ni)
            if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
                matCtot[:, Ni] = formalCMat[:, Ni]
            else:
                matCtot[0, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[0, Ni], t, ts*dt)
                matCtot[1, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[1, Ni], t, ts*dt)
                matCtot[2, Ni] = immergedInterfaces.matCInter(
                    matCtot[0, Ni], matCtot[1, Ni])
        # print(matCtot)
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            if xstep <= x_0:
                if (xstep+xmin/dx) <= alpha[interS]/dx-K/2:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]
                    Utn[1, xstep] = Un[1, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                else:
                    xj = int(alpha[interS]/dx)
                    if xstep == xj-int(K/2)+1:
                        matCL = matCtot[0, interS]
                        matCR = matCtot[1, interS]
                        x = np.linspace(xj-Npt+1, xj+Npt, 2*Npt)*dx
                        mR = immergedInterfaces.matR(
                            dim, order, alpha[interS], x, Npt, matCL, matCR)
                        mL = np.dot(matCtot[2, interS], mR)
                        # print(mL)
                        Uint = np.reshape(
                            U[:, xj-Npt+1:xj+Npt+1], (dim*2*Npt, 1), order='F')
                        UmR = immergedInterfaces.modifiedValuesR(
                            U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                        UmL = immergedInterfaces.modifiedValuesL(
                            U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                    if K/2*dx > alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin <= alpha[interS]:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]
                        Utn[1, xstep] = Un[1, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    elif xstep*dx+xmin > alpha[interS] and K/2*dx+alpha[interS] >= (xstep*dx+xmin):
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]
                        Utn[1, xstep] = Un[1, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    else:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]
                        Utn[1, xstep] = Un[1, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        if xstep == xj+int(K/2)+1:
                            interS = interS+1
            if xstep > x_0:
                if (xstep+xmin/dx) < alpha[interS]/dx-K/2:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]
                    Utn[1, xstep] = Un[1, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                else:
                    xj = int(alpha[interS]/dx)
                    if xstep == xj-int(K/2):
                        matCL = matCtot[0, interS]
                        matCR = matCtot[1, interS]
                        x = np.linspace(xj-Npt, xj+Npt-1, 2*Npt)*dx
                        mR = immergedInterfaces.matR(
                            dim, order, alpha[interS], x, Npt, matCL, matCR)
                        mL = np.dot(matCtot[2, interS], mR)
                        # print(mL)
                        Uint = np.reshape(
                            U[:, xj-Npt:xj+Npt], (dim*2*Npt, 1), order='F')
                        UmR = immergedInterfaces.modifiedValuesRbis(
                            U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                        UmL = immergedInterfaces.modifiedValuesLbis(
                            U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                    if K/2*dx >= alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin < alpha[interS]:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]
                        Utn[1, xstep] = Un[1, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    elif xstep*dx+xmin >= alpha[interS] and K/2*dx+alpha[interS] > (xstep*dx+xmin):
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]
                        Utn[1, xstep] = Un[1, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    else:
                        Us = np.zeros([2, 1])
                        # rhos=rhox[xstep]
                        # cms=cmx[xstep]
                        # A=inputReading.comA(rhos,cms)
                        for s in range(K+1):
                            cs = csTot[s, xstep]
                            # cs=CS(K,s,A,dt,dx)
                            sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                            Us[0, 0] = Us[0, 0]+sumUs[0]
                            Us[1, 0] = Us[1, 0]+sumUs[1]
                        Un[0, xstep] = U[0, xstep]-Us[0, 0]
                        Un[1, xstep] = U[1, xstep]-Us[1, 0]
                        Utn[0, xstep] = Un[0, xstep]
                        Utn[1, xstep] = Un[1, xstep]+dt * \
                            dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                        if xstep == xj+int(K/2)+1:
                            interS = interS+1
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]


def FD_sourceVII_Homo(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film):
    X, xmin, xmax, Nx, dx, ESIM, Npt = inputReading.geometric(config)
    order = 2*ESIM-1
    Nmat, rho, cm = inputReading.material(mat)
    Ninter = int(frontiere["nombre de frontieres "])
    alpha = inputReading.frontiere(frontiere)
    alpha[-1] = xmax+0.0001*xmax
    dim = 2
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    csTot = constructCS(K, Nx, rhox, cmx, dt, dx)
    formalCMat = preTMat(frontiere, rho, cm, order, modulation)
    t = sp.symbols('t')
    FMatrix = np.empty((4, Ninter), dtype=object)
    for Ni in range(Ninter):
        contacti = "Contact_"+str(Ni)
        if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
            FMatrix[:, Ni] = formalCMat[:, Ni]
        else:
            FMatrix[0, Ni] = sp.lambdify(t, formalCMat[0, Ni], modules="numpy")
            FMatrix[1, Ni] = sp.lambdify(t, formalCMat[1, Ni], modules="numpy")
    matCtot = np.empty((4, Ninter), dtype=object)
    # formalCMat=immergedInterfaces.formatC_SM_modulated(rho,cm,order,formalB,x,t,Lambda,Omega,DF,DM)#,alpha,ts,LambdaM,OmegaM,ampF,ampM)
    for ts in range(Nt):
        print(ts)
        interS = 0
        # Un = np.zeros([2,Nx])
        for Ni in range(Ninter):
            contacti = "Contact_"+str(Ni)
            if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
                matCtot[:, Ni] = formalCMat[:, Ni]
            else:
                matCtot[0, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[0, Ni], t, ts*dt)
                matCtot[1, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[1, Ni], t, ts*dt)
                matCtot[2, Ni] = immergedInterfaces.matCInter(
                    matCtot[0, Ni], matCtot[1, Ni])
        # print(matCtot)
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            if (xstep+xmin/dx) <= alpha[interS]/dx-K/2:
                Us = np.zeros([2, 1])
                # rhos=rhox[xstep]
                # cms=cmx[xstep]
                # A=inputReading.comA(rhos,cms)
                for s in range(K+1):
                    cs = csTot[s, xstep]
                    # cs=CS(K,s,A,dt,dx)
                    sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                    Us[0, 0] = Us[0, 0]+sumUs[0]
                    Us[1, 0] = Us[1, 0]+sumUs[1]
                Un[0, xstep] = U[0, xstep]-Us[0, 0]
                Un[1, xstep] = U[1, xstep]-Us[1, 0]
                Utn[0, xstep] = Un[0, xstep]+dt * \
                    dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                Utn[1, xstep] = Un[1, xstep]
            else:
                xj = int(alpha[interS]/dx)
                if xstep == xj-int(K/2)+1:
                    matCL = matCtot[0, interS]
                    matCR = matCtot[1, interS]
                    x = np.linspace(xj-Npt+1, xj+Npt, 2*Npt)*dx
                    mR = immergedInterfaces.matR(
                        dim, order, alpha[interS], x, Npt, matCL, matCR)
                    mL = np.dot(matCtot[2, interS], mR)
                    # print(mL)
                    Uint = np.reshape(U[:, xj-Npt+1:xj+Npt+1],
                                      (dim*2*Npt, 1), order='F')
                    UmR = immergedInterfaces.modifiedValuesR(
                        U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                    UmL = immergedInterfaces.modifiedValuesL(
                        U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                if K/2*dx > alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin <= alpha[interS]:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = Un[1, xstep]
                elif xstep*dx+xmin > alpha[interS] and K/2*dx+alpha[interS] >= (xstep*dx+xmin):
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = Un[1, xstep]
                else:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = Un[1, xstep]
                    if xstep == xj+int(K/2)+1:
                        interS = interS+1
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]


def FD_sourceVII_opti(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, modulation, film):
    X, xmin, xmax, Nx, dx, ESIM, Npt = inputReading.geometric(config)
    order = 2*ESIM-1
    Nmat, rho, cm = inputReading.material(mat)
    Ninter = int(frontiere["nombre de frontieres "])
    alpha = inputReading.frontiere(frontiere)
    alpha[-1] = xmax+0.01*xmax
    dim = 2
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    csTot = constructCS(K, Nx, rhox, cmx, dt, dx)
    formalCMat = preTMat(frontiere, rho, cm, order, modulation, ts=0, dt=1)
    t = sp.symbols('t')
    FMatrix = np.empty((4, Ninter), dtype=object)
    for Ni in range(Ninter):
        contacti = "Contact_"+str(Ni)
        if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
            FMatrix[:, Ni] = formalCMat[:, Ni]
        else:
            FMatrix[0, Ni] = sp.lambdify(t, formalCMat[0, Ni], modules="numpy")
            FMatrix[1, Ni] = sp.lambdify(t, formalCMat[1, Ni], modules="numpy")
    matCtot = np.empty((4, Ninter), dtype=object)
    # formalCMat=immergedInterfaces.formatC_SM_modulated(rho,cm,order,formalB,x,t,Lambda,Omega,DF,DM)#,alpha,ts,LambdaM,OmegaM,ampF,ampM)
    for ts in range(Nt):
        print(ts)
        interS = 0
        # Un = np.zeros([2,Nx])
        for Ni in range(Ninter):
            contacti = "Contact_"+str(Ni)
            if frontiere[contacti] == "parfait" or frontiere[contacti] == "masse-ressort":
                matCtot[:, Ni] = formalCMat[:, Ni]
            else:
                matCtot[0, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[0, Ni], t, ts*dt)
                matCtot[1, Ni] = immergedInterfaces.valueModMatrix(
                    FMatrix[1, Ni], t, ts*dt)
                matCtot[2, Ni] = immergedInterfaces.matCInter(
                    matCtot[0, Ni], matCtot[1, Ni])
        for xstep in range(int(K/2), Nx-int(K/2)):
            if (xstep*dx+xmin) <= alpha[interS]-K/2*dx:
                # time0=time.time()
                Us = np.zeros([2, 1])
                # rhos=rhox[xstep]
                # cms=cmx[xstep]
                # A=inputReading.comA(rhos,cms)
                time0 = time.time()
                s_values = np.arange(K + 1)
                cst = np.stack(csTot[s_values, xstep])
                cst = cst.reshape(K+1, 2, 2)
                Ured = U[:, xstep-int(K/2)+s_values]
                print("Forme de cst:", cst.shape)
                print("Forme de Ured:", Ured.shape)
                # Ured=Ured.reshape(-1, 2)
                # if cst.shape[2] != Ured.shape[0]:
                #     raise ValueError("Les dimensions ne correspondent pas.")
                # .sum(axis=0)#np.tensordot(cst, Ured, axes=([0], [0]))# np.einsum('ijk,ki->ij', cst, Ured).sum(axis=0)
                Us_sum = np.matmul(cst.transpose(1, 2, 0), Ured)
                print(time.time()-time0)
                # np.einsum('ijk,ik->ij', cst, Ured).sum(axis=0, keepdims=True).T
                # #Us = np.dot(Ured,cst)
                time0 = time.time()
                Us = np.zeros([2, 1])
                for s in range(K+1):
                    cs = csTot[s, xstep]
                    # cs=CS(K,s,A,dt,dx)
                    sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                    # print(sumUs)
                    Us[0, 0] = Us[0, 0]+sumUs[0]
                    Us[1, 0] = Us[1, 0]+sumUs[1]
                # print(Us)
                print(time.time()-time0, 'time')
                # Un[0,xstep]=U[0,xstep]-Us[0,0]
                # Un[1,xstep]=U[1,xstep]-Us[1,0]
                Utn[0, xstep] = U[0, xstep]-Us_sum[0]+dt * \
                    dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                Utn[1, xstep] = U[1, xstep]-Us_sum[1]
                # timef=time.time()-time0
                # print(timef)
            else:
                xj = int(alpha[interS]/dx)
                if xstep == xj-int(K/2)+1:
                    matCL = matCtot[0, interS]
                    matCR = matCtot[1, interS]
                    x = np.linspace(xj-Npt+1, xj+Npt, 2*Npt)*dx
                    mR = immergedInterfaces.matR(
                        dim, order, alpha[interS], x, Npt, matCL, matCR)
                    mL = np.dot(matCtot[2, interS], mR)
                    Uint = np.reshape(U[:, xj-Npt+1:xj+Npt+1],
                                      (dim*2*Npt, 1), order='F')
                    UmR = immergedInterfaces.modifiedValuesR(
                        U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                    UmL = immergedInterfaces.modifiedValuesL(
                        U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                if K/2*dx > alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin <= alpha[interS]:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    s_values = np.arange(K + 1)
                    cst = np.stack(csTot[s_values, xstep])
                    cst = cst.reshape(K+1, 2, 2)
                    Ured = UmR[:, xstep-int(K/2)+s_values]
                    # Ured=Ured.reshape(-1, 2)
                    # if cst.shape[2] != Ured.shape[0]:
                    #     raise ValueError("Les dimensions ne correspondent pas.")
                    # np.einsum('ijk,ik->ij', cst, Ured).sum(axis=0, keepdims=True).T
                    Us = np.einsum('ijk,ki->ij', cst, Ured)
                    Us_sum = np.sum(Us, axis=0)
                    # for s in range(K+1):
                    #     cs=csTot[s,xstep]
                    #     #cs=CS(K,s,A,dt,dx)
                    #     sumUs=np.dot(cs,UmR[:,xstep+s-int(K/2)])
                    #     Us[0,0]=Us[0,0]+sumUs[0]
                    #     Us[1,0]=Us[1,0]+sumUs[1]
                    # Un[0,xstep]=U[0,xstep]-Us[0,0]
                    # Un[1,xstep]=U[1,xstep]-Us[1,0]
                    Utn[0, xstep] = U[0, xstep]-Us_sum[0]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = U[1, xstep]-Us_sum[1]
                elif xstep*dx+xmin > alpha[interS] and K/2*dx+alpha[interS] >= (xstep*dx+xmin):
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    s_values = np.arange(K + 1)
                    cst = np.stack(csTot[s_values, xstep])
                    cst = cst.reshape(K+1, 2, 2)
                    Ured = UmL[:, xstep-int(K/2)+s_values]
                    # Ured=Ured.reshape(-1, 2)
                    # if cst.shape[2] != Ured.shape[0]:
                    #     raise ValueError("Les dimensions ne correspondent pas.")
                    # np.einsum('ijk,ik->ij', cst, Ured).sum(axis=0, keepdims=True).T
                    Us = np.einsum('ijk,ki->ij', cst, Ured)
                    Us_sum = np.sum(Us, axis=0)
                    # for s in range(K+1):
                    #     cs=csTot[s,xstep]
                    #     #cs=CS(K,s,A,dt,dx)
                    #     sumUs=np.dot(cs,UmL[:,xstep+s-int(K/2)])
                    #     Us[0,0]=Us[0,0]+sumUs[0]
                    #     Us[1,0]=Us[1,0]+sumUs[1]
                    # Un[0,xstep]=U[0,xstep]-Us[0,0]
                    # Un[1,xstep]=U[1,xstep]-Us[1,0]
                    Utn[0, xstep] = U[0, xstep]-Us_sum[0]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = U[1, xstep]-Us_sum[1]
                else:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    s_values = np.arange(K + 1)
                    cst = np.stack(csTot[s_values, xstep])
                    cst = cst.reshape(K+1, 2, 2)
                    Ured = U[:, xstep-int(K/2)+s_values]
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
                    Utn[0, xstep] = U[0, xstep]-Us_sum[0]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = U[1, xstep]-Us_sum[1]
                    if xstep == xj+int(K/2)+1:
                        interS = interS+1
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]


def FD_sourceVH(U0, Nt, Nx, K, rhox, cmx, dt, dx, x_0, sce_t, config, mat, frontiere, film):
    X, xmin, xmax, Nx, dx, ESIM, Npt = inputReading.geometric(config)
    order = 2*ESIM-1
    Nmat, rho, cm = inputReading.material(mat)
    alpha = inputReading.frontiere(frontiere)
    alpha[-1] = xmax+0.0001*xmax
    matCtot = preTMat(frontiere, rho, cm, order)
    dim = 2
    U = U0.copy()
    Utn = np.zeros([2, int(Nx)])
    V = np.zeros([int(Nt), int(Nx)])
    Sig = np.zeros([int(Nt), int(Nx)])
    csTot = constructCS(K, Nx, rhox, cmx, dt, dx)
    for ts in range(Nt):
        print(ts)
        interS = 0
        Un = np.zeros([2, Nx])
        for xstep in range(int(K/2), Nx-int(K/2)):
            if (xstep*dx+xmin) <= alpha[interS]-K/2*dx:
                Us = np.zeros([2, 1])
                # rhos=rhox[xstep]
                # cms=cmx[xstep]
                # A=inputReading.comA(rhos,cms)
                for s in range(K+1):
                    cs = csTot[s, xstep]
                    # cs=CS(K,s,A,dt,dx)
                    sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                    Us[0, 0] = Us[0, 0]+sumUs[0]
                    Us[1, 0] = Us[1, 0]+sumUs[1]
                Un[0, xstep] = U[0, xstep]-Us[0, 0]
                Un[1, xstep] = U[1, xstep]-Us[1, 0]
                Utn[0, xstep] = Un[0, xstep]+dt * \
                    dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                Utn[1, xstep] = Un[1, xstep]
            else:
                xj = int(alpha[interS]/dx)
                if xstep == xj-int(K/2)+1:
                    matCL = matCtot[0, interS]
                    matCR = matCtot[1, interS]
                    x = np.linspace(xj-Npt+1, xj+Npt, 2*Npt)*dx
                    mR = immergedInterfaces.matR(
                        dim, order, alpha[interS], x, Npt, matCL, matCR)
                    mL = np.dot(matCtot[2, interS], mR)
                    Uint = np.reshape(U[:, xj-Npt+1:xj+Npt+1],
                                      (dim*2*Npt, 1), order='F')
                    UmR = immergedInterfaces.modifiedValuesR(
                        U, Uint, mR, xj, Npt, dim, order, alpha[interS], x)
                    UmL = immergedInterfaces.modifiedValuesL(
                        U, Uint, mL, xj, Npt, dim, order, alpha[interS], x)
                if K/2*dx > alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin <= alpha[interS]:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, UmR[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = Un[1, xstep]
                elif xstep*dx+xmin > alpha[interS] and K/2*dx+alpha[interS] >= (xstep*dx+xmin):
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, UmL[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = Un[1, xstep]

                else:
                    Us = np.zeros([2, 1])
                    # rhos=rhox[xstep]
                    # cms=cmx[xstep]
                    # A=inputReading.comA(rhos,cms)
                    for s in range(K+1):
                        cs = csTot[s, xstep]
                        # cs=CS(K,s,A,dt,dx)
                        sumUs = np.dot(cs, U[:, xstep+s-int(K/2)])
                        Us[0, 0] = Us[0, 0]+sumUs[0]
                        Us[1, 0] = Us[1, 0]+sumUs[1]
                    Un[0, xstep] = U[0, xstep]-Us[0, 0]
                    Un[1, xstep] = U[1, xstep]-Us[1, 0]
                    Utn[0, xstep] = Un[0, xstep]+dt * \
                        dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
                    Utn[1, xstep] = Un[1, xstep]
                    if xstep == xj+int(K/2)+1:
                        interS = interS+1
        U = Utn[0:2, :].copy()
        V[ts, :] = U[0, :]
        Sig[ts, :] = U[1, :]
    if film == "no":
        return U
    else:
        return U, V[::film], Sig[::film]


# def FD_sourceSII(U0,Nt,Nx,K,rhox,cmx,dt,dx,x_0,sce_t,config,mat,frontiere,modulation,film):
#     X,xmin,xmax,Nx,dx,ESIM,Npt=inputReading.geometric(config)
#     order=2*ESIM-1
#     Nmat,rho,cm=inputReading.material(mat)
#     alpha=inputReading.frontiere(frontiere)
#     alpha[-1]=xmax+0.0001*xmax
#     matCtot=preTMat(frontiere,rho,cm,order)
#     dim=2
#     U=U0.copy()
#     Utn=np.zeros([2,int(Nx)])
#     V=np.zeros([int(Nt),int(Nx)])
#     Sig=np.zeros([int(Nt),int(Nx)])
#     csTot=constructCS(K,Nx,rhox,cmx,dt,dx)
#     for ts in range(Nt):
#         print(ts)
#         interS=0
#         Un = np.zeros([2,Nx])
#         matCtot=preTMat(frontiere,rho,cm,order,modulation,ts*dt,dt)
#         for xstep in range(int(K/2),Nx-int(K/2)):
#             if (xstep*dx+xmin)<=alpha[interS]-K/2*dx:
#                 Us=np.zeros([2,1])
#                 # rhos=rhox[xstep]
#                 # cms=cmx[xstep]
#                 #A=inputReading.comA(rhos,cms)
#                 for s in range(K+1):
#                     cs=csTot[s,xstep]
#                     #cs=CS(K,s,A,dt,dx)
#                     sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
#                     Us[0,0]=Us[0,0]+sumUs[0]
#                     Us[1,0]=Us[1,0]+sumUs[1]
#                 Un[0,xstep]=U[0,xstep]-Us[0,0]
#                 Un[1,xstep]=U[1,xstep]-Us[1,0]
#                 Utn[0,xstep]=Un[0,xstep]
#                 Utn[1,xstep]=Un[1,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
#             else:
#                 xj=int(alpha[interS]/dx)
#                 if xstep==xj-int(K/2)+1:
#                     matCL=matCtot[0,interS]
#                     matCR=matCtot[1,interS]
#                     x=np.linspace(xj-Npt+1, xj+Npt,2*Npt)*dx
#                     mR=immergedInterfaces.matR(dim,order,alpha[interS],x,Npt,matCL,matCR)
#                     mL=np.dot(matCtot[2,interS],mR)
#                     Uint=np.reshape(U[:,xj-Npt+1:xj+Npt+1],(dim*2*Npt,1),order='F')
#                     UmR=immergedInterfaces.modifiedValuesR(U,Uint,mR,xj,Npt,dim,order,alpha[interS],x)
#                     UmL=immergedInterfaces.modifiedValuesL(U,Uint,mL,xj,Npt,dim,order,alpha[interS],x)
#                 if K/2*dx>alpha[interS]-(xstep*dx+xmin) and xstep*dx+xmin<=alpha[interS]:
#                     Us=np.zeros([2,1])
#                     # rhos=rhox[xstep]
#                     # cms=cmx[xstep]
#                     #A=inputReading.comA(rhos,cms)
#                     for s in range(K+1):
#                         cs=csTot[s,xstep]
#                         #cs=CS(K,s,A,dt,dx)
#                         sumUs=np.dot(cs,UmR[:,xstep+s-int(K/2)])
#                         Us[0,0]=Us[0,0]+sumUs[0]
#                         Us[1,0]=Us[1,0]+sumUs[1]
#                     Un[0,xstep]=U[0,xstep]-Us[0,0]
#                     Un[1,xstep]=U[1,xstep]-Us[1,0]
#                     Utn[0,xstep]=Un[0,xstep]
#                     Utn[1,xstep]=Un[1,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
#                 elif xstep*dx+xmin>alpha[interS] and K/2*dx+alpha[interS]>=(xstep*dx+xmin):
#                     Us=np.zeros([2,1])
#                     #rhos=rhox[xstep]
#                     #cms=cmx[xstep]
#                     #A=inputReading.comA(rhos,cms)
#                     for s in range(K+1):
#                         cs=csTot[s,xstep]
#                         #cs=CS(K,s,A,dt,dx)
#                         sumUs=np.dot(cs,UmL[:,xstep+s-int(K/2)])
#                         Us[0,0]=Us[0,0]+sumUs[0]
#                         Us[1,0]=Us[1,0]+sumUs[1]
#                     Un[0,xstep]=U[0,xstep]-Us[0,0]
#                     Un[1,xstep]=U[1,xstep]-Us[1,0]
#                     Utn[0,xstep]=Un[0,xstep]
#                     Utn[1,xstep]=Un[1,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]

#                 else:
#                     Us=np.zeros([2,1])
#                     # rhos=rhox[xstep]
#                     # cms=cmx[xstep]
#                     #A=inputReading.comA(rhos,cms)
#                     for s in range(K+1):
#                         cs=csTot[s,xstep]
#                         #cs=CS(K,s,A,dt,dx)
#                         sumUs=np.dot(cs,U[:,xstep+s-int(K/2)])
#                         Us[0,0]=Us[0,0]+sumUs[0]
#                         Us[1,0]=Us[1,0]+sumUs[1]
#                     Un[0,xstep]=U[0,xstep]-Us[0,0]
#                     Un[1,xstep]=U[1,xstep]-Us[1,0]
#                     Utn[0,xstep]=Un[0,xstep]
#                     Utn[1,xstep]=Un[1,xstep]+dt*dirac_source(xstep*dx, x_0, dx)*sce_t[ts+1]
#                     if xstep == xj+int(K/2)+1:
#                         interS=interS+1
#         U=Utn[0:2,:].copy()
#         V[ts,:]=U[0,:]
#         Sig[ts,:]=U[1,:]
#     if film=="no":
#         return U
#     else:
#         return U,V[::film],Sig[::film]
