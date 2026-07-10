#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 15:46:50 2024

@author: michael
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
# Some code
# from matplotlib.path import Path
# import matplotlib.patches as patches

def customark():
    from svg_pltmarker import get_marker_from_svg
    return get_marker_from_svg(filepath='../WAPOD/Marker.svg')
    
def format_fig(x,y,title=None,legend=None,sizef=16):
    plt.xlabel(x,fontsize=sizef)
    plt.ylabel(y,fontsize=sizef)
    if title!=None:
        plt.title(title,fontsize=sizef)
    if title!=None:
        plt.legend(fontsize=sizef)
    plt.tight_layout()
    return

def sismo(x,t,data1,numLines,pFilm=1,Nt=None):
    if Nt==None:
        Nt=np.size(data1[0,:])
    maxDatas=np.max([np.max(data1[:Nt:int(Nt/numLines),:]),-np.min(data1[:Nt:int(Nt/numLines),:])])
    datas1=data1[:Nt,:]
    plt.figure()
    for j in range(numLines):
        plt.plot(x,(datas1[j*int(Nt/numLines),:]+2*j*1.02*maxDatas)*pFilm*int(Nt/numLines),'#1f77b4',linewidth=2)
    plt.yticks([i * 2 * 1.02 * maxDatas * pFilm * int(Nt/numLines) for i in range(0,numLines,2)],np.round(t[:Nt:2*pFilm*int(Nt/numLines)], 3))
    format_fig(r"$X$ (m)",r"$T$ (s)")
    
def sismoNorm(x,t,data1,numLines,pFilm=1,Nt=None,norm=1,first="yes"):
    if Nt==None:
        Nt=np.size(data1[0,:])
    maxDatas=np.max([np.max(data1[:Nt:int(Nt/numLines),:]),-np.min(data1[:Nt:int(Nt/numLines),:])])/norm
    datas1=data1[:Nt,:]
    Xmax=x[-1]
    plt.figure()
    if first=="yes":
        for j in range(numLines):
            plt.plot(x-Xmax/2,(datas1[j*int(Nt/numLines),:]/norm+2*j*1.02)*pFilm*int(Nt/numLines),'#1f77b4',linewidth=2)
        plt.yticks([i * 2 * 1.02 * pFilm * int(Nt/numLines) for i in range(0,numLines,2)],np.round(t[:Nt:2*pFilm*int(Nt/numLines)], 3))
    else:
        for j in range(1,numLines):
            plt.plot(x-Xmax/2,(datas1[j*int(Nt/numLines),:]/norm+2*j*1.02)*pFilm*int(Nt/numLines),'#1f77b4',linewidth=2)
        plt.yticks([i * 2 * 1.02 * pFilm * int(Nt/numLines) for i in range(0,numLines,2)],np.round(t[:Nt:2*pFilm*int(Nt/numLines)], 3))
    format_fig(r"$X$ (m)",r"$T$ (s)")
    
def sismo2(x,t,data1,data2,numLines,pFilm=1,Nt=None):
    if Nt==None:
        Nt=np.size(data1[0,:])
    maxDatas=np.max([np.max(data1[:Nt:int(Nt/numLines),:]),-np.min(data1[:Nt:int(Nt/numLines),:])])
    datas1=data1[:Nt,:]
    datas2=data2[:Nt,:]
    plt.figure()
    for j in range(numLines):
        plt.plot(x,(datas2[j*int(Nt/numLines),:]+2*j*1.02*maxDatas)*pFilm*int(Nt/numLines),'#1f77b4',linewidth=2)
        plt.plot(x,(datas1[j*int(Nt/numLines),:]+2*j*1.02*maxDatas)*pFilm*int(Nt/numLines),':r',linewidth=3)
    plt.yticks([i * 2 * 1.02 * maxDatas * pFilm * int(Nt/numLines) for i in range(0,numLines,2)],np.round(t[:Nt:2*pFilm*int(Nt/numLines)], 3))
    format_fig(r"$X$ (m)",r"$T$ (s)")
    
# def create_marker(radius=0.05):
#     # Creer un Path pour un cercle
#     circle_points = 100  # Nombre de points pour approximer le cercle
#     theta = np.linspace(0, 2 * np.pi, circle_points)
    
#     # Coordonnees du cercle
#     x = radius * np.cos(theta)
#     y = radius * np.sin(theta)
    
#     # Creer un Path à partir des points du cercle
#     verts = list(zip(x, y))
#     codes = [Path.MOVETO] + [Path.LINETO] * (len(verts) - 1)
    
#     # Retourner un Path representant le cercle
#     return Path(verts, codes)

def cosine_taper(signal, alpha=0.05,zp=20):
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

def specfkR(U,dx,dt):
    Nt=np.size(U[:,0])
    Nx=np.size(U[0,:])
    # Time padding
    zt=2**13
    zx=2**13
    Uzp=np.zeros([zt,zx])
    Uzp[:Nt,int(zx/2-Nx/2):int(zx/2+Nx/2)]=U
    #Fourier in time domain 
    FFTx = np.fft.fft(Uzp, axis=1)
    #Fourier in space
    FFTt = np.fft.rfft(FFTx, axis=0)
    FFT = np.flip(np.fft.fftshift(FFTt, axes=1), axis=1)
    FFK = np.absolute(FFT)

    # Get the frequency and K vectors
    frqv = np.fft.rfftfreq(zt, dt)
    wavv = np.fft.fftfreq(zx, dx)
    # Centering
    wavv = np.fft.fftshift(wavv)
    return frqv,wavv,FFK

def specfk(U,dx,dt,st=4,sx=4):
    Nt=np.size(U[::st,0])
    Nx=np.size(U[0,::sx])
    """
    FK spectrum of traces using the numpy.fft functions.
    """
    # Time padding
    zpt=2**13
    zps=2**13
    Uzp=np.zeros([zpt,zps])
    Uzp[:Nt,int(zps/2-Nx/2):int(zps/2+Nx/2)]=U[::st,::sx]
    #Fourier in time domain 
    FFTt = np.fft.rfft(Uzp[::], axis=0)
    #Fourier in space
    FFTx = np.fft.fft(FFTt, axis=1)
    FFT = np.flip(np.fft.fftshift(FFTx, axes=1), axis=1)
    FFK = np.absolute(FFT)

    # Get the frequency and K vectors
    frqv = np.fft.rfftfreq(zpt, st*dt)
    wavv = np.fft.fftfreq(zps, sx*dx)
    # Centering
    wavv = np.fft.fftshift(wavv)
    return frqv,wavv,FFK

def specfkfreqneg(U,dx,dt,st=4,sx=4):
    Nt=np.size(U[::st,0])
    Nx=np.size(U[0,::sx])
    """
    FK spectrum of traces using the numpy.fft functions.
    """
    # Time padding
    zpt=2**14
    zps=2**14
    Uzp=np.zeros([zpt,zps])
    Uzp[:Nt,int(zps/2-Nx/2):int(zps/2+Nx/2)]=U[::st,::sx]
    #Fourier in time domain 
    FFTt = np.fft.fft(Uzp[::], axis=0)
    #Fourier in space
    FFTx = np.fft.fft(np.fft.fftshift(FFTt, axes=0), axis=1)
    #FFT = np.flip(np.fft.fftshift(FFTx, axes=1), axis=1)
    FFK = np.absolute(FFTx)

    # Get the frequency and K vectors
    frqv = np.fft.fftfreq(zpt, st*dt)
    wavv = np.fft.fftfreq(zps, sx*dx)
    # Centering
    wavv = np.fft.fftshift(wavv)
    frqv = np.fft.fftshift(frqv)
    return frqv,wavv,FFK

def specfkfreqnegB(U,dx,dt,st=8,sx=2):
    Nt=np.size(U[::st,0])
    Nx=np.size(U[0,::sx])
    """
    FK spectrum of traces using the numpy.fft functions.
    """
    # Time padding
    zpt=2**15
    zps=2**12
    Uzp=np.zeros([zpt,zps])
    Uzp[:Nt,int(zps/2-Nx/2):int(zps/2+Nx/2)]=U[::st,::sx]
    #Fourier in time domain 
    FFTt = np.fft.fft(Uzp[::], axis=0)
    #Fourier in space
    FFTx = np.fft.fft(np.fft.fftshift(FFTt, axes=0), axis=1)
    #FFT = np.flip(np.fft.fftshift(FFTx, axes=1), axis=1)
    FFK = np.absolute(FFTx)

    # Get the frequency and K vectors
    frqv = np.fft.fftfreq(zpt, st*dt)
    wavv = np.fft.fftfreq(zps, sx*dx)
    # Centering
    wavv = np.fft.fftshift(wavv)
    frqv = np.fft.fftshift(frqv)
    return frqv,wavv,FFK

def specfkfreqneg3(U,dx,dt,st=4,sx=4):
    Nt=np.size(U[::st,0])
    Nx=np.size(U[0,::sx])
    """
    FK spectrum of traces using the numpy.fft functions.
    """
    # Time padding
    zpt=2**14
    zps=2**13
    Uzp=np.zeros([zpt,zps])
    Uzp[:Nt,int(zps/2-Nx/2):int(zps/2+Nx/2)]=U[::st,::sx]
    #Fourier in time domain 
    FFTt = np.fft.fft(Uzp[::], axis=0)
    #Fourier in space
    FFTx = np.fft.fft(FFTt, axis=1)
    FFT = np.flip(np.fft.fftshift(FFTx, axes=0), axis=0)
    FFT = np.flip(np.fft.fftshift(FFT, axes=1), axis=1)
    FFK = np.absolute(FFT)

    # Get the frequency and K vectors
    frqv = np.fft.fftfreq(zpt, st*dt)
    wavv = np.fft.fftfreq(zps, sx*dx)
    # Centering
    wavv = np.fft.fftshift(wavv)
    frqv = np.fft.fftshift(frqv)
    return frqv,wavv,FFK

def dispFK(datas,dt,dx,tau=1,c_0=1,norm='No'):
    frqv,wavv,FFK=specfk(datas,dx,dt)
    plt.figure()
    plt.pcolor(2*np.pi*wavv[::]*tau*c_0,2*np.pi*frqv[::]*tau,np.log(FFK[::,::]))
    plt.colorbar()
    if tau==1 and c_0==1 and norm=='No':
        format_fig(r"$k$ (rad/m)",r"$\omega$ (rad/s)")
    else:
        plt.xlim([0,np.pi])
        plt.ylim([0,2*np.pi])
        format_fig(r"$kc_0\tau_m$",r"$\omega\tau_m$ (rad)")
    return

def dispFKB(datas,dt,dx,tau=1,c_0=1,norm='No'):
    frqv,wavv,FFK=specfkfreqnegB(datas,dx,dt)
    plt.figure()
    plt.pcolor(2*np.pi*wavv[::]*tau*c_0,2*np.pi*frqv[::]*tau,np.log(FFK[::,::]))
    plt.colorbar()
    if tau==1 and c_0==1 and norm=='No':
        format_fig(r"$k$ (rad/m)",r"$\omega$ (rad/s)")
    else:
        plt.xlim([0,np.pi])
        plt.ylim([0,2*np.pi])
        format_fig(r"$kc_0\tau_m$",r"$\omega\tau_m$ (rad)")
    return

def dispFKshannon(datas,dt,dx,tau=1,c_0=1,norm='No'):
    frqv,wavv,FFK=specfk(datas,dx,dt)
    plt.figure()
    plt.pcolor(2*np.pi*wavv[::]*tau*c_0,2*np.pi*frqv[::]*tau,np.log(FFK[::,::]))
    plt.colorbar()
    if tau==1 and c_0==1 and norm=='No':
        format_fig(r"$k$ (rad/m)",r"$\omega$ (rad/s)")
    else:
        plt.xlim([0,np.pi])
        plt.ylim([0,2*np.pi])
        format_fig(r"$kc_0\tau_m$",r"$\omega\tau_m$ (rad)")
    return

def dispFKopti(datas,dt,dx,fmax=1000,kmax=0.15,norm='No'):
    dtopti=1/fmax
    if dt<dtopti:
        nskipt=int(dtopti//dt)
    else:
        nskipt=1
    dxopti=1/kmax/2
    if dx<dxopti:
        nskipx=int(dxopti//dx)
    else:
        nskipx=1
    frqv,wavv,FFK=specfk(datas[::,::],dx,dt,st=nskipt,sx=nskipx)
    plt.figure()
    #plt.pcolor(2*np.pi*wavv[numk//2-Nk:numk//2+Nk:],2*np.pi*frqv[:Nf:],np.log(FFK[:Nf:,numk//2-Nk:numk//2+Nk:]),cmap="inferno")
    plt.pcolor(2*np.pi*wavv[:],2*np.pi*frqv[::],np.log(FFK[::,:]),cmap="inferno")
    plt.colorbar()
    if norm=='No':
        format_fig(r"$k$ (rad/m)",r"$\omega$ (rad/s)")
    return frqv,wavv,FFK

def dispFKRed(datas,dt,dx,fm,kmax,norm='No'):
    frqv,wavv,FFK=specfk(datas,dx,dt)
    numk=np.size(wavv)
    df=frqv[2]-frqv[1]
    Nf=int(fm/df)+1
    dk=wavv[2]-wavv[1]
    Nk=int(kmax/dk)+1
    plt.figure()
    plt.pcolor(wavv[numk//2:numk//2+Nk:],frqv[:Nf:],np.log(FFK[:Nf:,:Nk:]),cmap="Blues")
    plt.colorbar()
    plt.xlim([0,kmax])
    format_fig(r"$k$ (1/m)",r"$f$ (Hz)")
    # else:
    #     plt.xlim([0,np.pi])
    #     plt.ylim([0,2*np.pi])
    #     format_fig(r"$kc_0\tau_m$",r"$\omega\tau_m$ (rad)")
    return

def dispFK_tperiodic(datas,dt,dx,fm,kmax,norm='No'):
    frqv,wavv,FFK=specfk(datas,dx,dt)
    numk=np.size(wavv)
    df=frqv[2]-frqv[1]
    Nf=int(fm/df)+1
    dk=wavv[2]-wavv[1]
    Nk=int(kmax/dk)+1
    plt.figure()
    plt.pcolor(wavv[numk//2:numk//2+Nk:],frqv[:Nf:],np.log(FFK[:Nf:,:Nk:]),cmap="Blues")
    plt.colorbar()
    plt.xlim([0,kmax])
    format_fig(r"$k$ (1/m)",r"$f$ (Hz)")
    # else:
    #     plt.xlim([0,np.pi])
    #     plt.ylim([0,2*np.pi])
    #     format_fig(r"$kc_0\tau_m$",r"$\omega\tau_m$ (rad)")
    return

def dispFK_tperiodic2(datas,dt,dx,fm,kmax,norm='No'):
    frqv,wavv,FFK=specfkfreqneg(datas,dx,dt)
    numk=np.size(wavv)
    numf=np.size(frqv)
    df=frqv[2]-frqv[1]
    Nf=int(fm/df)+1
    dk=wavv[2]-wavv[1]
    Nk=int(kmax/dk)+1
    plt.figure()
    plt.pcolor(wavv[numk//2:numk//2+Nk:],frqv[numf//2-Nf//2:numf//2+Nf//2:],np.log(FFK[numf//2-Nf//2:numf//2+Nf//2:,:Nk:]),cmap="Blues",rasterized=True,edgecolors='face')
    plt.colorbar()
    plt.xlim([0,kmax])
    format_fig(r"$k$ (1/m)",r"$f$ (Hz)")
    # else:
    #     plt.xlim([0,np.pi])
    #     plt.ylim([0,2*np.pi])
    #     format_fig(r"$kc_0\tau_m$",r"$\omega\tau_m$ (rad)")
    return


def dispFK_tperiodic3(datas,dt,dx,fm,kmax,norm='No'):
    frqv,wavv,FFK=specfkfreqneg3(datas,dx,dt)
    numk=np.size(wavv)
    numf=np.size(frqv)
    df=frqv[2]-frqv[1]
    Nf=int(fm/df)+1
    dk=wavv[2]-wavv[1]
    Nk=int(kmax/dk)+1
    plt.figure()
    plt.pcolor(wavv[numk//2-Nk:numk//2+Nk:],frqv[numf//2-Nf//2:numf//2+Nf//2:],np.log(FFK[numf//2-Nf//2:numf//2+Nf//2:,numk//2-Nk:numk//2+Nk:]),cmap="Blues")
    plt.colorbar()
    plt.xlim([-kmax,kmax])
    format_fig(r"$k$ (1/m)",r"$f$ (Hz)")
    # else:
    #     plt.xlim([0,np.pi])
    #     plt.ylim([0,2*np.pi])
    #     format_fig(r"$kc_0\tau_m$",r"$\omega\tau_m$ (rad)")
    return

def spectk(U,dx):
    Nt=np.size(U[::4,0])
    Nx=np.size(U[0,:])
    """
    FK spectrum of traces using the numpy.fft functions.
    """
    # Time padding
    zp=2**14
    Uzp=np.zeros([Nt,zp])
    Uzp[:,int(zp/2-Nx/2):int(zp/2+Nx/2)]=U[::4,:]
    #Fourier in time domain 
    #Fourier in sapce
    FFTx = np.fft.fft(Uzp, axis=1)
    FFT = np.flip(np.fft.fftshift(FFTx, axes=1), axis=1)
    FFK = np.absolute(FFT)

    # Get the frequency and K vectors
    wavv = np.fft.fftfreq(zp, dx)
    # Centering
    wavv = np.fft.fftshift(wavv)
    return wavv,FFK

def dispTK(datas,t,dx,tau=1,c_0=1,norm='No',xlim=[0,np.pi]):
    wavv,FFK=spectk(datas,dx)
    plt.figure()
    plt.pcolor(wavv[::]*tau*c_0*2*np.pi,t[:-1:4],np.log(FFK[::,::]),vmin=-6, vmax=6)
    plt.colorbar()
    if tau==1 and c_0==1 and norm=='No':
        format_fig(r"$k$ (1/m)",r"$t$ (s)")
    else:
        format_fig(r"$kc_0\tau_m$",r"$t$ (s)")
    return    
    
# def dispFK():
def VtoU(datas,dt,u0=0):
    Nt=np.size(datas[:,0])
    Nx=np.size(datas[0,:])
    Ufilm=np.zeros([Nt,Nx])
    for i in range(0,Nt-1):
        Ufilm[i+1,:]=Ufilm[i,:]+dt/2*(datas[i+1,:]+datas[i,:])
    return Ufilm


def VtoU2(datas,dt,u0=0):
    Nt=np.size(datas[:,0])
    Nx=np.size(datas[0,:])
    Ufilm=np.zeros([Nt,Nx])
    for i in range(0,Nt-2):
        Ufilm[i+2,:]=Ufilm[i,:]+2*dt/6*(datas[i,:]+4*datas[i+1,:]+datas[i+2,:])
    return Ufilm

def VtoU3(datas,dt,u0=0):
    Nt=np.size(datas[:,0])
    Nx=np.size(datas[0,:])
    Ufilm=np.zeros([Nt,Nx])
    for i in range(0,Nt-3):
        Ufilm[i+3,:]=Ufilm[i,:]+3*dt/8*(datas[i,:]+3*datas[i+1,:]+3*datas[i+2,:]+datas[i+3,:])
    return Ufilm

def VtoU4(datas,dt,u0=0):
    Nt=np.size(datas[:,0])
    Nx=np.size(datas[0,:])
    Ufilm=np.zeros([Nt,Nx])
    for i in range(0,Nt-4):
        Ufilm[i+4,:]=Ufilm[i,:]+4*dt/90*(7*datas[i,:]+32*datas[i+1,:]+12*datas[i+2,:]+32*datas[i+3,:]+7*datas[i+4,:])
    return Ufilm



def UtoV(datas,dt):
    Nt=np.size(datas[:,0])
    Nx=np.size(datas[0,:])
    Vfilm=np.zeros([Nt,Nx])
    for i in range(1,Nt-1):
        Vfilm[i,:]=(datas[i+1,:]-datas[i-1,:])/2/dt
    return Vfilm

def UtoS(datas,dx):
    Nt=np.size(datas[:,0])
    Nx=np.size(datas[0,:])
    Sfilm=np.zeros([Nt,Nx])
    for i in range(1,Nx-1):
        Sfilm[:,i]=(datas[:,i+1]-datas[:,i-1])/2/dx
    return Sfilm



def txt2np(file):
    with open(file, "r") as fichier:
        lignes = fichier.readlines()
    colonne1 = []
    colonne2 = []
    for ligne in lignes:
        mots = ligne.split()
        if len(mots) >= 2:
            try:
                valeur_colonne1 = float(mots[0])
                valeur_colonne2 = float(mots[1])
                colonne1.append(valeur_colonne1)
                colonne2.append(valeur_colonne2)
            except ValueError:
                pass
    LWS0x = np.array(colonne1)
    LWS0y = np.array(colonne2) 
    return LWS0x,LWS0y

