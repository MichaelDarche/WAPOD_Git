#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 11:12:54 2026

@author: michael
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

import sys
import os

sys.path.append(os.path.abspath(".."))

import functions as func

configuration=func.readfile('Demarrer.txt')
mat=func.readfile('Milieu.txt')
frontiere=func.readfile('Frontiere.txt')     
sourcef=func.readfile('Source.txt')   

X,xmin,xmax,Nx,dx,ESIM,Npt=func.geometric(configuration)
fs,CFL,x_0=func.source(sourcef)
Nmat,rho,cm=func.material(mat)
alpha=func.frontiere(frontiere)
rhox,Celx=func.affectMaterials(alpha,rho,cm,xmin,xmax,Nx,dx)
fm = 80.0

nSveRes='Vfilm'+str(fs)+'fm'+str(fm)+'.txt'
nSvetp='tp'+str(fs)+'fm'+str(fm)+'.txt'
 
Vfilm=np.loadtxt(nSveRes, delimiter='\t')
tp=np.loadtxt(nSvetp,delimiter='\t')

dt=tp[1]-tp[0]
Ufilm=func.VtoU(Vfilm,dt)

Ntfig=int(np.floor(0.3/dt))
skipSpacePoints=4

vals = []
for i in range(alpha.size-1):
    if frontiere["Contact_"+str(i)] != "parfait":
        v = func.fonctionModMP(Ntfig*dt, frontiere, i)
        vals.append(v)
    else:
        vals.append(0)
vals = np.array(vals)
cmap = plt.cm.PiYG  
norm = mpl.colors.Normalize(vmin=2*vals.min(), vmax=2*vals.max())
pfilm=5

psmax=np.where(Ufilm[Ntfig-1,:Nx//2]==max(Ufilm[Ntfig-1,:Nx//2]))
xzoommin=X[psmax[0][0]-200]
xzoommax=X[psmax[0][0]+600]
zoomlim=[xzoommin,xzoommax]
##### Total displacement field
plt.figure()
#Color of the modulated interfaces
for i in range(alpha.size):
    if frontiere["Contact_"+str(i)]!="parfait":
        plt.axvline(x=alpha[i],ymin=0,ymax=1,color=cmap(vals[i]),alpha=np.abs(vals[i]/2),linewidth=0.5)
sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
#Displacement field
plt.scatter(X[::1],Ufilm[Ntfig,::1],marker='o',label=r'$U_h$',color='red',zorder=1000,s=30)
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.title("$t= $"+str(round(tp[Ntfig], 3))+" s",fontsize=18)
plt.legend(fontsize=16)
plt.tight_layout()
plt.show()

##### Zoom on a part
plt.figure()
for i in range(alpha.size):
    if frontiere["Contact_"+str(i)]!="parfait":
        plt.axvline(x=alpha[i],ymin=0,ymax=1,color=cmap(vals[i]),alpha=np.abs(vals[i]/1),linewidth=0.8)
sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar =plt.colorbar(sm, ax=plt.gca())
cbar.set_label(r'$\varepsilon_{M,\ell}~\sin (\Omega_{M,\ell} T+\Phi_\ell)$',size=16)
sm.set_clim(-1.0, 1.0)
plt.scatter(X[::skipSpacePoints],Ufilm[Ntfig,::skipSpacePoints],marker='o',label=r'$U_h$',color='red',s=30,zorder=1000)
y0 = plt.ylim()[0]
# Labels 
for args in zip([470,500,530], ['^','*','^'],
                ['#1f77b4','black','#ff7f0e'],['Left sensor','Source','Left sensor']):
    plt.scatter(args[0], y0, marker=args[1], s=100,
                color=args[2],label=args[3] ,zorder=1005, clip_on=False)
plt.xlabel(r'$X$ (m)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.title("$t= $"+str(round(tp[Ntfig], 3))+" s",fontsize=18)
plt.xlim([440,560])
plt.legend(fontsize=16)
plt.tight_layout()
plt.show()

#Non-reciprocal features
plt.figure()
plt.plot(tp[:-1],Ufilm[:,3760],label=r'$U_{left}$',linewidth=3)
plt.plot(tp[:-1],Ufilm[:,4240],label=r'$U_{right}$',linewidth=3)
plt.xlabel(r'$T$ (s)',fontsize=16)
plt.ylabel(r'$U$ (m)',fontsize=16)
plt.legend(fontsize=16)
plt.tight_layout()
plt.show()
