#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 14:59:46 2024

@author: michael
"""

import numpy as np
import matplotlib.pyplot as plt

rho = 1200 ;
mu = 9.408*1e9 ;
h = 10 ;    




# M1 = 38999.999999994
# M2 = 27000.000000006
# K1 = 3076926923.0650887
# K2 = 2352944117.6539793


rho_av = rho ; 
mu_av = mu ; 


kadim1 = K1*h/mu_av 
kadim2 = K2*h/mu_av  
Madim1 = M1/h/rho_av 
Madim2 = M2/h/rho_av 



rho_eff1 = rho_av+M1/h;
rho_eff2 = rho_av+M2/h; 
mu_eff1 = K1*h/(K1*h/mu_av+1); 
mu_eff2 =  K2*h/(K2*h/mu_av+1); 



phi = 0.75 ;  


c1 = np.sqrt(mu_eff1/rho_eff1);
c2 = np.sqrt(mu_eff2/rho_eff2);


c_eta = c1*phi+c2*(1-phi) 
ks = 2*np.pi           
tau = 10*ks*h/c_eta    

Z1 = np.sqrt(rho_eff1*mu_eff1) 
Z2 = np.sqrt(rho_eff2*mu_eff2)


L1 = phi*tau 
L2 = (1-phi)*tau 


def D_micro(k,omega):
    Dmicro=np.cos(omega*tau)-np.cos(c1*k*L1)*np.cos(c2*k*L2)+1/2*(Z1/Z2+Z2/Z1)*np.sin(c1*k*L1)*np.sin(c2*k*L2)
    return Dmicro

Nome=150
Nk=150

Omega_plot = np.linspace(-np.pi,np.pi,Nome)/tau 
k_map = np.linspace(0,4*np.pi,Nk)/tau/c_eta;


DD_micro = np.zeros([Nome,Nk])
DD_det = np.zeros([Nome,Nk])


for ind in range(Nk):
    for ind2 in range(Nome):
        DD_micro[ind2,ind] = D_micro(k_map[ind],Omega_plot[ind2])
        #DD_det(ind2,ind) = Det_micro(k_map(ind),Omega_plot(ind2)



plt.figure()
plt.pcolor(k_map*c_eta*tau,Omega_plot*tau,np.log10(np.abs(DD_micro)))

