# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 15:20:02 2025

@author: abdia
"""

import numpy as np
from numpy import pi, array, tan,linspace
import csv
from sympy import symbols
import matplotlib.pyplot as plt

def F_delta(delta_0,c1,c2):
    X=linspace(0,delta_0,1001)
    deltas,deltak=[],[]
    for x in X:
        deltas.append(c1*x)
        deltak.append(c2*(delta_0-x))
    return X, deltas, deltak

def diagram(X,lista1,lista2):
    plt.plot(X,lista1,label="deltas")
    plt.plot(X,lista2,label="deltak")
    plt.legend()
    plt.ylim(0,2.5e4)
    plt.show()

def maximalist(input_list):
    mops=0
    for l in input_list:
        if abs(float(l))>mops:
            mops=abs(float(l))
    return mops

def c(E,A,L):
    return(E*A/L)

def sigma(F,A):
    return(F/A)

def Moment(Fax,P,µ,d_2,rm):
    M=Fax*(0.16*P+0.58*µ*d_2+µ*rm)
    return M
#Läser filvärden
# =============================================================================

with open('H62_UniformFeedLoad.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    file_values=array(list(spamreader))[1:]
    t_list=(file_values[:,0])
    Fz1_list=(file_values[:,1])

# =============================================================================
# =============================================================================
n= 12 #antal skruvar
hydrosetvikt=2800 #kg
axelvikt=11e3 #kg
L_k = 0.076 #m
P24=3
dw=d=24e-3 #m
dh=d2=22.051e-3 #m
d1=20.752e-3 #m
E_värde=205*1e9 #Pa = 205 GPa
g=9.82 #m/s^27-
# =============================================================================

D_a =L_k+dw
x=(L_k/D_a)**(1/5)
Fg=(hydrosetvikt+axelvikt)*g
Fskruv=Fg/12
Askruv=(pi/16)*(d1+d2)**2
Aekv=pi/4*(dw**2-dh**2)+pi/8*(D_a-dw)*dw*((x+1)**2-1)    

Fs,Fk,F0, F_N=symbols("F_s F_k F_0 F_N")
E_k=symbols("E_k")
µ_min=0.1
µ_max=0.18
r_m= (dw+dh)/4


Fz_max=maximalist(Fz1_list)
F_Nmax=Fg+Fz_max

Fs_max=(F_Nmax)/12
sigmas_max=Fs_max/Askruv
sigmas_a=(Fs_max+Fskruv)/(2*Askruv)


c_s,c_k=c(E_värde,Askruv,L_k),c(E_värde,Aekv,L_k)

D_a=L_k+dw
x=(L_k/D_a)**(1/5)
Aekv=pi/4*(dw**2-dh**2)+pi/8*(D_a-dw)*dw*((x+1)**2-1)
δ_s=Fskruv/c_s
δ_k=Fskruv/c_k

δ_0=δ_s+abs(δ_k)

Fs_max=F_Nmax/12
F0_max=1.05*(c_k/(c_s+c_k)*(Fg))

Z,Fs_list,Fk_list=F_delta(δ_0,c_s,c_k)

diagram(Z, Fs_list, Fk_list)
# =============================================================================
M_totmin=Moment(F0_max,P24,µ_min,d2,r_m)
M_totmax=Moment(F0_max,P24,µ_max,d2,r_m)
print(M_totmax)

F0_maxmin=M_totmin/(0.16*P24)
# =============================================================================