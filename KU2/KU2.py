# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 15:20:02 2025

@author: abdia
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, array
import csv
from sympy import symbols


    

Fz1_list,t_list=[],[]
Fax,Fs,Fk,F0, F_N=symbols("F_ax F_s F_k F_0 F_N")
E_k=symbols("E_k")
µ,µ_b,r_m= symbols("µ µ_b,r_m")
# =============================================================================
n= 12 #antal skruvar
hydrosetvikt=2800 #kg
axelvikt=11e3 #kg
L_k = 0.076 #m
P=3
dw=d=24e-3 #m
dh=d2=22.051e-3
d1=20.752e-3
D_a =L_k+dw
x=(L_k/D_a)**(1/5)
# =============================================================================
E=205*1e9 #Pa = 205 GPa
g=9.82 #m/s^2
Fz=(hydrosetvikt+axelvikt)*g
Fskruv=Fz/12
Askruv=pi*d**2/4
Aekv=pi/4*(dw**2-dh**2)+pi/8*(D_a-dw)*dw*((x+1)**2-1)

#kraftdata=pd.read_excel("file://H62_UniformFeedLoad.xlsx","H62_data",0,"Fz1",1)
#Varuns värden
# =============================================================================

with open('H62_UniformFeedLoad.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')

    file_values=array(list(spamreader))[1:]
    t_list.append(file_values[:,0])
    Fz1_list.append(file_values[:,1])
# =============================================================================

c_s=E*Askruv/L_k
c_k=E*Aekv/L_k
D_a=L_k+dw
x=(L_k/D_a)**(1/5)
Aekv=pi/4*(dw**2-dh**2)+pi/8*(D_a-dw)*dw*((x+1)**2-1)
δ_s=Fskruv/c_s
δ_k=Fk/c_k





c_k=Aekv*E/L_k
ck_lutning=np.arctan(c_k)
cs_lutning=np.arctan(c_s)
Z=(0,10,1001)

for t in t_list:
    pass
# =============================================================================
# M_tot=Fax*(0.16*P+0.58*µ*d2+µ_b*r_m)
# print(M_tot)
# =============================================================================

plt.plot