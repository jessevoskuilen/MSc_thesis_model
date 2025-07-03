# -*- coding: utf-8 -*-
"""
Created on Thu May 22 11:09:10 2025

@author: vosku
"""
import numpy as np
mf = 0.8
Ebin = [2350,2*2250,3*2350,4*2350]
alpha = 0.43
mark = [0,0,0]
T_tpd = 50
nu = 5*10**12
for T in range(10,180):
    for i in range(len(Ebin)-1):
        Es = (1-alpha)*Ebin[i]
        Eact = alpha*Ebin[i]
        kp = 4*(((Ebin[i+1]-Es)/(Ebin[i]-Es))/(1+((Ebin[i+1]-Es)/(Ebin[i]-Es)))**2)*nu*np.exp(-Eact/T)
        if (mark[i]==0 and kp>(1/T_tpd)):
            mark[i]=1
            print('Temperature: ', T, 'Jump', i+1, 'to', i+2)