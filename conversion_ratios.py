# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 00:02:31 2025

@author: vosku
"""
Kr_over_Ar = 338/330
Xe_over_Ar = 606/330
Kr_over_Ar /=  (3.89e-6 / 2.04e-9) * (748/1119)
Xe_over_Ar /= (3.89e-6 / 2.19e-10) * (608/1119)

print(Kr_over_Ar,Xe_over_Ar)