# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 15:00:39 2025

@author: vosku
"""

import matplotlib.pyplot as plt
import numpy as np


plt.close('all')


array1 = np.array([50.952, 51.373, 50.259, 53.055, 48.262, 51.805, 49.521, 41.693, 33.670, 34.761, 
                   29.696, 25.743, 29.385, 19.739, 16.373, 20.606, 19.341, 20.806, 18.054, 17.563, 
                   25.020, 16.292, 14.430, 11.337, 4.103, 3.771])

array3 = np.array([53.992, 54.964, 49.595, 58.072, 52.612, 50.010, 41.268, 31.426, 28.510, 24.706, 
                    27.133, 21.131, 21.473, 14.804, 13.618, 14.183, 13.716, 16.699, 16.810, 15.595, 
                    12.045, 9.572, 8.075, 4.229, 4.229, 0.378])

array2 = np.array([57.778, 54.063, 48.083, 53.944, 50.917, 47.874, 41.589, 36.588, 44.390, 33.948, 
                    29.346, 29.991, 24.262, 26.946, 18.313, 19.472, 23.919, 20.414, 22.894, 30.340, 
                    28.273, 21.686, 6.409, 6.409, 0.176, 0.000])

# X-axis values (excluding first point for the experiment data)
x_labels = np.array([30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155])
exp_x_labels = np.array([35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155])
# Experiment data (excluding the first value to match dataset size)
experiment_values = np.array([0.48, 0.505, 0.52, 0.54, 0.535, 0.53, 0.56, 0.575, 0.59, 0.59, 0.615, 0.63, 0.65, 0.675, 0.69, 0.725, 0.64, 0.67, 0.705, 0.725, 0.75, 0.765, 0.82, 0.83, 0.73])

# Plot
plt.figure(figsize=(8, 5))
plt.plot(exp_x_labels, experiment_values, 'r--s', label='Experiment', markersize = 7)
plt.plot(x_labels, 0.9*(1-(array2/100)), 'r--s', label='Eincr = 465',  color= 'green', markersize = 7)
plt.plot(x_labels, 0.9*(1-(array1/100)), 'r--s', label='Eincr = 965', color= 'blue', markersize = 7)
plt.plot(x_labels, 0.9*(1-(array3/100)), 'r--s', label='Eincr = 1465',  color= 'magenta', markersize = 7)

plt.xlabel('Deposition Temperature (K)')
plt.ylabel('Density (g/cm3)')
#plt.title('Model vs Experiment - Porosity at Deposition Temperature')
plt.legend()
plt.grid(True)
plt.show()
