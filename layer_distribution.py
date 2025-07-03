# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 12:47:32 2025

@author: vosku
"""

import matplotlib.pyplot as plt
import numpy as np
import os

plt.close('all')

# File naming pattern
path = r"C:\Users\vosku\source\repos\MC_sim_Renger\MC_sim_Renger\grid_data"
file_pattern = "grid00078.dat"  # Replace with your file naming convention
# Create a PyVista Plotter

species_colors = {
    10: [0.0, 0.0, 1.0],  # Blue for species 10
    38: [1.0, 0.0, 0.0],  # Red for species 38
    # Add other species IDs and colors as needed
}

file_path = os.path.join(path, file_pattern.format(file_pattern))

contents = np.loadtxt(file_path, dtype=float)

    
layer = np.zeros((int(max(contents[:,2])+1)))
layer_H2O = np.zeros((int(max(contents[:,2])+1)))
layer_Ar = np.zeros((int(max(contents[:,2])+1)))
for i in range(len(contents)):
    layer[int(contents[i,2])] += 1
    if contents[i,3] == 10:
        layer_H2O[int(contents[i,2])] += 1
    if contents[i,3] == 38:
        layer_Ar[int(contents[i,2])] += 1

plt.plot(layer, label = 'Total')
plt.plot(layer_H2O, label = 'H2O')
plt.plot(layer_Ar, label = 'Ar')
plt.legend()
plt.show()