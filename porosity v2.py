# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 13:51:46 2025

@author: vosku
"""

import numpy as np
import os

# File naming pattern
path = r"C:\Users\vosku\source\repos\MC_sim_Renger5\MC_sim_Renger\grid_data"
file_pattern = "grid00{:03d}.dat"  # Replace with your file naming convention
Frame = 1
base_height = 60
height = 200

def read_and_process_file(file_path, grid_file_path):
    """
    Reads a .dat file and extracts data points and neighbor counts.
    """
    contents = np.loadtxt(file_path, dtype=float)
    grid_contents = np.loadtxt(grid_file_path, dtype=float)
    return contents, grid_contents

# Process and plot each frame

file_path = os.path.join(path, file_pattern.format(Frame))

# Read and process the data
contents = np.loadtxt(file_path, dtype=float)    
max_values = np.zeros((base_height,base_height))
for j in range(base_height):
    for k in range(base_height):
        for i in range(len(contents)):
            if contents[i,0]==j and contents[i,1]==k:
                if max_values[j,k] < contents[i,2]:
                    max_values[j,k] = contents[i,2]

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Generate grid
x = np.arange(base_height)
y = np.arange(base_height)
X, Y = np.meshgrid(x, y)

# Flatten the arrays
X_flat = X.flatten()
Y_flat = Y.flatten()
Z_flat = max_values.T.flatten()  # transpose to match meshgrid orientation

# Filter out zero values
mask = Z_flat > 0
X_plot = X_flat[mask]
Y_plot = Y_flat[mask]
Z_plot = Z_flat[mask]

# Plotting
#fig = plt.figure(figsize=(10, 7))
#ax = fig.add_subplot(111, projection='3d')
# Use scatter to plot only non-zero points
#scatter = ax.scatter(X_plot, Y_plot, Z_plot, c=Z_plot, cmap='viridis', marker='o')

# Labels and title
#ax.set_xlabel('X')
#ax.set_ylabel('Y')
#ax.set_zlabel('Max Value')
#ax.set_title(f'3D Scatter Plot of max_values > 0 (Frame {Frame})')

# Colorbar
#fig.colorbar(scatter, ax=ax, shrink=0.5, aspect=10)
 
#plt.show()

empty = np.zeros((base_height,base_height))
for i in range(base_height):
    for j in range(base_height):
            empty[i,j] = (max_values[i,j] - np.mod(max_values[i,j],4))/4 
            if max_values[i,j] > 0:
                empty[i,j] = empty[i,j] + 1

empty = sum(sum(empty))
species = len(contents)
porosity = 100*(1-(species/empty))

print('n empty: ', empty, 'spec: ', species, 'Porosity: ',f"{porosity:.3f}", '%')

        