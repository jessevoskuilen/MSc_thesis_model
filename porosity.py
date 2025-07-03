# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 13:51:46 2025

@author: vosku
"""

import numpy as np
import os

# File naming pattern
path = r"C:\Users\vosku\source\repos\MC_sim_Renger\MC_sim_Renger\grid_data"
file_pattern = "empty00{:03d}.dat"  # Replace with your file naming convention
grid_file_pattern = "grid00{:03d}.dat"  # Replace with your file naming convention
num_frames = 1
base_height = 40
height = 200

def read_and_process_file(file_path, grid_file_path):
    """
    Reads a .dat file and extracts data points and neighbor counts.
    """
    contents = np.loadtxt(file_path, dtype=float)
    grid_contents = np.loadtxt(grid_file_path, dtype=float)
    return contents, grid_contents

# Process and plot each frame
for frame in range(1, num_frames+1):
    file_path = os.path.join(path, file_pattern.format(frame))
    grid_file_path = os.path.join(path, grid_file_pattern.format(frame))
    
    # Read and process the data
    contents = np.loadtxt(file_path, dtype=float)
    grid_contents = np.loadtxt(grid_file_path, dtype=float)
    
    ind = 0
    data_points = []
    color_array = []
    count = 0
    for i in range(height):
        for j in range(base_height):
            for k in range(base_height):
                if contents[ind]>0:
                    count += 1
                ind += 1
    print('nempty : ', count, 'spec : ',len(grid_contents[:,0]), 'porosity : ', 100*(count/(count+len(grid_contents[:,0]))), '%')

