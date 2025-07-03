# -*- coding: utf-8 -*-
"""
Created on Tue May 20 12:22:07 2025

@author: vosku
"""

import pyvista as pv
import numpy as np
import os

# File paths
path = r"C:\Users\vosku\source\repos\MC_sim_Renger\MC_sim_Renger\grid_data"
empty_file = "empty00001.dat"
grid_file = "grid00001.dat"

# Species color definitions
empty_species_colors = {
    1: [0.0, 0.0, 0.0],
    2: [0.0, 0.0, 0.0],
    3: [0.0, 0.0, 0.0],
    4: [0.0, 0.0, 0.0],
    5: [1.0, 0.647, 0.0],
    6: [1.0, 1.0, 0.0],
}

grid_species_colors = {
    10: [0.0, 0.0, 1.0],
    38: [0.0, 0.0, 1.0],
}

# -------------------- Load and Process empty00001.dat --------------------
file_path_empty = os.path.join(path, empty_file)
contents_empty = np.loadtxt(file_path_empty, dtype=float)

ind = 0
empty_points = []
empty_colors = []

for i in range(200):
    for j in range(40):
        for k in range(40):
            val = contents_empty[ind]
            if val > 0:
                if val > 6:
                    val = 6
                empty_points.append([k, j, i])
                empty_colors.append(empty_species_colors[int(val)])
            ind += 1

empty_points = np.array(empty_points)
empty_colors = np.array(empty_colors)
empty_cloud = pv.PolyData(empty_points)
empty_cloud.point_data['Color'] = empty_colors

# -------------------- Load and Process grid00001.dat --------------------
file_path_grid = os.path.join(path, grid_file)
contents_grid = np.loadtxt(file_path_grid, dtype=float)

grid_points = contents_grid[:, :3]
species_id = contents_grid[:, 3]

grid_colors = np.zeros((len(species_id), 3))
for i, species in enumerate(species_id):
    if species in grid_species_colors:
        grid_colors[i] = grid_species_colors[species]

grid_cloud = pv.PolyData(grid_points)
grid_cloud.point_data['Color'] = grid_colors

# -------------------- Load and Process gid00001-1.dat --------------------
gid_file = "grid00001-1.dat"
file_path_gid = os.path.join(path, gid_file)
contents_gid = np.loadtxt(file_path_gid, dtype=float)

gid_points = contents_gid[:, :3]
gid_species_id = contents_gid[:, 3]

gid_colors = np.tile([1.0, 0.0, 0.0], (len(gid_species_id), 1))

gid_cloud = pv.PolyData(gid_points)
gid_cloud.point_data['Color'] = gid_colors


# -------------------- Plot Both Datasets --------------------
plotter = pv.Plotter()
plotter.add_points(gid_cloud, scalars='Color', rgb=True, render_points_as_spheres=True, point_size=10, opacity=1.0)
plotter.add_points(empty_cloud, scalars='Color', rgb=True, render_points_as_spheres=True, point_size=15, opacity=1.0)
plotter.add_points(grid_cloud, scalars='Color', rgb=True, render_points_as_spheres=True, point_size=15, opacity=1.0)

plotter.show()
