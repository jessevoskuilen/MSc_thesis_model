import pyvista as pv
import numpy as np
import os
# File naming pattern


# File naming pattern
path = r"C:\Users\vosku\source\repos\MC_sim_Renger3\MC_sim_Renger\grid_data"
file_pattern = "empty00001.dat"  # Replace with your file naming convention
# Create a PyVista Plotter
plotter = pv.Plotter(off_screen=True)  # Off-screen rendering for movie creation

# Start the movie
output_movie = "AA_ice_structure_evolution.mp4"
movie_path = os.path.join(path, output_movie)
plotter.open_movie(movie_path, framerate=5)  # Adjust FPS as needed

species_colors = {
    1: [0.0,0.0,0.0],
    2: [0.0, 1.0, 0.0],  # blue for 2 neighbours
    3: [0.0, 0.0, 1.0],  # green for 3 neighbours
    4: [1.0, 1.0, 0.0],  # yellow for 4 neighbours
    5: [1.0, 0.647, 0.0],# orange for 5 neighbours
    6: [1.0, 1.0, 0.0],   #red for 6 neighbours
}

file_path = os.path.join(path, file_pattern.format(file_pattern))

contents = np.loadtxt(file_path, dtype=float)

ind = 0
data_points = []
color_array = []
for i in range(200):
    for j in range(40):
        for k in range(40):
            if contents[ind] > 0:
                data_points.append([k,j,i])
                if contents[ind] > 6:
                    contents[ind] = 6
                color_array.append(species_colors[int(contents[ind])]) # Assign the color from the dictionary
            ind += 1

data_points, color_array = np.array(data_points),np.array(color_array)
# Create a PyVista PointCloud object
point_cloud = pv.PolyData(data_points)

point_cloud.point_data['Color'] = color_array  # Assign to point_data for color

plotter = pv.Plotter()
plotter.add_points(point_cloud, scalars='Color', rgb=True, render_points_as_spheres=True, point_size=10, opacity=1.0)
plotter.show()