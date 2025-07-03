import pyvista as pv
import numpy as np
import os

# File naming pattern
path = r"C:\Users\vosku\source\repos\MC_sim_Renger\MC_sim_Renger\grid_data"
file_pattern = "grid00001.dat"  # Replace with your file naming convention
# Create a PyVista Plotter
plotter = pv.Plotter(off_screen=True)  # Off-screen rendering for movie creation

# Start the movie
output_movie = "AA_ice_structure_evolution.mp4"
movie_path = os.path.join(path, output_movie)
plotter.open_movie(movie_path, framerate=5)  # Adjust FPS as needed

species_colors = {
    10: [0.0, 0.0, 1.0],  # Blue for species 10
    38: [1.0, 0.0, 0.0],  # Red for species 38
    # Add other species IDs and colors as needed
}

file_path = os.path.join(path, file_pattern.format(file_pattern))

contents = np.loadtxt(file_path, dtype=float)
data_points = []
species_id = []
time = []
temperature = []
layer_count = np.zeros([int(max(contents[:,2])+1),1])
print(len(contents[:,0]))
for i in range(len(contents[:,0])):
    layer_count[int(contents[i,2])] += 1
 
base_width = 100
max_layer = (base_width*0.25)*(base_width*0.25)

for i in range(len(layer_count)):
    print("layer ", i, ":", (100*layer_count[i]/max_layer)[0], " % filled")
for i in range(contents.shape[0]):
    data_points.append([contents[i, 0], contents[i, 1], contents[i, 2]])
    species_id.append(contents[i, 3])
    time.append(contents[i, 4])
    temperature.append(contents[i,5])
   
# Create a PyVista PointCloud object
point_cloud = pv.PolyData(data_points)

# Map the species IDs to colors
color_array = np.zeros((len(species_id), 3))  # Initialize color array
for i, species in enumerate(species_id):
    if species in species_colors:
        color_array[i] = species_colors[species]  # Assign the color from the dictionary
        
point_cloud.point_data['Color'] = color_array

plotter = pv.Plotter()

plotter.add_points(point_cloud, scalars='Color', rgb=True, render_points_as_spheres=True, point_size=15, opacity=1.0)
#special_point = pv.PolyData([[39,5,7]])
#plotter.add_points(special_point, color='green', render_points_as_spheres=True, point_size=20, opacity=1.0)

plotter.show()