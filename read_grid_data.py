import pyvista as pv
import numpy as np
import os

# File naming pattern
path = r"C:\Users\vosku\source\repos\MC_sim_Renger4\MC_sim_Renger\grid_data"
file_pattern = "grid00{:03d}.dat"  # Replace with your file naming convention
num_frames = 114 # Total number of frames (update based on your dataset)
focus_height = 0
# Create a PyVista Plotter
plotter = pv.Plotter(off_screen=True)  # Off-screen rendering for movie creation

# Start the movie
output_movie = "AA_ice_structure_evolution.mp4"
movie_path = os.path.join(path, output_movie)
plotter.open_movie(movie_path, framerate=5)  # Adjust FPS as needed

species_colors = {
    10: [0.0, 0.0, 1.0],  # Blue for species 10
    38: [1.0, 0.0, 0.0],  # Red for species 39
    # Add other species IDs and colors as needed
}

def read_and_process_file(file_path):
    """
    Reads a .dat file and extracts data points and neighbor counts.
    """
    contents = np.loadtxt(file_path, dtype=float)
    data_points = []
    species_id = []
    time = []
    temperature = []
    
    for i in range(contents.shape[0]):
        data_points.append([contents[i, 0], contents[i, 1], contents[i, 2]])
        species_id.append(contents[i, 3])
        time.append(contents[i, 4])
        temperature.append(contents[i,5])

    return np.array(data_points), species_id, time, temperature

# Process and plot each frame
for frame in range(1, num_frames + 1):
    file_path = os.path.join(path, file_pattern.format(frame))
    if not os.path.exists(file_path):
        print(f"File {file_path} does not exist. Skipping.")
        continue
    
    print(f"Processing {file_path}...")
    
    # Read and process the data
    data_points, species_id, time, temperature = read_and_process_file(file_path)
    
    # Create a PyVista PointCloud object
    point_cloud = pv.PolyData(data_points)

    # Map the species IDs to colors
    color_array = np.zeros((len(species_id), 3))  # Initialize color array
    for i, species in enumerate(species_id):
        if species in species_colors:
            color_array[i] = species_colors[species]  # Assign the color from the dictionary

    # Assign the color array to the points directly (correct method)
    point_cloud.point_data['Color'] = color_array  # Assign to point_data for color

    # Clear the previous frame
    plotter.clear()

    # Add points for the current frame with RGB color and opacity
    plotter.add_points(point_cloud, scalars='Color', rgb=True, render_points_as_spheres=True, point_size=7, opacity=1.0)

    # Calculate the percentage of red vs green atoms (species ID 39 and 10)
    Ar_atoms = np.sum(np.array(species_id) == 38)
    H2O_atoms = np.sum(np.array(species_id) == 10)
    total_atoms = len(species_id)
    Ar_percentage = (Ar_atoms / total_atoms) * 100
    H2O_percentage = (H2O_atoms / total_atoms) * 100
    
    # Extract x, y, z coordinates
    x, y, z = data_points[:, 0], data_points[:, 1], data_points[:, 2]
    
    # Get unique (x, y) pairs and their indices
    unique_xy, inverse_indices = np.unique(data_points[:, :2], axis=0, return_inverse=True)
    
    # Create an array to store the max z-values
    max_z_values = np.zeros(unique_xy.shape[0])
    
    # Efficiently accumulate the max z-values for each unique (x, y)
    np.maximum.at(max_z_values, inverse_indices, z)

    # Compute the median height of the top surface
    avg_top_height = np.mean(max_z_values)
    print(avg_top_height)
    porosity = (1-(total_atoms / (avg_top_height*(0.25*0.25*(1+max(data_points[:,0]**2)))))) *100
    print(porosity)
    time_stamp = max(time)
    temp = temperature[0]
    height = max(data_points[:,2])
    base_height = max(data_points[:,0])
    # Add the percentage text on the frame
    percentage_text = f"Ar atoms: {Ar_percentage:.2f}%\n H2O atoms: {H2O_percentage:.2f} % \n Porosity: {porosity:.2f} % \n T: {temp:.2f} K \n H2O density: 2e8 cm-3 \n Kr density: 2e8 cm-3 \n Stack height: {height:.2f} mly"
    plotter.add_text(percentage_text, position='right', font_size=15, color='black')
    
    time_text = f"Time {time_stamp:.2f} [s]"
    plotter.add_text(time_text, position='upper_left', font_size=15, color='black')
    
    if 0.5*max(data_points[:,2]) > focus_height:
        focus_height = 0.5*max(data_points[:,2])
    # Adjust camera position and focal point
    plotter.camera_position = [((1.5*base_height + 2*height), (1.5*base_height + 2*height), base_height), (30, 30, 0.5*height), (0, 0, 1)]  # More zoomed-out view
    plotter.camera.focal_point = [30, 30, 0.5*height]  # Focus on the center of the scene
    

    # Write the current frame to the movie
    plotter.write_frame()

# Finalize and close the movie
plotter.close()
print(f"Movie saved as {output_movie}")
