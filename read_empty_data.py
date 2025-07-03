import pyvista as pv
import numpy as np
import os

# File naming pattern
path = r"C:\Users\vosku\source\repos\MC_sim_Renger\MC_sim_Renger\grid_data"
file_pattern = "empty00{:03d}.dat"  # Replace with your file naming convention
num_frames = 92  # Total number of frames (update based on your dataset)
focus_height = 0
base_height, height = 40,150
# Create a PyVista Plotter
plotter = pv.Plotter(off_screen=True)  # Off-screen rendering for movie creation

# Start the movie
output_movie = "AA_empty_evolution.mp4"
movie_path = os.path.join(path, output_movie)
plotter.open_movie(movie_path, framerate=5)  # Adjust FPS as needed

species_colors = {
    2: [0.0, 0.0, 1.0],  # blue for 2 neighbours
    3: [0.0, 1.0, 0.0],  # green for 3 neighbours
    4: [1.0, 1.0, 0.0],  # yellow for 4 neighbours
    5: [1.0, 0.647, 0.0],# orange for 5 neighbours
    6: [1.0, 1.0, 0.0],  # red for 6 neighbours
    7: [0.0, 0.0, 0.0],   # black for 7 neighbours or above
}

def read_and_process_file(file_path):
    """
    Reads a .dat file and extracts data points and neighbor counts.
    """
    contents = np.loadtxt(file_path, dtype=float)

    return contents

# Process and plot each frame
for frame in range(4, num_frames + 1):
    file_path = os.path.join(path, file_pattern.format(frame))
    
    print(f"Processing {file_path}...")
    
    # Read and process the data
    contents = np.loadtxt(file_path, dtype=float)

    ind = 0
    data_points = []
    color_array = []
    count = 0
    for i in range(height):
        for j in range(base_height):
            for k in range(base_height):
                if contents[ind]>0:
                    count +=1
                if contents[ind] > 1:
                    data_points.append([k,j,i])
                    colour = int(contents[ind])
                    if colour > 7:
                        colour = 7
                    color_array.append(species_colors[colour]) # Assign the color from the dictionary
                ind += 1
    print(count)
    data_points, color_array = np.array(data_points),np.array(color_array)

    # Create a PyVista PointCloud object
    point_cloud = pv.PolyData(data_points)

    point_cloud.point_data['Color'] = color_array  # Assign to point_data for color

    # Clear the previous frame
    plotter.clear()

    # Add points for the current frame with RGB color and opacity
    plotter.add_points(point_cloud, scalars='Color', rgb=True, render_points_as_spheres=True, point_size=7, opacity=0.5)
    
    # Adjust camera position and focal point
    plotter.camera_position = [(200,200,75), (30,30,0.5*max(data_points[:,2])), (0, 0, 1)]  # More zoomed-out view

    # Write the current frame to the movie
    plotter.write_frame()

# Finalize and close the movie
plotter.close()
print(f"Movie saved as {output_movie}")
