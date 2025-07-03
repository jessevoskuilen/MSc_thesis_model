# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 12:05:02 2025

@author: vosku
"""

import os
from scipy.signal import savgol_filter
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import matplotlib.lines as mlines
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib.patches import Polygon
import matplotlib.patches as mpatches

plt.close('all')

def padded_savgol(data, window_length=5, polyorder=0):
    pad = window_length // 2
    padded = np.pad(data, pad_width=pad, mode='reflect')
    smoothed = savgol_filter(padded, window_length, polyorder)
    return smoothed[pad:-pad]

def read_and_process_file(file_path):
    contents = np.loadtxt(file_path, dtype=float)
    species_id, temperature, time = [], [], []

    for i in range(contents.shape[0]):
        species_id.append(contents[i, 3])
        time.append(contents[i, 4])
        temperature.append(contents[i, 5])

    return species_id, time[0], temperature[0]

def main():
    path = r"C:\Users\vosku\source\repos\MC_sim_Renger3\MC_sim_Renger\grid_data"
    file_pattern = "grid0{:04d}.dat"
    num_frames = 10

    Ar_data, Kr_data, Xe_data = [], [], []
    time_data = []
    temp_data = []

    for frame in range(1, num_frames + 1):
        file_path = os.path.join(path, file_pattern.format(frame))
        if not os.path.exists(file_path):
            print(f"File {file_path} does not exist. Skipping.")
            continue

        print(f"Processing {file_path}...")
        species_id, time_sec, temperature = read_and_process_file(file_path)
        species_id = np.array(species_id)

        Ar = np.sum(species_id == 38)
        Kr = np.sum(species_id == 39)
        Xe = np.sum(species_id == 40)

        Ar_data.append(Ar)
        Kr_data.append(Kr)
        Xe_data.append(Xe)
        time_data.append(time_sec / (1 * 3.154e7))  # seconds to years
        temp_data.append(temperature)

    Ar_array = np.array(Ar_data)
    Kr_array = np.array(Kr_data)
    Xe_array = np.array(Xe_data)
    time_array = np.array(time_data)
    mask = Ar_array > 0
    Kr_over_Ar = Kr_array[mask] / Ar_array[mask]
    Xe_over_Ar = Xe_array[mask] / Ar_array[mask]
    Time_years = time_array[mask]
    Kr_over_Ar /=  (3.89e-6 / 2.04e-9) * (748/1119)
    Xe_over_Ar /= (3.89e-6 / 2.19e-10) * (608/1119)

    output_file = 'Dense_ISM_f10h40_50ML_nptg120_1000Y_100_12_1_1_1.txt'
    data_to_save = np.vstack((Time_years, Kr_over_Ar, Xe_over_Ar)).T
    np.savetxt(output_file, data_to_save, header='Time_years Kr_over_Ar Xe_over_Ar', fmt='%.6e', comments='')

def extract_fh_from_filename(filename):
    match = re.search(r"Dense_ISM_f(\d+)h(\d+)", filename)
    return int(match.group(1)), int(match.group(2)) if match else (None, None)

def load_data(filename):
    data = np.loadtxt(filename, skiprows= 1 if (match := re.search(r'f(\d+)h(\d+)', filename)) and match.group(1) == match.group(2) else 2
)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return data[:, 0], data[:, 1], data[:, 2]


def plot_retention_by_formation(filenames):
    # --- Plot setup ---
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams.update({
        'figure.figsize': (9, 8),
        'axes.titlesize': 16,
        'axes.labelsize': 16,
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'font.family': 'DejaVu Sans',
    })

    fig, ax = plt.subplots()

    # --- Reference Data ---
    kr_ar_solar = 2.04e-9 / 3.89e-6
    xe_ar_solar = 2.19e-10 / 3.89e-6
    kr_ar_67p, xe_ar_67p = 0.058, 0.013
    kr_ar_mt, xe_ar_mt = 0.5*(1.8*(0.5*(1/40+1/130))),0.5*(1/40+1/130)
    kr_ar_earth, xe_ar_earth = 2.05e-2, 7e-4
    kr_ar_mars, xe_ar_mars = 2.2e-2, 1.003e-3
    kr_ar_titan, xe_ar_titan = (1/2.06)/10,(1/2.06)/10

    # Reference points
    ref_points = [
        (kr_ar_solar, xe_ar_solar, 'Solar', 'yellow', 'k', '*', 100),
        (kr_ar_earth, xe_ar_earth, 'Earth', 'green', 'green', '^', 80),
        (kr_ar_mars, xe_ar_mars, 'Mars', 'red', 'red', '^', 80),
        (kr_ar_67p, xe_ar_67p, '67P', 'blue', 'blue', 'o', 80),
        (kr_ar_mt, xe_ar_mt, 'Chondrites', 'magenta', 'magenta', 's', 80),
        (kr_ar_titan, xe_ar_titan, 'Titan', 'black', 'black', 'd', 1)
    ]
    
    for kr, xe, label, fc, ec, marker, size in ref_points:
        ax.scatter(kr, xe, color=fc, edgecolor=ec, marker=marker, s=size, zorder=10)
        x,y=0.95,0.74
        if label == 'Solar':
            x,y=1.03,0.68
        elif label == 'Mars' or label == 'Earth':
            y=0.78
        ax.text(kr * x, xe * y, label, fontsize=16, ha='right', color=ec if fc != 'white' else 'black')
    
    ax.arrow(kr_ar_titan, xe_ar_titan, -0.015, 0, head_width=0.004, head_length=0.004, fc='black', ec='black')  # Left arrow
    ax.arrow(kr_ar_titan, xe_ar_titan, 0, -0.015, head_width=0.004, head_length=0.004, fc='black', ec='black')  # Down arrow

    # Add error bars only for 67P
    ax.errorbar(kr_ar_67p, xe_ar_67p, xerr=0.013, yerr=0.003, fmt='none', ecolor='blue', elinewidth=1, capsize=3, zorder=9)
    ax.errorbar(kr_ar_mt, xe_ar_mt, xerr=kr_ar_mt-0.5*(0.5*(1/40+1/130)), yerr=xe_ar_mt-(1/130), fmt='none', ecolor='magenta', elinewidth=1, capsize=3, zorder=9)

    # --- Define color mapping for heating temperatures ---
    heating_temps = [10, 20, 30, 40, 50, 60]
    heating_colors = plt.cm.coolwarm(np.array(np.linspace(0.1,1,6)))
    heating_color_map = dict(zip(heating_temps, heating_colors))
    
    # --- Define linestyle mapping for formation temperatures ---
    formation_temps = [10, 20, 30, 40, 50,60]
    formation_styles = ['-', '--', 'dotted', 'dashdot', (0,(3,1,1,1,1,1)), '-']  # solid, dashed, dash-dot
    formation_style_map = dict(zip(formation_temps, formation_styles))
    
    polygon_info = []  # List to store ellipse details for the inset
    
    # --- Process data ---
    idx=0
    for filename in filenames:
        idx+=1
        f_temp, h_temp = extract_fh_from_filename(filename)
        if f_temp is None or h_temp is None:
            print(f"Could not extract temperatures from {filename}, skipping.")
            continue
            
        # Skip if not in our predefined temperature sets
        if f_temp not in formation_temps or h_temp not in heating_temps:
            print(f"Temperature combination f{f_temp}h{h_temp} not in defined sets, skipping.")
            continue
            
        time, kr_ar, xe_ar = load_data(filename)
        

        color = heating_color_map[f_temp]
        linestyle = formation_style_map[h_temp]
        
        ratio = 0.05
        if  idx < 17:
             ratio = 0.2
        
        def create_oriented_ellipse(kr, xe, width_log=0.1, aspect_ratio=ratio, n_points=50):
            """
            Creates an ellipse oriented along the data trend
            """
            # Remove duplicates and work with unique points
            seen = set()
            kr_unique, xe_unique = [], []
            for i in range(len(kr)):
                if (kr[i], xe[i]) not in seen:
                    seen.add((kr[i], xe[i]))
                    kr_unique.append(kr[i])
                    xe_unique.append(xe[i])
            
            kr_unique = np.array(kr_unique)
            xe_unique = np.array(xe_unique)
            
            if idx < 8: #make overlapping ovals visible with small move
                offset_x = (idx-1)*0.02*np.mean(kr_unique)
                offset_y = (idx-1)*0.02*np.mean(xe_unique)
            elif idx < 17:
                offset_x = np.mod((idx-8),3)*0.1*np.mean(kr_unique)
                offset_y = np.mod((idx-8),3)*0.1*np.mean(xe_unique)
            else:
                offset_x,offset_y = 0,0
            
            # Work in log space
            log_kr = np.log10(kr_unique+offset_x)
            log_xe = np.log10(xe_unique+offset_y)
            
            # Find the principal direction using simple linear fit
    
            p = np.polyfit(log_kr, log_xe, 1)
            angle = np.arctan(p[0])  # angle of the trend line
            
            # Center of the data
            center_kr = np.mean(log_kr) 
            center_xe = np.mean(log_xe) 
            
            # Data span
            data_span = np.sqrt((log_kr.max() - log_kr.min())**2 + 
                               (log_xe.max() - log_xe.min())**2)
                
            # Create rotated ellipse
            t = np.linspace(0, 2*np.pi, n_points)
            vertices = []
            
            for i in range(n_points):
                # Ellipse in local coordinates
                x_local = (data_span/2 + width_log) * np.cos(t[i])
                y_local = (data_span/2 + width_log) * aspect_ratio * np.sin(t[i])
                
                # Rotate to align with data
                x_rot = x_local * np.cos(angle) - y_local * np.sin(angle)
                y_rot = x_local * np.sin(angle) + y_local * np.cos(angle)
                
                # Translate to center and convert to linear space
                kr_point = 10**(center_kr + x_rot)
                xe_point = 10**(center_xe + y_rot)
                
                vertices.append([kr_point, xe_point])
            
            return np.array(vertices)
        
        # Use one of the methods
        vertices = create_oriented_ellipse(kr_ar, xe_ar)
        
        hatch_patterns = {
            10: '---',       # solid fill
            20: '**',
            30: '//',
            40: '...',    # dotted hatch
        }
        hatch = hatch_patterns.get(h_temp, '')
        
        poly = Polygon(vertices,
                       facecolor=color,
                       edgecolor='black',
                       linewidth=1,
                       hatch=hatch,
                       alpha=1.0)  # Adjust alpha if hatches make fill too dark
        ax.add_patch(poly)
            
        polygon_info.append({
            "facecolor": color,
            "edgecolor": 'black',
            "linestyle": linestyle,
            'hatch': hatch,
            "linewidth": 1.5,
            "alpha": 1.0,
            "vertices": vertices,  # ADD THIS LINE to store polygon vertices
            "kr_ar": kr_ar,        # ADD THIS LINE to store original data
            "xe_ar": xe_ar         # ADD THIS LINE to store original data
        })
            

        if f_temp == 10 and h_temp == 40:  # Check if this is the f10h40 dataset
            # Calculate center and angle
            log_kr = np.log10(kr_ar)
            log_xe = np.log10(xe_ar)
            center_kr = 10**np.mean(log_kr)
            center_xe = 10**np.mean(log_xe)
            
            # Get angle from the ellipse orientation

            p = np.polyfit(log_kr, log_xe, 1)
            angle_deg = np.degrees(np.arctan(p[0]))
    
    # --- Axis Settings ---
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("Kr / Ar")
    ax.set_ylabel("Xe / Ar")
    ax.xaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.grid(True, alpha=0.3, linestyle=':')

    formation_temps = [10, 20, 30, 40, 50, 60]
    formation_colors = plt.cm.coolwarm(np.array(np.linspace(0.1,1,6)))
    formation_color_map = dict(zip(formation_temps, formation_colors))

    heating_temps = [10, 20, 30, 40]
    heating_styles = ['-', '--', 'dotted', 'dashdot']
    heating_style_map = dict(zip(heating_temps, heating_styles))

    # Handles and labels for heating ("row 1")
    formation_handles = [
        mlines.Line2D([], [], color='w', marker='o', markerfacecolor=formation_color_map[temp],
                      markeredgecolor='k', markersize=10, linewidth=0) for temp in formation_temps
    ]
    formation_labels = [str(temp) for temp in formation_temps]

    # Handles and labels for formation ("row 2")
    heating_handles = [
        mpatches.Patch(facecolor='white', edgecolor='k',
                       hatch=hatch_patterns[temp], linewidth=1.5)
        for temp in heating_temps
    ]
    heating_labels = [str(temp) for temp in heating_temps]

    fig = plt.gcf()

    # HEATING LEGEND (row 1)
    legend1 = fig.legend(
        [mlines.Line2D([], [], alpha=0, label="Formation (K):")] + [mlines.Line2D([], [], alpha=0, label="Formation (K):")] + formation_handles,
        ["Formation Temperature (K):"]  + [" "] + formation_labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.0),
        ncol=len(formation_handles)+1,
        frameon=True,
        handletextpad=0.3,
        columnspacing=0.3,
        fontsize=16
    )

    # FORMATION LEGEND (row 2)
    legend1 = fig.legend(
        [mlines.Line2D([], [], alpha=0, label="Heating (K):")] + [mlines.Line2D([], [], alpha=0, label="Heating (K):")] + heating_handles,
        ["Heating Temperature (K):"]  + [" "] + heating_labels,
        loc="lower center",
        bbox_to_anchor=(0.485, -0.05),
        ncol=len(heating_handles)+1,
        frameon=False,
        handletextpad=0.3,
        columnspacing=0.3,
        fontsize=16
    )

    plt.gca().add_artist(legend1)  # Ensures both legends are visible

    plt.subplots_adjust(bottom=0.22)  # Make room for stacked legends
    
    zoom_point = (kr_ar_solar*1.11, xe_ar_solar*1.135)
    zoom_width = zoom_point[0] * 0.26        # 50% zoom window width
    zoom_height = zoom_point[1] * 0.24       # 120% zoom window height

    axins = zoomed_inset_axes(ax, zoom=6, loc='lower right', borderpad=0.3)

    for kr, xe, label, fc, ec, marker, size in ref_points:
        axins.scatter(kr, xe, color=fc, edgecolor=ec, marker=marker, s=size+100, linewidths=1, zorder=10)
    axins.scatter(kr_ar_earth, xe_ar_earth, color='green', edgecolor='green', marker='^', s=80, zorder=11)
   
    axins.set_xlim(zoom_point[0] - zoom_width, zoom_point[0] + zoom_width)
    axins.set_ylim(zoom_point[1] - zoom_height, zoom_point[1] + zoom_height)
    axins.set_xscale('log')
    axins.set_yscale('log')

    axins.axes.get_xaxis().set_visible(False)
    axins.axes.get_yaxis().set_visible(False)

    # Draw *ellipses* if their center is inside the inset box!
    xlim = axins.get_xlim()
    ylim = axins.get_ylim()
    count = 0
    for info in polygon_info:
        count += 1
        
        # Check if any of the data points are in the zoom window
        kr_data = info["kr_ar"]
        xe_data = info["xe_ar"]
        
        # Check if at least one point is in the zoom box
        in_zoom = False
        for i in range(len(kr_data)):
            if (xlim[0] <= kr_data[i] <= xlim[1]) and (ylim[0] <= xe_data[i] <= ylim[1]):
                in_zoom = True
                break
        
        if in_zoom:
            # Draw the polygon in the inset
            poly = Polygon(info["vertices"], 
                          facecolor=info["facecolor"], 
                          alpha=info["alpha"], 
                          hatch=info["hatch"],
                          edgecolor=info["edgecolor"], 
                          linestyle='-', 
                          linewidth=info["linewidth"])
            axins.add_patch(poly)
            
  
    for spine in ['top', 'right', 'bottom', 'left']:
        axins.spines[spine].set_visible(True)
        axins.spines[spine].set_linewidth(1.0)
        axins.spines[spine].set_edgecolor("black")
    mark_inset(ax, axins, loc1=2, loc2=3, fc="none", ec="0.6")
    plt.subplots_adjust(left=0.10, right=0.99, top=0.99, bottom=0.19)
    
    plt.show()
# Call the final plotting function the same way:
main()
plot_retention_by_formation(
    [
     'Dense_ISM_f10h40_50ML_nptg120_1000Y_100_12_1_1_1.txt',
     'Dense_ISM_f20h40_50ML_nptg120_1000Y_100_12_1_1_1.txt',
     'Dense_ISM_f40h40_50ML_nptg120_1000Y_100_12_1_1_1.txt',
     ]
)