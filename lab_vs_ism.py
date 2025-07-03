
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 12:05:02 2025

@author: vosku
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import griddata

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
    path = r"C:\Users\vosku\source\repos\MC_sim_Renger4\MC_sim_Renger\grid_data"
    file_pattern = "grid0{:04d}.dat"
    num_frames = 1

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
        print(Ar,Kr,Xe)
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
    output_file = '85K_lab_50ML_nptg120_1000Y_100_12_1_1_1.txt'
    data_to_save = np.vstack((Time_years, Kr_over_Ar, Xe_over_Ar)).T
    np.savetxt(output_file, data_to_save, header='Time_years Kr_over_Ar Xe_over_Ar', fmt='%.6e', comments='')

def plot_multiple_files_clean_final(filenames, labels=None):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogLocator
    from matplotlib.colors import LogNorm
    from matplotlib.cm import ScalarMappable
    from matplotlib.collections import LineCollection
    import numpy as np

    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams.update({
        'figure.figsize': (8,8),
        'font.weight' : 'normal',
        'axes.titlesize': 18,
        'axes.labelsize': 18,
        'xtick.labelsize': 18,
        'ytick.labelsize': 18,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'lines.linewidth': 2,
        'grid.alpha': 0.0,
        'font.family': 'DejaVu Sans',
    })

    def load_data(filename):
        data = np.loadtxt(filename, skiprows=1)
        if data.ndim == 1:
            data = data[np.newaxis, :]
        return data[:, 0], data[:, 1], data[:, 2]

    fig, ax = plt.subplots()

    # --- Reference Data ---
    kr_ar_solar = 2.04e-9 / 3.89e-6
    xe_ar_solar = 2.19e-10 / 3.89e-6
    kr_ar_67p, xe_ar_67p = 0.058, 0.013
    kr_ar_mt, xe_ar_mt = 0.5*(1.8*(0.5*(1/40+1/130))),0.5*(1/40+1/130)
    kr_ar_earth, xe_ar_earth = 2.05e-2, 7e-4
    kr_ar_mars, xe_ar_mars = 2.2e-2, 1.003e-3
    kr_ar_titan, xe_ar_titan = (1/2.06)/10,(1/2.06)/10

    # Create colormap and normalization
    cmap = plt.cm.plasma
    norm = LogNorm(vmin=0.1, vmax=1e7)

    # --- Dataset Loop ---
    ISM_kr_ar,ISM_xe_ar=[],[]
    lab_kr_ar,lab_xe_ar=[],[]
    for idx, filename in enumerate(filenames):
        time, kr_ar, xe_ar = load_data(filename)
        label = labels[idx] if labels else f'Dataset {idx + 1}'

        # Create line segments colored by time
        points = np.array([kr_ar, xe_ar]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
        if idx< 6:
            ISM_kr_ar.append(kr_ar[0])
            ISM_xe_ar.append(xe_ar[0])
        else:
            lab_kr_ar.append(kr_ar[0])
            lab_xe_ar.append(xe_ar[0])

        # Start (blue diamond)
        if 2<idx < 6:
            ax.scatter(kr_ar[0], xe_ar[0], color='#1f77b4', marker='o', s=80, 
                      edgecolor='#1f77b4', linewidth=1.5, zorder=5,
                      label='ISM (this work)' if idx == 3 else None)
        elif idx == 6:
            ax.scatter(kr_ar[0], xe_ar[0], color='#1f77b4', marker='o', s=250, 
                      edgecolor='#ff7f0e', linewidth=3, zorder=5,
                      label=None)
        if idx>8:
            ax.scatter(kr_ar[0], xe_ar[0], color='#ff7f0e', marker='s', s=80, 
                      edgecolor='#ff7f0e', linewidth=1.5, zorder=5,
                      label='Laboratory (this work)' if idx == 9 else None)
      
        # Temperature labels
        if 6 > idx > 2:
             ax.text(kr_ar[0] * 0.63, xe_ar[0]*1.2, label, fontsize=18,
                    )
        elif idx > 8:
            ax.text(kr_ar[0] * 1.12, xe_ar[0]*0.8, label, fontsize=18,
                   )
    # --- Plot reference data ---
  #  ax.scatter(kr_ar_solar, xe_ar_solar, color='gold', marker='*', s=200, 
#              edgecolor='black', linewidth=1, zorder=6)
  #  ax.scatter(kr_ar_earth, xe_ar_earth, color='green', marker='^', s=100, zorder=6)
 #   ax.scatter(kr_ar_mars, xe_ar_mars, color='red', marker='^', s=100, zorder=6)
  #  ax.scatter(kr_ar_mt, xe_ar_mt, color='magenta', marker='s', s=100, zorder=6)
    ax.errorbar(kr_ar_67p, xe_ar_67p, xerr=0.013, yerr=0.003, fmt='o', 
               color='blue', ecolor='blue', elinewidth=1.5, capsize=4, zorder=6)
#    ax.errorbar(kr_ar_mt, xe_ar_mt, xerr=kr_ar_mt-0.5*(0.5*(1/40+1/130)), yerr=xe_ar_mt-(1/130), fmt='none', ecolor='magenta', elinewidth=1, capsize=3, zorder=9)

#    ax.arrow(kr_ar_titan, xe_ar_titan, -0.015, 0, head_width=0.004, head_length=0.004, fc='black', ec='black')  # Left arrow
#    ax.arrow(kr_ar_titan, xe_ar_titan, 0, -0.015, head_width=0.004, head_length=0.004, fc='black', ec='black')  # Down arrow

#    ax.text(kr_ar_solar * 0.64, xe_ar_solar * 1.3, 'Solar', fontsize=18)
#    ax.text(kr_ar_earth * 0.9, xe_ar_earth * 0.78, 'Earth', ha= 'right', fontsize=18, color='green')
#    ax.text(kr_ar_mars * 0.9, xe_ar_mars * 0.78, 'Mars', ha= 'right', fontsize=18, color='red')
    ax.text(kr_ar_67p * 0.95, xe_ar_67p * 0.74, '67P', ha= 'right', fontsize=18, color='blue')
 #   ax.text(kr_ar_titan * 0.98, xe_ar_titan * 0.74, 'Titan', ha= 'right', fontsize=18, color='black')
#    ax.text(kr_ar_mt * 0.95, xe_ar_mt * 0.74, 'Chondrites', ha= 'right', fontsize=18, color='magenta')

    
    alm_t = ['28K','40K','50K', '60K', '70K', '80K']
    alm_kr_ar = np.array([1.36,2.21,5.42,57.3,47.2,80.2])
    alm_xe_ar = np.array([1.50,3.40,5.1,107,587,2830])
    
    for i in range(len(alm_kr_ar)):
        ax.scatter(kr_ar_solar*alm_kr_ar[i],xe_ar_solar*alm_xe_ar[i], color = '#2ca02c', marker = '^',s=80, label = 'Almayrac et al., 2022' if i==0 else None)
        ax.text(kr_ar_solar*alm_kr_ar[i] * 0.6, xe_ar_solar*alm_xe_ar[i]*1.12, alm_t[i], fontsize=18, color = 'black'
              )
    ax.plot(kr_ar_solar*alm_kr_ar,xe_ar_solar*alm_xe_ar, color = '#2ca02c', linewidth=1)
    
    ax.plot(ISM_kr_ar,ISM_xe_ar, color = '#1f77b4', linewidth=1)
    ax.plot(lab_kr_ar,lab_xe_ar, color = '#ff7f0e', linewidth=1)
    
    bn_t = ['30K','38K','50K', '54K', '60K', '70K', '75K']
    bn_kr_ar = np.array([1,16,56,92,365,809,3714])
    bn_xe_ar = np.array([1,7,19,25,83,779,3680])
    
    for i in range(len(bn_kr_ar)):
        ax.scatter(kr_ar_solar*bn_kr_ar[i],xe_ar_solar*bn_xe_ar[i], color = '#d62728', marker = 'D',s=80, label = 'Bar-Nun & Owen, 1998' if i==0 else None)
        if i > 0:
            ax.text(kr_ar_solar*bn_kr_ar[i] * 1.05, xe_ar_solar*bn_xe_ar[i]*0.72, bn_t[i], fontsize=18, color = 'black'
                  )
    ax.plot(kr_ar_solar*bn_kr_ar,xe_ar_solar*bn_xe_ar, color = '#d62728', linewidth=1)
        
    # Axis settings
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("Kr / Ar", fontsize=18)
    ax.set_ylabel("Xe / Ar", fontsize=18)
    ax.xaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.grid(True, alpha=0.0, linestyle='--')

    # Legend for markers
    handles, labels_list = ax.get_legend_handles_labels()
    unique = dict(zip(labels_list, handles))
    ax.legend(unique.values(), unique.keys(), loc='upper center', bbox_to_anchor=(0.45, -0.10),
               fancybox=True, frameon=True, shadow=False, ncol=2, fontsize=18)

        
    plt.tight_layout()
    plt.tick_params(axis='both', which='major', length=6, width=1.5, direction='out')
    plt.tick_params(axis='both', which='minor', length=4, width=1, color='gray', direction='out')   
    plt.minorticks_on()
    plt.subplots_adjust(left=0.13, right=0.97, top=0.98, bottom=0.2)
    plt.show()

# Call the function
main()
plot_multiple_files_clean_final(
    [
        '10K_50ML_nptg120_1000Y_100_12_1_1_1.txt',
        '20K_50ML_nptg120_1000Y_100_12_1_1_1.txt',
        '30K_50ML_nptg120_1000Y_100_12_1_1_1.txt',
        '40K_50ML_nptg120_1000Y_100_12_1_1_1.txt',
        '50K_50ML_nptg120_1000Y_100_12_1_1_1.txt',
        '60K_50ML_nptg120_1000Y_100_12_1_1_1.txt',
        '30K_lab_50ML_nptg120_1000Y_100_12_1_1_1.txt',
        '40K_lab_50ML_nptg120_1000Y_100_12_1_1_1.txt',
        '50K_lab_50ML_nptg120_1000Y_100_12_1_1_1.txt',
        '60K_lab_50ML_nptg120_1000Y_100_12_1_1_1.txt',
        '65K_lab_50ML_nptg120_1000Y_100_12_1_1_1.txt',
        '70K_lab_50ML_nptg120_1000Y_100_12_1_1_1.txt',
        '75K_lab_50ML_nptg120_1000Y_100_12_1_1_1.txt',
        '80K_lab_50ML_nptg120_1000Y_100_12_1_1_1.txt',
        '85K_lab_50ML_nptg120_1000Y_100_12_1_1_1.txt',
    ],
    labels=['10K', '20K', '30K', '40K','50K','60K','30K', '40K', '50K', '60K', '65K', '70K', '75K', '80K', '85K']
)