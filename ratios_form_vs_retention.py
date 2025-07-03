# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 12:05:02 2025

@author: vosku
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

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
    num_frames = 62

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

    Kr_over_Ar /= (Kr_over_Ar[0]) * (3.89e-6 / 2.04e-9) * ((927/1270)/Kr_over_Ar[0])
    Xe_over_Ar /= (Xe_over_Ar[0]) * (3.89e-6 / 2.19e-10) * ((703/1270)/Xe_over_Ar[0])

    output_file = '40K_to_100K_50ML_nptg60_2Kmin_100_10_5_5_5.txt'
    data_to_save = np.vstack((Time_years, Kr_over_Ar, Xe_over_Ar)).T
    np.savetxt(output_file, data_to_save, header='Time_years Kr_over_Ar Xe_over_Ar', fmt='%.6e', comments='')

def plot_multiple_files_clean_final(filenames, labels=None):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogLocator
    import numpy as np

    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams.update({
        'figure.figsize': (6, 5),
        'axes.titlesize': 14,
        'axes.labelsize': 14,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'lines.linewidth': 2,
        'grid.alpha': 0.8,
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
    kr_ar_earth, xe_ar_earth = 2.05e-2, 7e-4
    kr_ar_mars, xe_ar_mars = 2.2e-2, 1.003e-3

    ax.scatter(kr_ar_solar, xe_ar_solar, color='white', marker='*', s=30, zorder=5)
    ax.scatter(kr_ar_earth, xe_ar_earth, color='green', marker='^', s=50)
    ax.scatter(kr_ar_mars, xe_ar_mars, color='red', marker='^', s=50)
    ax.errorbar(kr_ar_67p, xe_ar_67p, xerr=0.013, yerr=0.003, fmt='o', color='blue', ecolor='blue', elinewidth=1.5, capsize=4)

    ax.text(kr_ar_solar * 1.1, xe_ar_solar * 0.8, 'Solar', fontsize=14)
    ax.text(kr_ar_earth * 1.1, xe_ar_earth * 0.9, 'Earth', fontsize=14, color='green')
    ax.text(kr_ar_mars * 1.1, xe_ar_mars * 0.9, 'Mars', fontsize=14, color='red')
    ax.text(kr_ar_67p * 1.07, xe_ar_67p * 0.7, '67P', fontsize=14, color='blue')

    # --- Dataset Loop ---
    for idx, filename in enumerate(filenames):
        time, kr_ar, xe_ar = load_data(filename)
        label = labels[idx] if labels else f'Dataset {idx + 1}'

        # Plot full path
        ax.plot(kr_ar, xe_ar, color='gray', alpha=0.4, linewidth=1)

        # Start (blue diamond)
        ax.scatter(kr_ar[0], xe_ar[0], color='blue', marker='D', s=50, edgecolor='k', linewidth=0.5, label='Formation (T = 0 yr)' if idx == 0 else None)

        # End (orange diamond)
        if idx<5:
            ax.scatter(kr_ar[-1], xe_ar[-1], color='orange', marker='D', s=50, edgecolor='k', linewidth=0.5, label='T = 1000 yr' if idx == 0 else None)

        # Optional: label certain points

            ax.text(kr_ar[0] * 1.1, xe_ar[0]*1.15, '10-50K', fontsize=14)
            ax.text(kr_ar[-1] * 1.13, xe_ar[-1]*0.85, label, fontsize=14)     


    # Axis settings
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("Kr / Ar")
    ax.set_ylabel("Xe / Ar")
    ax.xaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.grid(False)

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), loc='upper center', bbox_to_anchor=(0.47, -0.15),
               fancybox=True, frameon=True, shadow=False, ncol=2, fontsize=14)

    plt.tight_layout()
    plt.tick_params(axis='both', which='major', length=6, width=1.5, direction='out')
    plt.tick_params(axis='both', which='minor', length=4, width=1, color = 'gray', direction='out')   
    plt.minorticks_on()
    plt.show()




# Call the final plotting function the same way:
#main()
plot_multiple_files_clean_final(
    [
        '40K_to_100K_50ML_nptg60_2Kmin_100_10_5_5_5.txt'
    ],
    labels=['40K']
)
