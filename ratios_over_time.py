# -*- coding: utf-8 -*-
"""
Created on Thu May 29 10:57:46 2025
@author: vosku
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from matplotlib.collections import LineCollection
import matplotlib.cm as cm

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
    path = r"C:\Users\vosku\source\repos\MC_sim_Renger5\MC_sim_Renger\grid_data"
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

    output_file = '50K_50ML_nptg60_1000Y_100_10_5_5_5.txt'
    data_to_save = np.vstack((Time_years, Kr_over_Ar, Xe_over_Ar)).T
    #np.savetxt(output_file, data_to_save, header='Time_years Kr_over_Ar Xe_over_Ar', fmt='%.6e', comments='')

def plot_multiple_files_with_plasma_lines(filenames, labels=None):
    import matplotlib as mpl
    from matplotlib.ticker import LogLocator

    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams.update({
        'figure.figsize': (10, 8),
        'axes.titlesize': 18,
        'axes.labelsize': 14,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'lines.linewidth': 2,
        'grid.alpha': 0.8,
        'font.family': 'DejaVu Sans',
    })

    def load_data(filename):
        data = np.loadtxt(filename, skiprows=1)
        return data[:, 0], data[:, 1], data[:, 2]

    cmap = plt.get_cmap('plasma')
    markers = ['o', 's', 'D', '^', 'v']
    all_times = []

    # Load all data first
    datasets = []
    for fname in filenames:
        time, kr_ar, xe_ar = load_data(fname)
        datasets.append((time, kr_ar, xe_ar))
        all_times.extend(time)

    norm = mpl.colors.Normalize(vmin=min(all_times), vmax=max(all_times))

    for idx, (time, kr_ar, xe_ar) in enumerate(datasets):
        # Optional jitter to visually separate nearly overlapping lines
        jitter = 1e-5 * idx
        kr_ar += jitter

        points = np.array([kr_ar, xe_ar]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=cmap, norm=norm)
        lc.set_array(time[1:])
        lc.set_linewidth(3)
        lc.set_alpha(0.6)
        #plt.gca().add_collection(lc)

        # Scatter points with increased size and opacity
        plt.scatter(
            kr_ar[-1], xe_ar[-1], c=time[-1], cmap=cmap, norm=norm,
            s=35, alpha=0.8, marker=markers[idx % len(markers)],
            edgecolors='k', linewidths=0.4,
            label=labels[idx] if labels else f'Dataset {idx + 1}'
        )

        # Add label at end of line
        plt.text(
            kr_ar[-1] * 1.05, xe_ar[-1], labels[idx] if labels else f'Dataset {idx + 1}',
            fontsize=10, fontweight='bold', color='black'
        )

    # Reference points
    kr_ar_solar = 2.04e-9 / 3.89e-6
    xe_ar_solar = 2.19e-10 / 3.89e-6
    kr_ar_67p, xe_ar_67p = 0.058, 0.013
    kr_ar_earth, xe_ar_earth = 2.05e-2, 1.105e-3
    kr_ar_mars, xe_ar_mars = 2.1875e-2, 1.6875e-4
    kr_ar_67p_err = [0.013, 0.013]
    xe_ar_67p_err = [0.003, 0.003]

    plt.scatter(kr_ar_solar, xe_ar_solar, color='black', marker='*', s=200)
    plt.scatter(kr_ar_earth, xe_ar_earth, color='green', marker='^', s=100)
    plt.scatter(kr_ar_mars, xe_ar_mars, color='red', marker='^', s=100)

    plt.errorbar(
        kr_ar_67p, xe_ar_67p,
        xerr=np.array([[kr_ar_67p_err[0]], [kr_ar_67p_err[1]]]),
        yerr=np.array([[xe_ar_67p_err[0]], [xe_ar_67p_err[1]]]),
        fmt='o', color='blue', ecolor='blue', elinewidth=1.5, capsize=4,
    )

    # Reference labels
    plt.text(kr_ar_solar * 1.1, xe_ar_solar / 1.1, 'Solar', fontsize=12, fontweight='bold')
    plt.text(kr_ar_earth * 1.1, xe_ar_earth / 1.1, 'Earth', fontsize=12, fontweight='bold', color='green')
    plt.text(kr_ar_mars * 1.1, xe_ar_mars / 1.1, 'Mars', fontsize=12, fontweight='bold', color='red')
    plt.text(kr_ar_67p * 1.1, xe_ar_67p * 1.1, '67P', fontsize=12, fontweight='bold', color='blue')

    # Axis settings
    plt.xlabel("Kr / Ar")
    plt.ylabel("Xe / Ar")
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)

    # Improve tick formatting
    plt.gca().xaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
    plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, numticks=10))

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm)
    cbar.set_label("Time (years)")

    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    #main()
    plot_multiple_files_with_plasma_lines(
        [
            '10K_50ML_nptg60_1000Y_100_10_5_5_5.txt',
            '20K_50ML_nptg60_1000Y_100_10_5_5_5.txt',
            '30K_50ML_nptg60_1000Y_100_10_5_5_5.txt',
            '40K_50ML_nptg60_1000Y_100_10_5_5_5.txt',
            '50K_50ML_nptg60_1000Y_100_10_5_5_5.txt',
        ],
        labels=[
            '10K', '20K', '30K', '40K', '50K'
        ]
    )