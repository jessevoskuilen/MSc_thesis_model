# -*- coding: utf-8 -*-
"""
Created on Sun Jun  1 21:27:30 2025
@author: vosku
"""

import os
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

plt.rcParams.update({
    'figure.figsize': (6, 9),
    'axes.titlesize': 18,
    'axes.labelsize': 18,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'legend.fontsize': 18,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'lines.linewidth': 1.0,
    'grid.alpha': 0.0,
    'lines.markersize': 5,
    'font.family': 'DejaVu Sans',
})

def read_and_process_file(file_path):
    contents = np.loadtxt(file_path, dtype=float)
    species_id, temperature, time = [], [], []

    for i in range(contents.shape[0]):
        species_id.append(contents[i, 3])
        time.append(contents[i, 4])
        temperature.append(contents[i, 5])

    return np.array(species_id), time[0], temperature[0]

def main():
    path = r"C:\Users\vosku\source\repos\MC_sim_Renger5\MC_sim_Renger\grid_data"
    file_pattern = "grid0{:04d}.dat"
    num_frames = 1

    time_data = []
    Ar_counts, Kr_counts, Xe_counts = [], [], []

    for frame in range(1, num_frames + 1):
        file_path = os.path.join(path, file_pattern.format(frame))
        if not os.path.exists(file_path):
            print(f"File {file_path} does not exist. Skipping.")
            continue

        print(f"Processing {file_path}...")
        species_id, time_sec, _ = read_and_process_file(file_path)

        Ar = np.sum(species_id == 38)
        Kr = np.sum(species_id == 39)
        Xe = np.sum(species_id == 40)

        Ar_counts.append(Ar)
        Kr_counts.append(Kr)
        Xe_counts.append(Xe)
        time_data.append(time_sec / (1 * 3.154e7))  # Convert seconds to years

    # Save to txt
    output_file = 'Ar_Kr_Xe_vs_time_1K_per_1M.txt'
    data_to_save = np.vstack((time_data, Ar_counts, Kr_counts, Xe_counts)).T
    # np.savetxt(output_file, data_to_save, header='Time_years Ar Kr Xe', fmt='%.6e', comments='')

    print(f"Data saved to {output_file}")

def plot_species_vs_time(filenames, labels=None):
    if labels is None:
        labels = [f"Dataset {i+1}" for i in range(len(filenames))]

    fig, axs = plt.subplots(1, 3, figsize=(18, 5))

    for i, file in enumerate(filenames):
        if not os.path.exists(file):
            print(f"File {file} not found. Skipping.")
            continue

        data = np.loadtxt(file, skiprows=1)
        time, Ar, Kr, Xe = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

        # Convert to percent of starting value
        Ar_percent = (Ar / Ar[0]) * 100
        Kr_percent = (Kr / Kr[0]) * 100
        Xe_percent = (Xe / Xe[0]) * 100

        axs[0].plot(time, Ar_percent, label=labels[i], alpha=0.4)
        axs[1].plot(time, Kr_percent, label=labels[i], alpha=0.4)
        axs[2].plot(time, Xe_percent, label=labels[i], alpha=0.4)

        # Optional markers for 50K time point
        if len(time) > 35:
            axs[0].scatter(time[35], Ar_percent[35], marker='D', s=30, zorder=5)
            axs[1].scatter(time[35], Kr_percent[35], marker='D', s=30, zorder=5)
            axs[2].scatter(time[35], Xe_percent[35], marker='D', s=30, zorder=5)

    from matplotlib.ticker import LogLocator, AutoMinorLocator
    
    for ax in axs:
        ax.set_xscale('log')
        ax.set_xlabel("Time (years)")
        ax.set_ylabel("Retention (%)")
        ax.grid(False)
    
        # ----- X-axis (log scale) -----
        ax.xaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=5))
        ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=100))
        ax.tick_params(axis='x', which='minor', length=4, width=1.5, color='grey', direction='out')
    
        # ----- Y-axis (linear scale) -----
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(axis='y', which='minor', length=4, width=1.5, color='grey', direction='out')
    
        # Optional: turn on major ticks styling too
        ax.tick_params(axis='both', which='major', length=6, width=1.5, direction='out')
    
    # Add legend
    from matplotlib.lines import Line2D
    marker_50K = Line2D([0], [0], marker='D', color='black', linestyle='None', markersize=8, label='T = 50K')

    line_handles, line_labels = axs[0].get_legend_handles_labels()
    handles = [marker_50K] + line_handles
    labels = ['T = 50K'] + line_labels

    from collections import OrderedDict
    unique = OrderedDict(zip(labels, handles))

    fig.legend(
        unique.values(), unique.keys(),
        loc='lower center',
        bbox_to_anchor=(0.5, 0.0),
        ncol=6,
        frameon=True
    )

    plt.tight_layout(rect=[0, 0.1, 1, 1])
    plt.show()


if __name__ == '__main__':
    main()
    plot_species_vs_time(
        [
            'Ar_Kr_Xe_vs_time_1K_per_10Y.txt',
            'Ar_Kr_Xe_vs_time_1K_per_1Y.txt',
            'Ar_Kr_Xe_vs_time_1K_per_1D.txt',
            'Ar_Kr_Xe_vs_time_1K_per_1H.txt',
            'Ar_Kr_Xe_vs_time_1K_per_1M.txt'
        ],
        labels=[
            '1K/10Y',
            '1K/1Y',
            '1K/day',
            '1K/hour',
            '1K/minute'
        ]
    )
