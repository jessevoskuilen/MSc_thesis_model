# -*- coding: utf-8 -*-
"""
Created on Wed May  7 14:32:15 2025

@author: vosku
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

plt.close('all')

plt.rcParams.update({
    'figure.figsize': (16,5),
    'axes.titlesize': 18,
    'axes.labelsize': 18,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'legend.fontsize': 18,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'lines.linewidth': 1.0,
    'grid.alpha': 0.3,
    'scatter.marker': 's',
    'lines.markersize': 5,
    'font.family': 'DejaVu Sans',
})

Ar =[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 1, 6, 5, 4, 5, 11, 9, 10, 24, 27, 22, 23, 23, 18, 12, 2, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 3, 2, 2, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 2, 2, 0, 1, 0, 0, 2, 4, 2, 0, 2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]
t_Ar = [11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0, 111.0, 112.0, 113.0, 114.0, 115.0, 116.0, 117.0, 118.0, 119.0, 120.0, 121.0, 122.0, 123.0, 124.0, 125.0, 126.0, 127.0, 128.0, 129.0, 130.0, 131.0, 132.0, 133.0, 134.0, 135.0, 136.0, 137.0, 138.0, 139.0, 140.0, 141.0, 142.0, 143.0, 144.0, 145.0, 146.0, 147.0, 148.0, 149.0, 150.0, 151.0, 152.0, 153.0, 154.0, 155.0, 156.0, 157.0, 158.0, 159.0, 160.0, 161.0, 162.0, 163.0, 164.0, 165.0, 166.0, 167.0, 167.0, 168.0, 169.0]
Kr = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 1, 2, 1, 3, 4, 2, 5, 5, 5, 10, 3, 7, 13, 12, 8, 17, 12, 14, 4, 2, 2, 1, 0, 1, 0, 1, 0, 0, 0, 0, 3, 0, 3, 0, 3, 1, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 2, 1, 0, 0, 1, 0, 0, 0, 0, 1, 2, 0, 1, 0, 1, 2, 0, 0, 1, 0, 1, 1, 1, 0, 0, 3, 0, 0, 0, 0, 0, 0]
t_Kr = [11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0, 111.0, 112.0, 113.0, 114.0, 115.0, 116.0, 117.0, 118.0, 119.0, 120.0, 121.0, 122.0, 123.0, 124.0, 125.0, 126.0, 127.0, 128.0, 129.0, 130.0, 131.0, 132.0, 133.0, 134.0, 135.0, 136.0, 137.0, 138.0, 139.0, 140.0, 141.0, 142.0, 143.0, 144.0, 145.0, 146.0, 147.0, 148.0, 149.0, 150.0, 151.0, 152.0, 153.0, 154.0, 155.0, 156.0, 157.0, 158.0, 159.0, 160.0, 161.0, 162.0, 163.0, 164.0, 165.0, 166.0, 167.0, 168.0, 169.0, 169.0]
Xe = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 2, 1, 2, 0, 4, 0, 3, 6, 2, 5, 5, 8, 6, 9, 5, 8, 12, 10, 6, 2, 6, 3, 5, 0, 1, 4, 2, 0, 0, 0, 2, 1, 0, 2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 2, 0, 3, 1, 1, 2, 0, 3, 2, 0, 0, 2, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0]
t_Xe = [11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0, 110.0, 111.0, 112.0, 113.0, 114.0, 115.0, 116.0, 117.0, 118.0, 119.0, 120.0, 121.0, 122.0, 123.0, 124.0, 125.0, 126.0, 127.0, 128.0, 129.0, 130.0, 131.0, 132.0, 133.0, 134.0, 135.0, 136.0, 137.0, 138.0, 139.0, 140.0, 141.0, 142.0, 143.0, 144.0, 145.0, 146.0, 147.0, 148.0, 149.0, 150.0, 151.0, 152.0, 153.0, 154.0, 155.0, 156.0, 157.0, 158.0, 159.0, 160.0, 161.0, 162.0, 163.0, 164.0, 165.0, 166.0, 166.0]

# Calculate cumulative volume remaining
Ar_volume = [sum(Ar) - sum(Ar[:i]) for i in range(len(Ar))]
Kr_volume = [sum(Kr) - sum(Kr[:i]) for i in range(len(Kr))]
Xe_volume = [sum(Xe) - sum(Xe[:i]) for i in range(len(Xe))]

# Normalize to 100%
# Convert lists to numpy arrays
Ar = np.array(Ar)
Kr = np.array(Kr)
Xe = np.array(Xe)

# Cumulative volume remaining
Ar_volume = np.cumsum(Ar[::-1])[::-1]
Kr_volume = np.cumsum(Kr[::-1])[::-1]
Xe_volume = np.cumsum(Xe[::-1])[::-1]

# Normalize
Ar_norm = 100 * Ar_volume / np.max(Ar_volume)
Kr_norm = 100 * Kr_volume / np.max(Kr_volume)
Xe_norm = 100 * Xe_volume / np.max(Xe_volume)

Ar_form = [251.,251.,251.,249.,224.,7.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
t_Form_Ar = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170]

Kr_form = [173.,173.,173.,173.,173.,173.,47.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
t_Form_Kr = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170]

Xe_form = [152.,152.,152.,152.,152.,152.,152.,89.,1.,0.,0.,0.,0.,0.,0.,0.,0.]
t_Form_Xe = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170]

# Interpolate 'Ar_form' to match the length of 'Ar_norm'
interp_func_Ar = interp1d(t_Form_Ar, [100 * val / max(Ar_form) for val in Ar_form], kind='linear', fill_value="interpolate")
Ar_form_interp = interp_func_Ar(t_Ar)

interp_func_Kr = interp1d(t_Form_Kr, [100 * val / max(Kr_form) for val in Kr_form], kind='linear', fill_value="interpolate")
Kr_form_interp = interp_func_Kr(t_Kr)

interp_func_Xe = interp1d(t_Form_Xe, [100 * val / max(Xe_form) for val in Xe_form], kind='linear', fill_value="interpolate")
Xe_form_interp = interp_func_Xe(t_Xe)

# Now you can fill the region between the lines
fig, axs = plt.subplots(1, 3, sharex=True)

import seaborn as sns
colors = sns.color_palette("Set2", n_colors=3)
# Plot Ar
axs[0].plot(t_Ar, Ar_norm, color=colors[0],linewidth=2, label='TPD')
axs[0].plot(t_Ar, Ar_form_interp, color='black', linewidth=0.5)
axs[0].fill_between(t_Ar, Ar_norm, Ar_form_interp, color=colors[0], alpha=0.6)
axs[0].scatter(t_Form_Ar, [100 * val / max(Ar_form) for val in Ar_form], color='black', marker='^', label='$T_d$')
axs[0].set_ylabel('Abundance [%]')
axs[0].set_xlabel('Temperature [K]')
axs[0].legend(loc='upper center', bbox_to_anchor=(0.47, -0.18),
            fancybox=True, frameon=True, shadow=False, ncol=2,columnspacing=0.5, fontsize=18)
axs[0].grid(False)

# Plot Kr
axs[1].plot(t_Kr, Kr_norm, color=colors[1],linewidth=2, label='TPD')
axs[1].plot(t_Kr, Kr_form_interp, color='black', linewidth=0.5)
axs[1].fill_between(t_Kr, Kr_norm, Kr_form_interp, color=colors[1], alpha=0.6)
axs[1].scatter(t_Form_Kr, [100 * val / max(Kr_form) for val in Kr_form], color='black', marker='^', label='$T_d$')
axs[1].set_ylabel('Abundance [%]')
axs[1].set_xlabel('Temperature [K]')
axs[1].legend(loc='upper center', bbox_to_anchor=(0.47, -0.18),
            fancybox=True, frameon=True, shadow=False, ncol=2,columnspacing=0.5, fontsize=18)
axs[1].grid(False)

# Plot Xe
axs[2].plot(t_Xe, Xe_norm, color=colors[2], linewidth=2, label='TPD')
axs[2].plot(t_Xe, Xe_form_interp, color='black', linewidth=0.5)
axs[2].fill_between(t_Xe, Xe_norm, Xe_form_interp, color=colors[2], alpha=0.6)
axs[2].scatter(t_Form_Xe, [100 * val / max(Xe_form) for val in Xe_form], color='black', marker='^', label='$T_d$')
axs[2].set_ylabel('Abundance [%]')
axs[2].set_xlabel('Temperature [K]')
axs[2].legend(loc='upper center', bbox_to_anchor=(0.47, -0.18),
            fancybox=True, frameon=True, shadow=False, ncol=2, columnspacing=0.5, fontsize=18)
axs[2].grid(False)

#plt.suptitle('Normalized Remaining Volume Over Time', fontsize=16)
plt.tight_layout()
plt.subplots_adjust(left=0.06, right=0.99, top=0.92, bottom=0.25)

# Now plot Kr/Ar vs Xe/Ar (separate plot for ratios)
Kr_Ar_ratio_Td = np.array(Kr_form[0:10]) / np.array(Ar_form[0:10])
Xe_Ar_ratio_Td = np.array(Xe_form[0:10]) / np.array(Ar_form[0:10])

# Select only datapoints at 10K intervals (10, 20, ..., 100)
desired_temperatures = np.arange(10, 101, 10)

# Helper function to extract values at specific temperatures
def extract_values(temperatures, values, desired_temps):
    return [values[i] for i, t in enumerate(temperatures) if t in desired_temps]

# Extract corresponding values
Ar_selected = extract_values(t_Ar, Ar_volume, desired_temperatures)
Kr_selected = extract_values(t_Kr, Kr_volume, desired_temperatures)
Xe_selected = extract_values(t_Xe, Xe_volume, desired_temperatures)

Kr_selected = np.array(Kr_selected)
Ar_selected = np.array(Ar_selected)
Xe_selected = np.array(Xe_selected)
# Calculate ratios
Kr_Ar_ratio = (Kr_selected/Ar_selected)
Xe_Ar_ratio = (Xe_selected/Ar_selected)

# Scatter plot of Xe/Ar vs Kr/Ar
fig4, axs4 = plt.subplots(figsize=(7, 5))
axs4.scatter(Kr_Ar_ratio_Td[:-2], Xe_Ar_ratio_Td[:-2], color='tab:purple', label='$T_d$')
axs4.plot(Kr_Ar_ratio_Td[:-2], Xe_Ar_ratio_Td[:-2], color='black')
axs4.scatter(Kr_Ar_ratio[:-3], Xe_Ar_ratio[:-3], color='tab:green', label='TPD')
axs4.plot(Kr_Ar_ratio[:-3], Xe_Ar_ratio[:-3], color='black')
axs4.set_xlabel('Kr/Ar')
axs4.set_ylabel('Xe/Ar')
axs4.set_xscale('log')
axs4.set_yscale('log')
axs4.legend(loc='upper center', bbox_to_anchor=(0.47, -0.2),
            fancybox=True, frameon=True, shadow=False, ncol=2, fontsize=18)
axs4.grid(True)
labels = [f"{(i+1)*10}K" for i in range(len(desired_temperatures))]  # ['10K', '20K', ..., '100K']

# Add text annotations
for i in range(3,len(desired_temperatures)-4):
    if i>4:
        axs4.text(Kr_Ar_ratio_Td[i], Xe_Ar_ratio_Td[i]*0.85, labels[i], fontsize=18)  # offset x a bit so text isn't on top of the point
    axs4.text(Kr_Ar_ratio[i], Xe_Ar_ratio[i]*0.85, labels[i], fontsize=18)  # offset x a bit so text isn't on top of the point
  
fig4.tight_layout()

# Optional: Customize tick appearance
gas = ['Ar','Kr','Xe']    
for i in range(3):
    axs[i].minorticks_on()  # turn on minor ticks for both axes
    axs[i].tick_params(axis='both', which='major', length=7)
    axs[i].tick_params(axis='both', which='minor', length=4, color='gray')
    axs[i].set_title(gas[i])
# Show all plots
plt.show()