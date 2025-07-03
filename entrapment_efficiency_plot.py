import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.close('all')

# --- Data Arrays ---
Ligt_ar_array = np.array([
    [200, 200, 200, 1000],
    [1005, 1105, 1205, 1104],
    [88, 87.5, 85, 89.5],
])
Ligt_kr_array = np.array([
    [200, 200, 200, 1000],
    [1005, 1105, 1205, 1104],
    [95.5, 86, 87.5, 92],
])    
Ligt_xe_array = np.array([
    [200, 200, 200, 1000],
    [1005, 1105, 1205, 1104],
    [99.5, 90, 88.5, 91.5],
])
sch_ar_array = np.array([
    [30, 79, 127, 118, 40, 42, 241, 231, 219, 35, 2.5, 2.9, 5.2],
    [239+7.9, 219+2.8, 61+0.5, 120+1.0, 47+1.2, 113+2.7, 28+0.1, 48+0.2, 97+0.4, 89+2.6, 37+15.0, 22+7.7, 122+23.5],
    [58, 69, 68, 70, 58, 60, 63, 69, 72, 74, 25, 28, 47]
])
sch_kr_array = np.array([
    [37, 8, 47, 50, 17, 12],
    [226+6.2, 185+24.3, 108+2.3, 53+1.1, 75+4.4, 111+9.3],
    [78, 53, 76, 75, 62, 62]
])
sch_xe_array = np.array([
    [30, 22, 17, 130, 86, 90, 7, 14, 4, 266, 260, 147],
    [33+1.1, 51+2.3, 75+4.4, 33+0.3, 43+0.5, 85+0.9, 75+10.1, 48+3.5, 80+19.0, 86+0.3, 161+0.6, 158+1.1],
    [46, 53, 55, 75, 78, 80, 79, 81, 79, 90, 91, 94]
])
Simon_ar_array = np.array([
    [1.8, 2.5, 2.0, 2.8,6.1,8.3,4.7,2.8,11,13],
    [56+31, 53+22, 51+25, 50+18,77,65,23,68,13.1,56],
    [36, 40, 35, 43,52,54,24,43,36,37]
])
Simon_kr_array = np.array([[], [], []])
Simon_xe_array = np.array([[], [], []])
my_ar_array = np.array([
    [14.8,30.9,48.7,15.0,31.3,47.3,14.4,31.7,48.9],
    [20,20,20,40,40,40,60,60,60],
    [11.8,20.3,21.1,17.1,17.4,32.9,14.6,24.3,28.5]
])
my_kr_array = np.array([
    [22.3,45.2,70.4,23.3,47.2,70.1,22.4,43.3,69.0],
    [20,20,20,40,40,40,60,60,60],
    [20.9,30.6,38.4,23.5,23.8,28.9,22.6,26.9,30.6]
])
my_xe_array = np.array([
    [26.3,51.6,75.9,27.6,54.0,77.8,27.2,56.1,83.2],
    [20,20,20,40,40,40,60,60,60],
    [61.8,72.4,62.5,86.2,85.9,90.6,89.1,91.9,91.6]
])
my2_ar_array = np.array([
    [14.8,30.9,48.7,15.0,31.3,47.3,15.3,31.7,43.6],
    [20,20,20,40,40,40,60,60,60],
    [40.4,40.6,39.1,36.5,41.8,43.2,41.4,47.4,45.4]
])
my2_kr_array = np.array([
    [22.3,45.2,70.4,23.3,47.2,70.1,22.4,46.2,65.9],
    [20,20,20,40,40,40,60,60,60],
    [41.9,50.3,50.9,44.2,47.6,48.9,48.0,49.1,45.4]
])
my2_xe_array = np.array([
    [26.3,51.6,73.8,27.6,54.0,86.4,26.7,56.1,83.2],
    [20,20,20,40,40,40,60,60,60],
    [57.7,67.1,70.1,77.5,77.0,84.7,82.4,86.2,84.2]
])

# --- Plot Setup ---
fig, axes = plt.subplots(3, 2, sharex='col', sharey=True,
                         gridspec_kw={'hspace': 0, 'wspace': 0})

species = ['Argon', 'Krypton', 'Xenon']
sch_data_list = [sch_ar_array, sch_kr_array, sch_xe_array]
my_data_list = [my_ar_array, my_kr_array, my_xe_array]
my2_data_list = [my2_ar_array, my2_kr_array, my2_xe_array]
simon_data_list = [Simon_ar_array, Simon_kr_array, Simon_xe_array]
ligterink_data_list = [Ligt_ar_array, Ligt_kr_array, Ligt_xe_array]

plt.rcParams.update({
    'figure.figsize': (8, 10),
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

# --- Plot the data ---
for i, (name, sch_data, my_data, my2_data, simon_data, ligterink_data) in enumerate(zip(
    species, sch_data_list, my_data_list, my2_data_list, simon_data_list, ligterink_data_list)):

    for label, data, color, marker in [('Schneiderman', sch_data, 'purple', 's'),
                                       ('This work', my_data, 'green', 'o'),
                                       ('This work - denser ice', my2_data, 'blue', 'p'),
                                       ('Simon', simon_data, 'red', '^'),
                                       ('Ligterink', ligterink_data, 'black', 'd')]:
        if data.shape[1] == 0:
            continue
        ratio, ml_total, efficiency = data
        ax_l = axes[i, 0]
        ax_r = axes[i, 1]

        ax_l.scatter(ratio, efficiency, color=color, s=30, marker=marker)
        ax_r.scatter(ml_total, efficiency, color=color, s=30, marker=marker)

        ax_l.set_xscale('log')
        ax_l.set_ylim(5, 105)
        ax_l.set_yticks([20, 40, 60, 80, 100])
        ax_l.set_ylabel(f'{name} [%]')
        ax_l.set_xlim([1.1, 1700])
        ax_l.grid(True, which='both', linestyle='--', linewidth=0.5, color='black')
        if i == 2:
            ax_l.set_xlabel('H$_2$O:noble gas')
        else:
            ax_l.tick_params(labelbottom=False)
        for spine in ax_l.spines.values():
            spine.set_visible(True)
            spine.set_color('black')
            spine.set_linewidth(1)

        ax_r.set_xscale('log')
        ax_r.set_ylim(5, 105)
        ax_r.set_yticks([20, 40, 60, 80, 100])
        ax_r.set_xlim([5, 1700])
        ax_r.grid(True, which='both', linestyle='--', linewidth=0.5, color='black')
        if i == 2:
            ax_r.set_xlabel('Thickness [ML]')
        else:
            ax_r.tick_params(labelbottom=False)
        ax_r.tick_params(labelleft=False)
        for spine in ax_r.spines.values():
            spine.set_visible(True)
            spine.set_color('black')
            spine.set_linewidth(1)

# --- Collective Trendlines ---
colors = {'all': 'gray', 'my': 'green', 'my2': 'blue'}
linestyles = {'all': '-.', 'my': '--', 'my2': ':'}

x_limits_ratio = [1.1, 1400]
x_limits_ml = [1.1, 1400]

for i, (sch_data, my_data, my2_data, simon_data, ligterink_data) in enumerate(zip(
    sch_data_list, my_data_list, my2_data_list, simon_data_list, ligterink_data_list)):

    groups = {
        'all': [sch_data, simon_data, ligterink_data],
        # 'my': [my_data],      # Uncomment to add your own fits
        # 'my2': [my2_data],
    }

    for key, data_list in groups.items():
        x_l, y_l, x_r, y_r = [], [], [], []
        for data in data_list:
            if data.shape[1] == 0:
                continue
            ratio, ml, eff = data
            x_l.extend(ratio)
            y_l.extend(eff)
            x_r.extend(ml)
            y_r.extend(eff)

        # Fit and plot for ratio (left panel)
        if len(x_l) > 1:
            x = np.array(x_l)
            y = np.array(y_l)
            coeffs = np.polyfit(np.log10(x), y, 1)
            x_fit = np.logspace(np.log10(x_limits_ratio[0]), np.log10(x_limits_ratio[1]), 200)
            y_fit = np.polyval(coeffs, np.log10(x_fit))
            axes[i, 0].plot(x_fit, y_fit, linestyle=linestyles[key], color=colors[key])

        # Fit and plot for ML (right panel)
        if len(x_r) > 1:
            x = np.array(x_r)
            y = np.array(y_r)
            coeffs = np.polyfit(np.log10(x), y, 1)
            x_fit = np.logspace(np.log10(x_limits_ml[0]), np.log10(x_limits_ml[1]), 200)
            y_fit = np.polyval(coeffs, np.log10(x_fit))
            axes[i, 1].plot(x_fit, y_fit, linestyle=linestyles[key], color=colors[key])
# --- Legend ---
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='This work', markerfacecolor='green', markersize=8),
    Line2D([0], [0], marker='p', color='w', label='This work - denser ice', markerfacecolor='blue', markersize=8),
    Line2D([0], [0], marker='s', color='w', label='Schneiderman, 2022', markerfacecolor='purple', markersize=8),
    Line2D([0], [0], marker='^', color='w', label='Simon et al., 2023', markerfacecolor='red', markersize=8),
    Line2D([0], [0], marker='d', color='w', label='Ligterink et al., 2024', markerfacecolor='black', markersize=8),
]

axes[1,0].tick_params(axis='both', which='major', length=6, width=1.5, direction='out')
axes[0,0].minorticks_on()

fig.legend(handles=legend_elements, loc='lower center', ncol=2, handletextpad=0.0, columnspacing=0.2, frameon=True, bbox_to_anchor=(0.55, 0.0))
# --- Final Layout ---
plt.subplots_adjust(wspace=0.05, hspace=0.05, left=0.12, right=0.99, top=0.99, bottom=0.20)
plt.show()
