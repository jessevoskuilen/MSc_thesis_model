import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

# Configure figure aesthetics
plt.rcParams.update({
    'figure.figsize': (8,5),
    'axes.titlesize': 18,
    'axes.labelsize': 18,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'legend.fontsize': 18,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.spines.left': True,
    'lines.linewidth': 1.5,
    'grid.alpha': 0.0,
    'scatter.marker': '^',
    'lines.markersize': 5,
    'font.family': 'DejaVu Sans',  # or 'serif', 'Arial', etc.
})


data_str = """
0.295 0.307 0.307 0.350 0.454 0.552 0.555 0.570 0.596 0.580 0.611 0.557 0.564 0.623 0.655 0.762 0.875 0.935 0.944 0.944 0.944 0.944 0.944 0.944 0.944 0.944
0.308 0.315 0.300 0.363 0.442 0.566 0.581 0.591 0.624 0.573 0.567 0.628 0.611 0.564 0.657 0.775 0.857 0.936 0.944 0.944 0.944 0.944 0.944 0.944 0.944 0.944
0.323 0.318 0.297 0.345 0.447 0.579 0.579 0.584 0.577 0.608 0.579 0.639 0.622 0.605 0.671 0.750 0.869 0.926 0.944 0.944 0.944 0.944 0.944 0.944 0.944 0.944
0.313 0.307 0.316 0.331 0.438 0.511 0.574 0.587 0.592 0.573 0.620 0.616 0.578 0.635 0.672 0.750 0.863 0.930 0.944 0.944 0.944 0.944 0.944 0.944 0.944 0.944
0.308 0.293 0.300 0.345 0.467 0.530 0.586 0.582 0.606 0.594 0.597 0.583 0.579 0.630 0.689 0.755 0.866 0.930 0.944 0.944 0.944 0.944 0.944 0.944 0.944 0.944
"""

data_values = np.array([list(map(float, row.split())) for row in data_str.strip().split('\n')])

cmap = plt.get_cmap('Blues')
colors = [cmap(i) for i in np.linspace(0.2, 1, 4)]  # Avoid too light colors

cmap = plt.get_cmap('Reds')
colors2 = [cmap(i) for i in np.linspace(0.2, 1, 4)]  # Avoid too light colors

# Compute mean and standard deviation across repetitions
means = np.mean(data_values, axis=0)
stds = np.std(data_values, axis=0)

# X-axis values
x_labels = np.array([30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155])
exp_x_labels = np.array([35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155])
experiment_values = np.array([0.48, 0.505, 0.52, 0.54, 0.535, 0.53, 0.56, 0.575, 0.59, 0.59, 0.615, 0.63, 0.65, 0.675, 0.69, 0.725, 0.64, 0.67, 0.705, 0.725, 0.75, 0.765, 0.82, 0.83, 0.73])

# Alpha scaling based on mean density
alphas = (means - means.min()) / (means.max() - means.min()) * 0.4 + 0.3  # alpha in range [0.5, 1]
alphas2 = (experiment_values - experiment_values.min()) / (experiment_values.max() - experiment_values.min()) * 0.4 + 0.3  # alpha in range [0.5, 1]

# Plot
plt.figure(figsize=(8, 5))

# Error bars
plt.errorbar(x_labels, means, yerr=stds, fmt='none', capsize=0, ecolor='black', zorder=1)
plt.plot(x_labels, means, color='black', linewidth = 0.5, alpha = 0.8)

# Scatter with blue color and alpha based on density
for i in range(len(x_labels)):
    if i == 0:
        plt.scatter(x_labels[i], means[i], color=colors[2], alpha=0.8, marker='o', s=40, zorder=2, label = 'Model')
    else:
        plt.scatter(x_labels[i], means[i], color=colors[2], alpha=0.8, marker='o', s=40, zorder=2)
    for j in range(4):
        plt.scatter(x_labels[i],data_values[j,i], color = colors[2], alpha = 0.4, marker = 'o')
# Add linear trendline for the model data
coeffs = np.polyfit(x_labels, means, 1)
trendline = np.poly1d(coeffs)

# Experiment data
plt.plot(exp_x_labels, experiment_values, color= 'black', linewidth = 0.5, alpha = 0.8)
for i in range(len(exp_x_labels)):
    if i==0:
        plt.scatter(exp_x_labels[i], experiment_values[i], color=colors2[1], alpha=0.8, marker='o', s=40, zorder=2, label = 'Experiment')
    else:
        plt.scatter(exp_x_labels[i], experiment_values[i], color=colors2[1], alpha=0.8, marker='o', s=40, zorder=2)

# Labels and legend
plt.xlabel('Deposition Temperature (K)')
plt.ylabel('Density (g/cm³)')
plt.legend()
plt.grid(False)

plt.tick_params(axis='x', which='major', length=6, width=1.5, direction='out')
plt.tick_params(axis='y', which='major', length=6, width=1.5, direction='out')
plt.xticks(np.arange(30, 160 + 1, 20))  # Customize as needed
plt.yticks(np.arange(0.4,1.0 + 0.01, 0.2))    # Y-ticks from 0 to 110, every 10
plt.minorticks_on()
plt.tick_params(which='minor', length=4, color='gray')
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles,
           labels,
           loc='lower center',
           bbox_to_anchor=(0.5, -0.38),  # Adjust Y further down if overlapping
           frameon=True,
           fontsize=18,
           ncol=2)  # Optional: display in 2 columns
plt.tight_layout()
plt.gcf().subplots_adjust(left=0.12, bottom=0.25, right=0.99, top=0.95)
plt.show()

