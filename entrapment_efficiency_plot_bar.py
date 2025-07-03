import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Setup
plt.close('all')
plt.rcParams.update({
    'figure.figsize': (16, 5),
    'axes.titlesize': 18,
    'axes.labelsize': 18,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'legend.fontsize': 18,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'font.family': 'DejaVu Sans',
})

# Colors
palette = sns.color_palette("Set2", n_colors=3)
color_map = {
    'Kr': palette[1],
    'Ar': palette[0],
    'Xe': palette[2],
}

# Data
gases = ['Ar', 'Kr', 'Xe']
ratios = ['10:1', '20:1', '30:1']
mls = [20, 40, 60]
n_per_group = len(ratios)
group_spacing = 0.5  # Adjust to control spacing
x = []
for i in range(len(mls)):
    start = i * (n_per_group + group_spacing)
    x.extend(start + np.arange(n_per_group))
x = np.array(x)
width = 0.3

eff_before = {
    'Ar': np.array([11.8, 20.3, 21.1, 17.1, 17.4, 32.9, 14.6, 24.3, 28.5]),
    'Kr': np.array([20.9, 30.6, 38.4, 23.5, 23.8, 28.9, 22.6, 26.9, 30.6]),
    'Xe': np.array([61.8, 72.4, 62.5, 86.2, 85.9, 90.6, 89.1, 91.9, 91.6]),
}
eff_after = {
    'Ar': np.array([40.4, 40.6, 39.1, 36.5, 41.8, 43.2, 41.1,47.4,45.4]),
    'Kr': np.array([41.9, 50.3, 50.9, 44.2, 47.6, 48.9, 48.0, 49.1,45.4]),
    'Xe': np.array([57.7, 67.1, 70.1, 77.5, 77.0, 84.7, 82.4,86.2,84.2]),
}

fig, axes = plt.subplots(1, 3, sharey=True)

for i, gas in enumerate(gases):
    ax = axes[i]
    ax.grid(False)
    before = eff_before[gas]
    after = eff_after[gas]

    ax.bar(x - width/2, before, width=width, label='porous ice', color=color_map[gas], alpha = 0.5)
    ax.bar(x + width/2, after, width=width, label='denser ice', color=color_map[gas])

    for group_idx in range(len(mls)):
        group_start = group_idx * n_per_group
        group_end = group_start + n_per_group
        group_indices = np.arange(group_start, group_end)

        # Get max height from both before and after in the group
        group_max_height = max(np.max(before[group_indices]), np.max(after[group_indices])) + 3

        for j in group_indices:
            delta = after[j] - before[j]
            if not np.isnan(delta) and (before[j] != 0 or after[j] != 0):
                color = 'green' if delta > 0 else 'red' if delta < 0 else 'black'
                ax.text(x[j], group_max_height, f"{delta:+.1f}", ha='center', fontsize=18, color=color, rotation=45)
 
    ax.set_title(f"{gas}")
    ax.set_xticks(x)
    ax.set_xticklabels(ratios * len(mls), rotation=45)
    ax.set_ylim(0, 105)

    ax.set_ylabel("Efficiency [%]")
    ax.yaxis.set_tick_params(labelleft=True)  # Force y-tick labels on all subplots


    # Add ML labels under each ratio group
    for group_idx, ml in enumerate(mls):
        group_center = np.mean(x[group_idx * n_per_group:(group_idx + 1) * n_per_group])
        ax.text(group_center, -20, f"{ml} ML", ha='center', va='top', fontsize=18)
    
    ax.tick_params(axis='y', which='major', length=6, width=1.5, direction='out')
    ax.tick_params(axis='y', which='minor', length=3, width=1.5, color='gray', direction='out')
    ax.minorticks_on()


plt.tight_layout(rect=[0, 0.05, 1, 1])  # No title
# Get handles and labels from each axis
legends = [ax.get_legend_handles_labels() for ax in axes]

# Add three separate legends below the figure
legend_positions = [0.15, 0.45, 0.75]  # Adjust x-positions for layout

for i, (handles, labels) in enumerate(legends):
    fig.legend(handles, labels,
               loc='lower center',
               ncol = 2,
               bbox_to_anchor=(legend_positions[i]+0.05+0.02*i, -0.015),  # Y below 0 to push below plot
               frameon=True,
               fontsize=18,
                  )
    
plt.subplots_adjust(bottom=0.28, left = 0.07, right = 0.99, top = 0.94)  # Make space for legends
plt.show()