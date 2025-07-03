import numpy as np
import matplotlib.pyplot as plt
import re

# Load and reshape your data
with open("note.txt", "r") as f:
    content = f.read()

matches = re.findall(r'neighbour distribution\s*=\s*([ \d\n\r\t]+)', content)

all_numbers = []
for match in matches:
    nums = np.fromstring(match, sep=' ', dtype=int)
    all_numbers.extend(nums)

neighbor_distribution = np.array(all_numbers)
if neighbor_distribution.size % 4 != 0:
    raise ValueError("Total elements not divisible by 4 for 4-column reshape.")

data_4col = neighbor_distribution.reshape(-1, 4)

# Convert to percentages
totals = data_4col.sum(axis=1, keepdims=True)
percentages = data_4col / totals * 100

# Plotting style
plt.close('all')
plt.rcParams.update({
    'figure.figsize': (8,5),
    'axes.titlesize': 18,
    'axes.labelsize': 18,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'legend.fontsize': 18,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'lines.linewidth': 1.0,
    'grid.alpha': 0.0,
    'scatter.marker': 's',
    'lines.markersize': 5,
    'font.family': 'DejaVu Sans',
})

# Stacked area plot
x = np.arange(len(percentages)) + 10

# Use evenly spaced colors from the YlGn colormap
cmap = plt.get_cmap('Blues')
colors = [cmap(i) for i in np.linspace(0.2, 1, 4)]  # Avoid too light colors

fig, ax = plt.subplots()
ax.stackplot(x, percentages.T, labels=[f'{i+1} N' for i in range(4)], colors=colors, alpha=0.95)

ax.set_xlabel('Temperature [K]')
ax.set_ylabel('Total [%]')
ax.set_ylim(0, 100)

# Legend at the bottom
ax.legend(loc='upper center', ncol=4, bbox_to_anchor=(0.5, -0.18), frameon=True)
ax.axvline(x=34, color='black', linestyle = 'dotted', linewidth = 1)
ax.axvline(x=62, color='black', linestyle = 'dotted', linewidth = 1)
ax.axvline(x=83, color='black', linestyle = 'dotted', linewidth = 1)

# Example: Custom ticks
plt.tick_params(axis='x', which='major', length=6, width=1.5, direction='out')
plt.tick_params(axis='y', which='major', length=6, width=1.5, direction='out')
plt.xticks(np.arange(0, 160 + 1, 20))  # Customize as needed
plt.yticks(np.arange(0, 100 + 1, 20))    # Y-ticks from 0 to 110, every 10
plt.minorticks_on()
plt.tick_params(which='minor', length=3, color='gray')
plt.tight_layout()
plt.gcf().subplots_adjust(left=0.12, bottom=0.25, right=0.99, top=0.95)
plt.show()
