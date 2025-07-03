import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

# --- Original Argon Data Arrays ---
sch_ar_array = np.array([
    [30, 79, 127, 118, 40, 42, 241, 231, 219, 35, 2.5, 2.9, 5.2],
    [239+7.9, 219+2.8, 61+0.5, 120+1.0, 47+1.2, 113+2.7, 28+0.1, 48+0.2, 97+0.4, 89+2.6, 37+15.0, 22+7.7, 122+23.5],
    [58, 69, 68, 70, 58, 60, 63, 69, 72, 74, 25, 28, 47]
])
Simon_ar_array = np.array([
    [1.8, 2.5, 2.0, 2.8,6.1,8.3,4.7,2.8,11,13],
    [56+31, 53+22, 51+25, 50+18,77,65,23,68,13.1,56],
    [36, 40, 35, 43,52,54,24,43,36,37]
])
Ligt_ar_array = np.array([
    [200, 200, 200, 1000],
    [1005, 1105, 1205, 1104],
    [88, 87.5, 85, 89.5],
])
my_ar_array = np.array([
    [14.8,30.9,48.7,15.0,31.3,47.3,14.4,31.7,48.9],
    [20,20,20,40,40,40,60,60,60],
    [11.8,20.3,21.1,17.1,17.4,32.9,14.6,24.3,28.5]
])
my2_ar_array = np.array([
    [14.8,30.9,48.7,15.0,31.3,47.3,15.3,31.7,43.6],
    [20,20,20,40,40,40,60,60,60],
    [40.4,40.6,39.1,36.5,41.8,43.2,41.4,47.4,45.4]
])

# --- H2O:Ar ratio vs retention efficiency trend ---
ratio = my_ar_array[0]
efficiency = my_ar_array[2]
coeffs_ratio = np.polyfit(np.log10(ratio), efficiency, 1)

def retention_efficiency_from_ratio(ratio_val):
    return np.polyval(coeffs_ratio, np.log10(ratio_val))

# --- Simulated Retention Curves Parameters ---
max_thickness = 1400
block_sizes = my_ar_array[1]
initial_ratios = my_ar_array[0]
base_retention_rates = my_ar_array[2] / 100

# --- Plot Setup ---
plt.figure(figsize=(7, 7))

# --- Plot Experimental Data ---
datasets = [
    ('This work', my_ar_array, 'green', 'o'),
    ('This work - denser ice', my2_ar_array, 'blue', 'p'),
    ('Schneiderman, 2022', sch_ar_array, 'purple', 's'),
    ('Simon et al., 2023', Simon_ar_array, 'red', '^'),
    ('Ligterink et al., 2024', Ligt_ar_array, 'black', 'd'),
]

x_all = []
y_all = []
x_trend = []
y_trend = []

for label, data, color, marker in datasets:
    if data.shape[1] == 0:
        continue
    ml, eff = data[1], data[2]
    x_all.extend(ml)
    y_all.extend(eff)
    plt.scatter(ml, eff, label=label, color=color, marker=marker, s=50)

    # Exclude "This work" and "This work - denser ice" from trendline
    if "This work" not in label:
        x_trend.extend(ml)
        y_trend.extend(eff)

# --- Trendline and 1σ Envelope (excluding "this work") ---
x = np.array(x_trend)
y = np.array(y_trend)
logx = np.log10(x)

coeffs = np.polyfit(logx, y, 1)
y_pred = np.polyval(coeffs, logx)
residuals = y - y_pred
sigma = np.std(residuals)
print(coeffs)
x_fit = np.logspace(np.log10(10), np.log10(1400), 200)
y_fit = np.polyval(coeffs, np.log10(x_fit))
plt.plot(x_fit, y_fit, linestyle='-.', color='gray', label='Trend (Literature Only)')
plt.fill_between(x_fit, y_fit - sigma, y_fit + sigma, color='gray', alpha=0.3)

# --- Simulated Retention ---
for i in range(len(initial_ratios)):
    block_size = block_sizes[i]
    start_ratio = initial_ratios[i]
    base_rate = base_retention_rates[i]

    thicknesses = np.arange(block_size, max_thickness + 1, block_size)
    blocks = len(thicknesses)

    incoming = 1.0
    retained_total = 0
    retention = []
    ratio_current = start_ratio

    for _ in range(blocks):
        efficiency_percent = retention_efficiency_from_ratio(ratio_current)
        adjusted_rate = base_rate * ((efficiency_percent / 100) / base_rate)
        if _ == 0:
            adjusted_rate = base_rate

        retained_here = incoming * adjusted_rate
        retained_total += retained_here
        incoming *= (1 - adjusted_rate)
        retention.append(retained_total * 100)
        ratio_current *= (1 - adjusted_rate)

    plt.plot(thicknesses, retention, linestyle='--', color='black', linewidth=1)

# Dummy handle for extrapolated retention model in legend
plt.plot([], [], linestyle='--', color='black', label='Model (Extrapolated)')

# --- Formatting ---
plt.xlabel('Thickness [ML]')
plt.ylabel('Argon [%]')
plt.xscale('log')
plt.ylim(0, 100)
plt.xlim(10, 1400)
plt.grid(False)

# --- Legend below plot ---
plt.legend(loc='upper center', bbox_to_anchor=(0.48, -0.15),
           fancybox=True, frameon=True, shadow=False, ncol=2)
plt.tick_params(axis='both', which='major', length=6, width=1.5, direction='out')
plt.tick_params(axis='both', which='minor', length=4, width=1, color = 'gray', direction='out')   

plt.minorticks_on()
plt.tight_layout()
plt.show()
