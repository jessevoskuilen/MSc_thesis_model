import numpy as np
import matplotlib.pyplot as plt

# Temperature values (K)
temperature = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150])

# Half-life values (years) from different sources
half_life_this_work = np.array([5.2e236, 2.1e109, 6.8e65, 5.2e43, 1.1e32, 1.6e25, 1.8e19, 2.8e14, 3.6e10, 2.3e7, 6.5e4, 2.9e2, 3.9e0, 3.9e-2, 1.05e-3])
half_life_fraser = np.array([8.8e231, 8.6e105, 8.5e63, 8.4e42, 2.1e30, 8.4e21, 8.4e15, 2.7e11, 8.4e7, 1.3e5, 6.8e2, 8.4e0, 2.0e-1, 8.4e-3, 5.3e-4])
half_life_brown = np.array([2.6e189, 1.6e85, 2.9e50, 1.2e33, 4.7e22, 5.3e15, 5.7e10, 1.1e7, 1.4e4, 6.7e1, 8.5e-1, 2.3e-2, 1.0e-3, 7.4e-5, 7.6e-6])

# Plot
plt.figure(figsize=(10, 6))
plt.plot(temperature, half_life_this_work, 'o-', label='This work')
plt.plot(temperature, half_life_fraser, 's-', label='Fraser et al. (2001)')
plt.plot(temperature, half_life_brown, 'd-', label='Brown et al. (2006)')

# Log scale for Y-axis
plt.yscale("log")

# Labels and title
plt.xlabel("Temperature (K)")
plt.ylabel("Half-life (years)")
plt.title("Half-life vs. Temperature")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)

# Show plot
plt.show()