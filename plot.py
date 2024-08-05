import numpy as np
import matplotlib.pyplot as plt

# Load original data
original_data = np.loadtxt('original_data.dat')
x_orig = original_data[:, 0]
y_orig = original_data[:, 1]
z_orig = original_data[:, 2]

# Load gridded data
gridded_data = np.loadtxt('gridded_data.dat')
x_grid = gridded_data[:, 0].reshape((100, 100))
y_grid = gridded_data[:, 1].reshape((100, 100))
z_grid = gridded_data[:, 2].reshape((100, 100))

# Plot original data
fig, ax = plt.subplots(1, 2, figsize=(14, 6))
scatter = ax[0].scatter(x_orig, y_orig, c=z_orig, cmap='viridis', marker='o')
ax[0].set_title('Original Data')
ax[0].set_xlabel('X')
ax[0].set_ylabel('Y')
fig.colorbar(scatter, ax=ax[0], label='Z')
scatter.set_clim(0, 3000)

# Plot gridded data
c = ax[1].pcolormesh(x_grid, y_grid, z_grid, cmap='viridis', shading='auto')
ax[1].set_title('Gridded Data')
ax[1].set_xlabel('X')
ax[1].set_ylabel('Y')
fig.colorbar(c, ax=ax[1], label='Z')
c.set_clim(0, 3000)

plt.tight_layout()
plt.savefig('interpolation.png')
plt.show()
