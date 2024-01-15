import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Create a grid of values
w_i, w_j = np.meshgrid(np.linspace(0, 10, 30), np.linspace(0, 10, 30))
d_ij = np.random.rand(30, 30)  # Random distances for illustration

# Weighted Product: f(w_i, w_j, d_ij) = w_i * w_j * d_ij
weighted_product = w_i * w_j * d_ij

# Plot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(w_i, w_j, weighted_product, cmap='viridis')

ax.set_xlabel('Weight w_i')
ax.set_ylabel('Weight w_j')
ax.set_zlabel('Weighted Product')
ax.set_title('Weighted Product Metric Visualization')

plt.show()

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Create a grid of values
w_i, w_j = np.meshgrid(np.linspace(0, 100, 300), np.linspace(0, 100, 300))
d_ij = np.random.rand(300, 300)  # Random distances for illustration
alpha = 0.5  # Example exponent

# Non-linear Combination: f(w_i, w_j, d_ij) = (w_i * w_j)^alpha * d_ij
non_linear_combination = (w_i * w_j)**alpha * d_ij

# Plot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(w_i, w_j, non_linear_combination, cmap='plasma')

ax.set_xlabel('Weight w_i')
ax.set_ylabel('Weight w_j')
ax.set_zlabel('Non-linear Combination')
ax.set_title('Non-linear Combination Metric Visualization')

plt.show()

