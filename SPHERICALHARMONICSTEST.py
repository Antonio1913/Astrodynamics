import numpy as np
from MEAN_PLANETARY_CONSTANTS import Earth as E
# import matplotlib.pyplot as plt
from TOOLS.PROPAGATION_TOOLS import sphericalharmonics
# from mpl_toolkits.mplot3d import Axes3D


# CREATING MESHGRID
dist = E.mu + 1000  # km
long_grid = np.linspace(-np.pi, np.pi, 360)
lat_grid = np.linspace(-np.pi/2+1e-4, np.pi/2-1e-4,180)

[mesh_long, mesh_lat] = np.meshgrid(long_grid, lat_grid)

# Calculate the 3D Cartesian coordinates of the sphere
x = dist * np.cos(mesh_long) * np.cos(mesh_lat)  # X-coordinates
y = dist * np.sin(mesh_long) * np.cos(mesh_lat)  # Y-coordinates
z = dist * np.sin(mesh_lat)                      # Z-coordinates

# Combine into a single 3D array of positions
position_Body = np.stack((x, y, z), axis=-1)

# Extract x, y, z coordinates from position_Body
x = position_Body[:, :, 0]  # Shape: (180, 360)
y = position_Body[:, :, 1]  # Shape: (180, 360)
z = position_Body[:, :, 2]  # Shape: (180, 360)

pos_body = np.vstack((x.flatten(), y.flatten(), z.flatten())).T
HarmonicValues = np.loadtxt(r'D:\ASTRODYNAMICS\EGM2008_Spherical_Harmonics')

accel_state = sphericalharmonics(pos_body, 2, HarmonicValues)

# Create a 3D plot
fig = plt.figure(figsize=(10, 10))
# ax = fig.add_subplot(111, projection='3d')
#
# # Plot the surface
# ax.plot_surface(x, y, z, cmap='viridis', edgecolor='none')
#
# # Set labels and title
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# ax.set_title('3D Sphere Visualization')
#
# # Show the plot
# plt.show()
