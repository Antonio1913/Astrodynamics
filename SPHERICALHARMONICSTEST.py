import numpy as np
from TOOLS.BODY_CONSTANTS import Earth as E
import matplotlib.pyplot as plt
import TOOLS as tls
import os
# from mpl_toolkits.mplot3d import Axes3D


# CREATING MESHGRID
dist = E.mu + 100  # km
long_grid = np.linspace(-np.pi, np.pi, 360)
lat_grid = np.linspace(-np.pi/2+1e-4, np.pi/2-1e-4, 180)

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

# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------

# LOADING THE EGM2008 COEFFICIENTS
egm_location = os.path.abspath(r"ASSETS\\EGM2008.txt")
HarmonicValues = np.loadtxt(egm_location, delimiter=',', skiprows=1)

accel_stateJ2 = tls.sphericalharmonics(pos_body, 2, HarmonicValues)
accel_stateJ2 = np.linalg.norm(accel_stateJ2, axis=1)

# accel_stateJ20 = sphericalharmonics(pos_body, 20, HarmonicValues)
# accel_stateJ20 = np.linalg.norm(accel_stateJ20, axis=1)

gravity_pointmass = (E.mu / dist**2)

# diff_gravity = accel_stateJ20 - accel_stateJ2

# Comparing J2 to point mass
accel = (accel_stateJ2 / gravity_pointmass) / 100

plt.pcolormesh(mesh_long*(180/np.pi), mesh_lat*(180/np.pi), accel.reshape(180, 360))
plt.colorbar()
plt.show()

# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# # Create a 3D plot
# fig = plt.figure(figsize=(10, 10).)
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
