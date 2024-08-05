

import numpy as np
from SOLN_TO_KEPLER.FIndC2C3 import findc2c3
import spiceypy as spice
from TOOLS.PLOTTING_TOOLS import orbitplot
from SPICE.SPICE_TOOLS import load_kernels, ephemdata, tvlist2array
from SOLN_TO_KEPLER.KEPLER import KEPLER
import matplotlib.pyplot as plt




def Lam_universe_var (r0_vec, r_vec, delta_t, t_m, mu):

# Defining Position Magnitudes
    r0_mag = np.linalg.norm(r0_vec)
    r_mag = np.linalg.norm(r_vec)

    cos_delta_v = np.dot(r0_vec, r_vec) / (r0_mag * r_mag)
    delta_v =  np.acos(cos_delta_v)
    sin_delta_v = t_m * np.sqrt(1 - (cos_delta_v*2))
    A = t_m * np.sqrt(r_mag * r0_mag * (1 + cos_delta_v))

    if A == 0:
        raise ValueError("No Solution Available Given the Inputs.")

    psi_n = 0
    c2 = 1/2
    c3 = 1/6
    psi_up = 4 * (np.pi**2)
    psi_low = -4 * np.pi

    for i in range(100):

        y_n = r0_mag + r_mag + (A * ((psi_n * c3) - 1) / np.sqrt(c2))

        # if A and y_n are out of range
        if A > 0 and y_n < 0.0:

            # increasing psi low value
            psi_low += np.pi

            # Adjusting y_n value
            y_n *= -1

        chi_n = np.sqrt(y_n / c2)
        delta_t_n = ((chi_n**3 * c3) + (A * np.sqrt(y_n))) / np.sqrt(mu)

        if delta_t_n <= delta_t:
            psi_low = psi_n
        else:
            psi_up = psi_n

        psi_nplus1 = (psi_up + psi_low) / 2

        c2, c3 = findc2c3(psi_nplus1)

        psi_n = psi_nplus1

        tolerance = 1 * 10 ** -6
        if abs(delta_t_n - delta_t) < tolerance:

            break

    f = 1 - (y_n / r0_mag)
    g_dot = 1 - (y_n / r_mag)
    g = A * np.sqrt(y_n / mu)

#   Defining Velocity Vectors
    v0_vec = (r_vec - (f * r0_vec)) / g
    v_vec = ((g_dot * r_vec) - r0_vec) / g

    return v0_vec, v_vec


load_kernels('solar_system_kernels.tm')

dates = ['2005 Dec 01 00:12:00', '2006 Mar 01 00:00:00']

et0 = spice.str2et(dates[0])
etf = spice.str2et(dates[1])
steps = 10000

time_vec = tvlist2array(et0, etf, steps)
delta_t = etf - et0
t_m = 1
mu = 1.327 * 10 ** 11 # Suns gravitational constant

# body_data= []
#
# body_data.append(ephemdata('EARTH BARYCENTER', time_vec, 'ECLIPJ2000', 'SUN'))
# body_data.append(ephemdata('VENUS BARYCENTER', time_vec, 'ECLIPJ2000', 'SUN'))
# names = ['EARTH BARYCENTER', 'VENUS BARYCENTER']
#
# orbitplot(body_data, names, animate=True)

Earth_vec = ephemdata('EARTH BARYCENTER', time_vec, 'ECLIPJ2000', 'SUN')
Venus_vec = ephemdata('VENUS BARYCENTER', time_vec, 'ECLIPJ2000', 'SUN')
r0_vec = Earth_vec[0, :3]
r_vec = Venus_vec[-1, :3]

v0_vec, v_vec = Lam_universe_var(r0_vec, r_vec, delta_t, t_m, mu)

spacecraft0 = r0_vec.tolist() + v0_vec.tolist()

spacecraft_states = [spacecraft0]
for i in range (1, len(time_vec)):
    changetime = time_vec[i] - time_vec[i-1]
    r_vec_new, v_vec_new = KEPLER(spacecraft_states[-1][:3], spacecraft_states[-1][3:6], changetime, mu)
    spacecraft_states.append(r_vec_new.tolist() + v_vec_new.tolist())

# print(spacecraft_states)
# print(v0_vec)
# print(v_vec)

# Extract position vectors for plotting
spacecraft_states = np.array(spacecraft_states)
positions = spacecraft_states[:, :3]

# Plot the trajectory
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], label='Spacecraft Trajectory')
ax.plot(Earth_vec[:, 0], Earth_vec[:, 1], Earth_vec[:, 2], label='Earth Trajectory')
ax.plot(Venus_vec[:, 0], Venus_vec[:, 1], Venus_vec[:, 2], label='Venus trajectory')


# Plot initial and final positions
ax.plot([r0_vec[0]], [r0_vec[1]], [r0_vec[2]], 'o', label='Initial Position (Earth)')
ax.plot([r_vec[0]], [r_vec[1]], [r_vec[2]], 'o', label='Final Position (Venus)')

# Set equal scaling
max_range = np.array([positions[:, 0].max() - positions[:, 0].min(),
                      positions[:, 1].max() - positions[:, 1].min(),
                      positions[:, 2].max() - positions[:, 2].min()]).max() / 2.0

mid_x = (positions[:, 0].max() + positions[:, 0].min()) * 0.5
mid_y = (positions[:, 1].max() + positions[:, 1].min()) * 0.5
mid_z = (positions[:, 2].max() + positions[:, 2].min()) * 0.5

ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

# Set labels
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')
ax.legend()

# Show the plot
plt.show()