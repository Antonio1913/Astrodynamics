

import numpy as np
import spiceypy as spice
from TOOLS.PLOTTING_TOOLS import orbitplot
from SPICE.SPICE_TOOLS import load_kernels, ephemdata, tvlist2array
import matplotlib.pyplot as plt
import KEPLER_TOOLS as KT



def Lam_universe_var (r0_vec, r_vec, delta_t, t_m, mu):

# Defining Position Magnitudes
    r0_mag = np.linalg.norm(r0_vec)
    r_mag = np.linalg.norm(r_vec)

    cos_delta_v = np.dot(r0_vec, r_vec) / (r0_mag * r_mag)
    delta_v = np.acos(cos_delta_v)
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

        c2, c3 = KT.findc2c3(psi_nplus1)

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




# Testing Lambert Code

# Loading Kernels
load_kernels('solar_system_kernels.tm')

# Timeline for Trajectory
dates = ['2005 Dec 01 00:12:00', '2009 Mar 01 00:00:00']

# Time in seconds
et0 = spice.str2et(dates[0])
etf = spice.str2et(dates[1])
steps = 10000

# Time Vector
time_vec = tvlist2array(et0, etf, steps)

# Change in time in seconds
delta_t = etf - et0

# Direction of Trajectory
t_m = -1

# Suns gravitational constant
mu = 1.327 * 10 ** 11

# Ephemeris Data For Planetary Bodies
Earth_vec = ephemdata('EARTH BARYCENTER', time_vec, 'ECLIPJ2000', 'SUN')
Venus_vec = ephemdata('JUPITER BARYCENTER', time_vec, 'ECLIPJ2000', 'SUN')

# Starting and Ending Position for Trajectory Analysis
r0_vec = Earth_vec[0, :3]
r_vec = Venus_vec[-1, :3]

# Running Lambert Code
v0_vec, v_vec = Lam_universe_var(r0_vec, r_vec, delta_t, t_m, mu)

# Assigning the Position and Velocity Vector for Spacecraft
spacecraft0 = r0_vec.tolist() + v0_vec.tolist()

# Assigning the First State of Spacecraft before Propagation
spacecraft_states = [spacecraft0]

# Propagating Spacecrafts Orbit
for i in range(1, len(time_vec)):
    changetime = time_vec[i] - time_vec[i-1]
    r_vec_new, v_vec_new = KT.KEPLER(spacecraft_states[-1][:3], spacecraft_states[-1][3:6], changetime, mu)
    spacecraft_states.append(r_vec_new.tolist() + v_vec_new.tolist())

# Extract position vectors for plotting
spacecraft_states = np.array(spacecraft_states)
positions = spacecraft_states[:, :3]

#Naming and Plotting the Planetary and Spacecraft Data
names = ['Earth', 'Venus', 'Spacecraft']
orbitplot([Earth_vec[:, :3], Venus_vec[:, :3], positions[:, :3]], names)

