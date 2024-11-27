import KEPLER_TOOLS as KT
from MEAN_PLANETARY_CONSTANTS import Earth as E
from TOOLS.PLOTTING_TOOLS import orbitplot
from TOOLS.PROPAGATION_TOOLS import OrbitProp, rkutta4, sphericalharmonics

import numpy as np
import pandas as pd

# Defining types of the outputs
Position = np.ndarray
Velocity = np.ndarray

# Defining the outputs from COE2RV
pos_sat: Position
vel_sat: Velocity

# Defining Unperturbed Orbit
pos_sat, vel_sat = KT.COE2RV(11067.790, 0.832, 87.87*(np.pi/180), 227.89*(np.pi/180), 53.38*(np.pi/180), 92.335*(np.pi/180), E.mu, 'none')


# state_sat_init = (pos_sat.tolist() + vel_sat.tolist())
state_sat_init = np.vstack((pos_sat, vel_sat))
pos_mag = np.linalg.norm(pos_sat.T)
period_sat = 2 * np.pi * (np.sqrt(1000**3 / E.mu))
num_orbit = 1
num_time_int = 1000
time_vec = np.linspace(0, period_sat * num_orbit, num_time_int)
dt = (period_sat * num_orbit) / np.size(time_vec)

pos_states, positions = OrbitProp(time_vec, state_sat_init, E.mu)
position_df = pd.DataFrame(positions)
orbitplot([positions], ['Satellite No-Perturbations'])
state_sat = np.zeros((6, num_time_int))
state_sat[:, 0] = np.array(state_sat_init).flatten()


HarmonicValues = np.loadtxt(r'D:\ASTRODYNAMICS\EGM2008_Spherical_Harmonics\EGM2008')









# def f_sphericalharmonics(y, 5, HarmonicValues):
#     return sphericalharmonics(t, y, HarmonicValues)


# for step in range(num_time_int - 1):
#     # state_sat[step + 1] = rkutta4(sphericalharmonics(state_sat[:, step].reshape(6, 1), 5, HarmonicValues), time_vec[step], state_sat[:, step].reshape(6, 1), dt)
#
#     state_sat[step + 1] = rkutta4(f_sphericalharmonics, time_vec[step], state_sat[:, step].reshape(6, 1), dt)
#

