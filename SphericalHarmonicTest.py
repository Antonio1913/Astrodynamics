import KEPLER_TOOLS as KT
from MEAN_PLANETARY_CONSTANTS import Earth as E
import numpy as np
import pandas as pd
from TOOLS import PLOTTING_TOOLS, PROPAGATION_TOOLS


# Defining types of the outputs
Position = np.ndarray
Velocity = np.ndarray

# Defining the outputs from COE2RV
pos_sat: Position
vel_sat: Velocity

# Defining Unperturbed Orbit
p = 1000
ecc = 0.2
a = p / np.sqrt(1 - ecc**2)
pos_sat, vel_sat = KT.COE2RV(p, ecc, 0*(np.pi/180), 227.89*(np.pi/180), 53.38*(np.pi/180), 92.335*(np.pi/180), E.mu, 53.38*(np.pi/180))

states, positions = PROPAGATION_TOOLS.Unperturbed_Orbit(pos_sat, vel_sat, a)
position_df = pd.DataFrame(positions.T)  # positions - ndarray: (3,1000)
PLOTTING_TOOLS.orbitplot(positions, ['Satellite No-Perturbations'])
# state_sat = np.zeros((6, num_time_int))
# state_sat[:, 0] = np.array(state_sat_init).flatten()


# HarmonicValues = np.loadtxt(r'D:\ASTRODYNAMICS\EGM2008_Spherical_Harmonics\EGM2008')
# harmonic_df = pd.DataFrame(HarmonicValues)
