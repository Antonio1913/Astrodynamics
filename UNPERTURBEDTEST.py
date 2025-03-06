from TOOLS.BODY_CONSTANTS import Earth as E
import numpy as np
import pandas as pd
from TOOLS import PLOTTING_TOOLS, PROPAGATION_TOOLS, KEPLER_TOOLS as KT

# Defining types of the outputs
Position = np.ndarray
Velocity = np.ndarray

# Defining the outputs from COE2RV
pos_sat: Position
vel_sat: Velocity

# Defining Unperturbed Orbit
p = 11067.790
ecc = 0.2
a = p / (1 - ecc**2)
pos_sat, vel_sat = KT.COE2RV(p, ecc, 40*(np.pi/180), 227.89*(np.pi/180), 53.38*(np.pi/180), 92.335*(np.pi/180), E.mu, "none")

states, positions = PROPAGATION_TOOLS.Unperturbed_Orbit(pos_sat, vel_sat, a)
position_df = pd.DataFrame(positions.T)  # positions - ndarray: (3,1000)
PLOTTING_TOOLS.orbitplot(positions, ['Satellite No-Perturbations'])


HarmonicValues = np.loadtxt(r'D:\ASTRODYNAMICS\EGM2008_Spherical_Harmonics\EGM2008')
first_state = positions[:, 0].T
first_state = first_state.reshape(3,1)

# accel_state = sphericalharmonics(first_state, 30, HarmonicValues)
# harmonic_df = pd.DataFrame(HarmonicValues)


