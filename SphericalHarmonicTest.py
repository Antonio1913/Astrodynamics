import KEPLER_TOOLS as KT
from BODY_VALUES.MEAN_PLANETARY_CONSTANTS import Earth as E
from TOOLS.PLOTTING_TOOLS import orbitplot
from TOOLS.PROPAGATION_TOOLS import OrbitProp
import numpy as np
# Defining Unperturbed Orbit
pos_sat, vel_sat = KT.COE2RV(1000, 0.1, 40*(np.pi/180), 0, 0, 0, E.mu, "none")

state_sat = pos_sat.tolist() + vel_sat.tolist()
pos_mag = np.linalg.norm(pos_sat.T)
period_sat = 2 * np.pi * (np.sqrt(1000**3 / E.mu))
time_vec = np.linspace(0, period_sat, 10000)

pos_states, positions = OrbitProp(time_vec, state_sat, E.mu)
print(pos_sat)
orbitplot([positions], ['Satellite No-Perturbations'])


HarmonicValues = np.loadtxt('D:\ASTRODYNAMICS\EGM2008_Spherical_Harmonics\EGM2008')



