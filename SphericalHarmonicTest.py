import KEPLER_TOOLS as KT
from BODY_VALUES.MEAN_PLANETARY_CONSTANTS import Earth as E
from TOOLS.PLOTTING_TOOLS import orbitplot
from TOOLS.PROPAGATION_TOOLS import OrbitProp, rkutta4, sphericalharmonics
import numpy as np
# Defining Unperturbed Orbit
pos_sat, vel_sat = KT.COE2RV(1000, 0.1, 40*(np.pi/180), 0, 0, 0, E.mu, "none")


state_sat_init = (pos_sat.tolist() + vel_sat.tolist())
pos_mag = np.linalg.norm(pos_sat.T)
period_sat = 2 * np.pi * (np.sqrt(1000**3 / E.mu))
num_orbit = 10
num_time_int = 10000
time_vec = np.linspace(0, period_sat * num_orbit, num_time_int)
dt = (period_sat * num_orbit) / np.size(time_vec)

pos_states, positions = OrbitProp(time_vec, state_sat_init, E.mu)
print(pos_sat)
# orbitplot([positions], ['Satellite No-Perturbations'])
state_sat = np.zeros((6, num_time_int))
state_sat[:, 0] = np.array(state_sat_init).flatten()

HarmonicValues = np.loadtxt('D:\ASTRODYNAMICS\EGM2008_Spherical_Harmonics\EGM2008')

def f_sphericalharmonics(y, 5, HarmonicValues):
    return sphericalharmonics(t, y, HarmonicValues)


for step in range(num_time_int - 1):
    # state_sat[step + 1] = rkutta4(sphericalharmonics(state_sat[:, step].reshape(6, 1), 5, HarmonicValues), time_vec[step], state_sat[:, step].reshape(6, 1), dt)

    state_sat[step + 1] = rkutta4(f_sphericalharmonics, time_vec[step], state_sat[:, step].reshape(6, 1), dt)


