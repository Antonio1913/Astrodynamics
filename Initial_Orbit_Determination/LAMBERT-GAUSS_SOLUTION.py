
import numpy as np
from BODY_VALUES.MEAN_PLANETARY_CONSTANTS import Earth as E

def lambert_Gauss(r0_vec, r_vec, delta_t, d_m):

# Defining position magnitudes
    r0 = np.linalg.norm(r0_vec)
    r = np.linalg.norm(r_vec)

    delta_v = np.acos(np.dot(r0_vec, r_vec) / (r0 * r))
    l = ((r0 + r) / (4 * np.sqrt(r0 * r) * np.cos(delta_v / 2))) - 0.5
    m = (E.Earth_mu * delta_t**2) / (2 * np.sqrt(r0 * r) * np.cos(delta_v / 2))**3

