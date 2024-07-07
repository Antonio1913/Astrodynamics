

import numpy as np

def RV2COE (r_ijk, v_ijk):
    h_vec = np.cross(r_ijk, v_ijk)
    k_vec = np.array([[0], [0], [1]])
    n_vec = np.cross(k_vec, h_vec)
    vel_mag =
    ecc_vec =
