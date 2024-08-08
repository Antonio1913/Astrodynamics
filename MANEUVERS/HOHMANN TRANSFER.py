
import numpy as np
from BODY_VALUES.MEAN_PLANETARY_CONSTANTS import Earth as E
from TOOLS.PROP


def HohmannTransfer(r_inital, r_final, mu):

    # Semi Major Axis of the Transfer Orbit
    a_trans = (r_inital + r_final) / 2

    # Initial Velocity and Initial Velocity Transfer
    v_int = np.sqrt(E.Earth_mu / r_inital)
    v_transa = np.sqrt(((2 * E.Earth_mu) / r_inital) - (E.Earth_mu / a_trans))

    # Final Velcoity and Final Velocity Transfer
    v_final = np.sqrt(E.Earth_mu / r_final)
    v_transb = np.sqrt(((2 * E.Earth_mu) / r_final) - (E.Earth_mu / a_trans))

    # Delta V for Initial Manuever
    delta_va = v_transa - v_int

    # Delta V for Final Manuever to Maintain targeted Altitude
    delta_vb = v_final - v_transb

    # Total Delta V
    delta_v = abs(delta_va) + abs(delta_vb)

    # Transfer Time of Orbit
    tau_trans = np.pi * np.sqrt(a_trans**3 / E.Earth_mu)

    return a_trans, tau_trans, delta_va, delta_vb


# r_int = np.array([[6500], [0], [0]])
# r_mag = np.linalg.norm(r_int)
# v_mag =
a = 6500
ecc = 0
incl = 0
ascending_node = 0
arg_perigee = 0
true_anomaly = 0
mu = E.Earth_mu


r_vec_IJK, v_vec_IJK = COE2RV(a, ecc, incl, ascending_node, arg_perigee, true_anomaly, mu, "lambda_true", 0)
