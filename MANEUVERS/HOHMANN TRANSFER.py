
import numpy as np
from BODY_VALUES.MEAN_PLANETARY_CONSTANTS import Earth as E
from TOOLS.PROPAGATION_TOOLS import OrbitProp
from TOOLS.PLOTTING_TOOLS import orbitplot
from SOLN_TO_KEPLER.COE2RV import COE2RV


def HohmannTransfer(r_initial, r_final, mu):

    # Semi Major Axis of the Transfer Orbit
    a_trans = (r_initial + r_final) / 2

    # Initial Velocity and Initial Velocity Transfer
    v_int = np.sqrt(E.Earth_mu / r_initial)
    v_transa = np.sqrt(((2 * E.Earth_mu) / r_initial) - (E.Earth_mu / a_trans))

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
pos_sat = r_vec_IJK.tolist() + v_vec_IJK.tolist()
r_mag = np.linalg.norm(r_vec_IJK)
period_sat = 2 * np.pi * (np.sqrt(r_mag**3 / mu))
time_vec = np.linspace(0, period_sat, 1000)

pos_states, positions = OrbitProp(time_vec, pos_sat, mu)
print(positions)


# orbitplot(positions, 'satellite')

r_initial = r_mag
r_final = r_mag + 6000

a_trans, tau_trans, delta_va, delta_vb = HohmannTransfer(r_initial, r_final, mu)

pos_sat1 = pos_sat.copy()
pos_sat1[4] += delta_va
period_sat1 = 2 * np.pi * np.sqrt(a_trans**3 / mu)
time_vec1 = np.linspace(0, period_sat1, 1000)
pos_states1, positions1 = OrbitProp(time_vec1, pos_sat1, mu)
# orbitplot(positions1, 'satellite')

orbitplot([positions, positions1], ['initial orbit', 'transfer'])

