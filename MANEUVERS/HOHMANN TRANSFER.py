
import numpy as np
from BODY_VALUES.MEAN_PLANETARY_CONSTANTS import Earth as E
from TOOLS.PROPAGATION_TOOLS import OrbitProp
from TOOLS.PLOTTING_TOOLS import orbitplot
from SOLN_TO_KEPLER.COE2RV import COE2RV
from SOLN_TO_KEPLER.Rv2COE import RV2COE


def HohmannTransfer(r_initial, r_final, mu):

    # Semi Major Axis of the Transfer Orbit
    a_trans = (r_initial + r_final) / 2

    # Initial Velocity and Initial Velocity Transfer
    v_int = np.sqrt(E.mu / r_initial)
    v_transa = np.sqrt(((2 * E.mu) / r_initial) - (E.mu / a_trans))

    # Final Velcoity and Final Velocity Transfer
    v_final = np.sqrt(E.mu / r_final)
    v_transb = np.sqrt(((2 * E.mu) / r_final) - (E.mu / a_trans))

    # Delta V for Initial Manuever
    delta_va = v_transa - v_int

    # Delta V for Final Manuever to Maintain targeted Altitude
    delta_vb = v_final - v_transb

    # Total Delta V
    delta_v = abs(delta_va) + abs(delta_vb)

    # Transfer Time of Orbit
    tau_trans = np.pi * np.sqrt(a_trans**3 / E.mu)

    return a_trans, tau_trans, delta_va, delta_vb


def PlotHohmannTransfer(a, lambdatrue, r_new, mu):
    r_vec_IJK, v_vec_IJK = COE2RV(a, 0, 0, 0, 0, 0, mu, "lambda_true", lambdatrue)
    pos_sat = r_vec_IJK.tolist() + v_vec_IJK.tolist()
    r_mag = np.linalg.norm(r_vec_IJK)
    period_sat = 2 * np.pi * (np.sqrt(r_mag**3 / mu))
    time_vec = np.linspace(0, period_sat, 1000)

    pos_states, positions = OrbitProp(time_vec, pos_sat, mu)
    print(positions)

    a1, ecc1, incl, ascending_node, arg_perigee, true_anomaly, arg_perigee_true, arg_latitude, lambda_true = RV2COE(pos_states[0, 0:3], pos_states[0, 3:6], mu)

    # orbitplot(positions, 'satellite')

    r_initial = r_mag
    r_final = r_mag + r_new

    a_trans, tau_trans, delta_va, delta_vb = HohmannTransfer(r_initial, r_final, mu)

    pos_sat1 = pos_sat.copy()
    pos_sat1[4] += delta_va
    period_sat1 = 2 * np.pi * np.sqrt(a_trans**3 / mu)
    time_vec1 = np.linspace(0, period_sat1, 1000)
    pos_states1, positions1 = OrbitProp(time_vec1, pos_sat1, mu)

    positions1_mag = np.linalg.norm(positions1, axis=1)
    max_pos1 = np.argmax(positions1_mag)
    positions1_transfer = positions1[0:max_pos1, :]

    # Second Burn To Maintain Transfer Orbit
    pos_state2 = pos_states1[max_pos1, :]
    pos_state2[4] -= delta_vb
    pos_state2_mag = np.linalg.norm(pos_state2[0:2, :])
    period_sat2 = 2 * np.pi * (np.sqrt(pos_state2_mag**3 / mu))
    time_vec2 = np.linspace(0, period_sat2, 1000)

    pos_states2, positions2 = OrbitProp(time_vec2, pos_state2, mu)


    orbitplot([positions, positions1_transfer, positions2], ['initial orbit', 'transfer', 'new orbit'])

    return orbitplot



# a = 6500
# lambdatrue = 0
# mu = E.mu
# r_new = 5000
# PlotHohmannTransfer(a, lambdatrue, r_new, mu)

