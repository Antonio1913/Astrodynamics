import numpy as np

from SOLN_TO_KEPLER.KEPLER import KEPLER


# pos_states = initial position state
# Sat_state = should be a list 1
def OrbitProp(time_vec, Sat_state, mu):

    if np.size(Sat_state) < 5:
        np.transpose(Sat_state)

    #Designated size of time_vec as a variable
    N = len(time_vec)

    Sat_states = [Sat_state]

    for i in range(1, N):
        changetime = time_vec[i] - time_vec[i - 1]
        r_vec_new, v_vec_new = KEPLER(Sat_states[-1][:3], Sat_states[-1][3:6], changetime, mu)
        Sat_states.append(r_vec_new.tolist() + v_vec_new.tolist())

    # Extract position vectors for plotting
    Sat_states = np.array(Sat_states)
    pos_sat = Sat_states[:, :3]
    return Sat_states, pos_sat