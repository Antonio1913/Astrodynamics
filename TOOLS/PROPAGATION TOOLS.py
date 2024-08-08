from SOLN_TO_KEPLER.KEPLER import KEPLER

def OrbitProp(time_vec, pos_states, mu):

    for i in range(1, len(time_vec)):
        changetime = time_vec[i] - time_vec[i - 1]
        r_vec_new, v_vec_new = KEPLER(pos_states[-1][:3], pos_states[-1][3:6], changetime, mu)
        pos_states.append(r_vec_new.tolist() + v_vec_new.tolist())

    # Extract position vectors for plotting
    pos_states = np.array(pos_states)
    positions = pos_states[:, :3]
    return pos_states, positions