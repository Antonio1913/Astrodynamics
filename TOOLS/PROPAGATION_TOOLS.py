import numpy as np
import KEPLER_TOOLS
import scipy as sc

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
        r_vec_new, v_vec_new = KEPLER_TOOLS.KEPLER(Sat_states[-1][:3], Sat_states[-1][3:6], changetime, mu)
        Sat_states.append(r_vec_new.tolist() + v_vec_new.tolist())

    # Extract position vectors for plotting
    Sat_states = np.array(Sat_states)
    pos_sat = Sat_states[:, :3]
    return Sat_states, pos_sat


# Spherical harmonic function will calculate the effect of Earths gravity at the inputted position vector.

#INPUT

#   pos_sat                     - position vector of orbiting body [3xN]
def sphericalharmonics(pos_sat, desired_degree, order, C, S):

    # Transforming inputted position into spherical coordinates
    pos_norm = np.linalg.norm(pos_sat)
    sat_lambda = np.asin(pos_sat[:, 3] / pos_norm)
    sat_phi = np.atan(pos_sat[:, 2], pos_sat[:, 1])

    # Calculating x value
    x = np.sin(sat_phi)
    xval = 1 - (x**2)

    # Pre-Allocating Arrays for Variables
    Plm = np.zeros([desired_degree + 2, desired_degree-1])
    Pbarlm = np.zeros([desired_degree + 2, desired_degree-1])

    #Beginning for Loops to Calculate Recursions
    for degrees in range(2, 1, desired_degree+1):

        # Degree correction to store values in correct column
        store = degrees - 2

        # Calculating Pl,l (degrees - 2 stores values into first column, -2 stores value into second to last column)
        Plm[-2, store] = (sc.special.factorial((2 * degrees) - 1) / (2**(degrees - 1) * sc.special.factorial(degrees - 1))) * (xval**(degrees / 2))

        # Calculating Pl,l-1 (degrees - 2 stores values into first column, -3 stores value into third to last column
        Plm[-3, store] = (x / xval**(1/2)) * Plm[degrees, -1]

        # Normalization
        # Pre-allocating delta variable with 0 and 1 values
        deltam0 = np.zeros([degrees, 1])
        deltam0[1] = 1

        # (degrees - 2 starts for loop at Pl,l-2 since Pl,l and Pl,l-1 have already been calculated)
        for order in range(degrees - 2, -1, -1):
            calc1 = 1 / ((degrees - order) * (degrees + order + 1))
            # Calculating Pl,m
            Plm[order, store] = calc1 * (((2 * (order + 1)) * (x / xval**(1/2)) * Plm[order+1, store]) - Plm[order+2, store])

        orders = np.arange(degrees+1)

        #Normalization
        Pbarlm = ((2 - deltam0) * ((2 * degrees) + 1) * ((sc.special.factorial(degrees - order)) / (sc.special.factorial(degrees + order))))**(1/2) * Plm

        # Normalizing Geo potential Coefficients
        cbar = (1 / ((2 - deltam0) * ((2 * degrees) + 1))) * ((sc.special.factorial(degrees + order)) / (sc.special.factorial(degrees - order)))**(1/2) * C
        sbar = (1 / ((2 - deltam0) * ((2 * degrees) + 1))) * ((sc.special.factorial(degrees + order)) / (sc.special.factorial(degrees - order)))**(1/2) * S

        return a_vec

