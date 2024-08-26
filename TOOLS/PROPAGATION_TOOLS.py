import numpy as np
import KEPLER_TOOLS
import scipy as sc
from BODY_VALUES.MEAN_PLANETARY_CONSTANTS import Earth as E

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

#   pos_sat                     - [3xN]position vector of orbiting body
def sphericalharmonics(pos_sat, desired_degree, Harmonic_values, mu = E.mu, r_bod = E.Radius):

    # Transforming inputted position into spherical coordinates
    pos_norm = np.linalg.norm(pos_sat)
    sat_phi = np.asin(pos_sat[:, 3] / pos_norm)
    sat_lambda = np.atan2(pos_sat[:, 2], pos_sat[:, 1])

    # Ratio between Bodies Radius and position
    pos_ratio = -r_bod / pos_norm

    # Calculating x value, [1xN]
    x = np.sin(sat_phi)
    xval = 1 - (x**2)
    xval2 = x / xval**(1/2)
    size_x = np.size(x)

    # Pre-Allocating Summation Arrays
    dudr_sum = 0
    dudphi_sum = 0
    dudlambda_sum = 0

    #Beginning for Loops to Calculate Recursions
    for Degree in range(2, 1, desired_degree+1):

        # Order Array, [1xDegree+1]
        Order = np.arange(Degree + 2).reshape(-1, 1)

        # Normalization
        # Pre-allocating delta variable with 0 and 1 values, [Degree + 2,1], (+2 to account for m =0 and m = l+1)
        deltam0 = np.zeros([Degree + 2, 1])
        deltam0[1] = 1

        # EQUATION 10, Normalization Scaling Factor for Plm
        Plm_scaling = ((2 - deltam0) * ((2 * Degree) + 1) * (sc.special.factorial(Degree - Order) / sc.special.factorial(Degree + Order))) ** (1/2)

        # Assigning C and S values from the Normalized Harmonic Values, [1, Degree+1]
        values_beg = np.sum(np.array([2, Degree + 1]))
        values_end = values_beg + Degree + 1

        C = np.array(Harmonic_values[values_beg:values_end, 2]).reshape(-1, 1)
        S = np.array(Harmonic_values[values_beg:values_end, 3]).reshape(-1, 1)

        # EQUATION 14 - Calculating Pl,l
        Pll = (sc.special.factorial((2 * Degree) - 1) / (2**(Degree - 1) * sc.special.factorial(Degree - 1))) * (xval**(Degree / 2))

        # Calculating Pl,l-1 (degrees - 2 stores values into first column, -3 stores value into third to last column
        Pllplus1 = xval2 * Pll

        # Pre_Allocating Array for all Order values for the Degree accounting for all position in trajectory.
        Plm_bar = np.zeros([Degree + 2, size_x])

        # Assigning Pll and Pllplus1 into Plm_bar
        Plm_bar[-1, :] = Plm_scaling[-1, :] * Pllplus1
        Plm_bar[-2, :] = Plm_scaling[-2, :] * Pll

        # Calculating sections of Equation 17 to make more readable
        calc1 = (Degree + Order + 1) * (Degree - Order)
        calc2 = (2 * (Order + 1)) * (1 / calc1)**(1/2) * xval2
        calc3 = (((Degree + Order + 2) * (Degree - Order - 1)) / calc1)**(1/2)

        # For Loop Calculates the remaining order values
        for order in range(Degree - 2, -1, -1):
            Plm_bar[order, :] = (calc2[order, :] * Plm_bar[order + 1, :]) - (calc3[order+2, :] * Plm_bar[order + 2, :])

        # Equation (8-25) du/dr
        # Calculating last section for dudr Equation (8-25)
        C_calc = C * np.cos(Order[0:Degree, :] * sat_lambda)
        S_calc = S * np.sin(Order[0:Degree, :] * sat_lambda)
        mtanphi = Order * np.tan(sat_phi)

        # Calculating Equation (8-25) - *************** Ensure np.sum takes the sum up the columns to make a [1xN]
        dudr = np.sum(pos_ratio**Degree * (Degree + 1) * Plm_bar[0:Degree, :] * (C_calc + S_calc))

        # Equation (8-25) du/dphi
        # Creating Variable for Plm_plus1
        Plm_plus1scaling = Plm_scaling[1:, :]
        Plm_plus1 = Plm_bar[1:, :]
        # Corrected Plmplus1 by multiplying by lm normalization factor and divide by lm+1 normalization factor
        Plm_plus1corrected = Plm_plus1 * (Plm_scaling[0:-1,:] / Plm_plus1scaling)

        # Calculating Equation (8-25) - *************** Ensure np.sum takes the sum up the columns to make a [1xN]
        dudphi = (pos_ratio**Degree) * (Plm_plus1corrected - (mtanphi *







        # Summation Calculations for All Degree and Order
        dudr_sum = dudr_sum + dudr
        dudphi_sum = dudphi_sum + dudphi
        dudlambda_sum = dudlambda_sum + dudlambda







    # Final Calculations for Equation (8-25)
    dUdR = (-mu / pos_norm**2) * dudr_sum
    dUdPhi = (mu / pos_norm) * dudphi_sum
    dUdLambda = (mu / pos_norm) * dudlambda_sum

    return a_vec

