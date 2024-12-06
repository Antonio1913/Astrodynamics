import numpy as np
from typing import Tuple
from numpy.typing import NDArray
import KEPLER_TOOLS as KP
import scipy as sc
from MEAN_PLANETARY_CONSTANTS import Earth as E


# pos_states = initial position state
# Sat_state = should be a list 1
def OrbitProp(time_vec: [NDArray[np.float64]], Sat_state: [NDArray[np.float64]], mu: 'float'= E.mu) \
        -> Tuple[NDArray[np.float64], NDArray[np.float64]]:

    if np.size(Sat_state) < 5:
        np.transpose(Sat_state)

    # Designated size of time_vec as a variable
    N = len(time_vec)
    dt = time_vec[1] - time_vec[0]

    # Pre Allocating size of Sat_states
    Sat_states = np.zeros((6, N))
    Sat_states[:, 0] = Sat_state
    dt_counter = 0

    for i in range(1, N):
        r_vec_new, v_vec_new = KP.KEPLER(Sat_states[:3, i-1], Sat_states[3:, i-1], dt, mu)
        Sat_states[:, i] = np.concatenate((r_vec_new, v_vec_new))
        dt_counter = dt_counter + dt

    # Extract position vectors for plotting
    # Sat_states = np.array(Sat_states)
    pos_sat = Sat_states[:3, :]
    return Sat_states, pos_sat


# Spherical harmonic function will calculate the effect of Earths gravity at the inputted position vector.

#INPUT

#   pos_sat                     - [3xN]position vector of orbiting body
def sphericalharmonics(state_sat, desired_degree, Harmonic_values, mu=E.mu, r_bod=E.Radius):

    # Extracting pos_sat and vel_sat from state_sat
    pos_sat = state_sat[0:3, :]
    vel_sat = state_sat[3:6, :]

    # Transforming inputted position into spherical coordinates
    pos_norm = np.linalg.norm(pos_sat)
    sat_phi = np.asin(pos_sat[2, :] / pos_norm)
    sat_lambda = np.atan2(pos_sat[1, :], pos_sat[0, :])

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
    for Degree in range(2, desired_degree+1):

        # Order Array, [1xDegree+1]
        Order = np.arange(Degree + 2).reshape(-1, 1)

        # Normalization
        # Pre-allocating delta variable with 0 and 1 values, [Degree + 2,1], (+2 to account for m =0 and m = l+1)
        deltam0 = np.zeros([Degree + 2, 1])
        deltam0[0] = 1

        # EQUATION 10, Normalization Scaling Factor for Plm
        Plm_scaling = ((2 - deltam0) * ((2 * Degree) + 1) * (sc.special.factorial(Degree - Order) / sc.special.factorial(Degree + Order))) ** (1/2)

        # Assigning C and S values from the Normalized Harmonic Values, [1, Degree+1]
        values_end = sum(range(3, Degree + 2))
        values_beg = values_end - Degree

        C = np.array(Harmonic_values[values_beg-1:values_end, 2]).reshape(-1, 1)
        S = np.array(Harmonic_values[values_beg-1:values_end, 3]).reshape(-1, 1)

        # EQUATION 14 - Calculating Pl,l
        Pll = (sc.special.factorial((2 * Degree) - 1) / (2**(Degree - 1) * sc.special.factorial(Degree - 1))) * (xval**(Degree / 2))

        # Calculating Pl,l-1 (degrees - 2 stores values into first column, -3 stores value into third to last column
        Pllminus1 = xval2 * Pll

        # Pre_Allocating Array for all Order values for the Degree accounting for all position in trajectory.
        Plm_bar = np.zeros([Degree + 2, size_x])

        # Assigning Pll and Pllplus1 into Plm_bar
        Plm_bar[-3, :] = Plm_scaling[-3, :] * Pll
        Plm_bar[-2, :] = Plm_scaling[-2, :] * Pllminus1

        # Calculating sections of Equation 17 to make more readable
        calc1 = (Degree + Order + 1) * (Degree - Order)
        calc2 = (2 * (Order + 1)) * (1 / calc1)**(1/2) * xval2
        calc3 = (((Degree + Order + 2) * (Degree - Order - 1)) / calc1)**(1/2)

        # For Loop Calculates the remaining order values
        for order in range(Degree - 2, -1, -1):
            Plm_bar[order, :] = (calc2[order, :] * Plm_bar[order + 1, :]) - (calc3[order, :] * Plm_bar[order + 2, :])

        # Equation (8-25) du/dr
        # Calculating last section for dudr Equation (8-25)
        C_calc = C * np.cos(Order[0:Degree+1, :] * sat_lambda)
        S_calc = S * np.sin(Order[0:Degree+1, :] * sat_lambda)
        mtanphi = Order[0:Degree+1] * np.tan(sat_phi)

        # Calculating Equation (8-25) - *************** Ensure np.sum takes the sum up the columns to make a [1xN]
        dudr = np.sum(pos_ratio**Degree * (Degree + 1) * Plm_bar[0:Degree+1, :] * (C_calc + S_calc))

        # Equation (8-25) du/dphi
        # Creating Variable for Plm_plus1
        Plm_plus1scaling = Plm_scaling[1:, :]
        Plm_plus1 = Plm_bar[1:, :]

        # Corrected Plmplus1 by multiplying by lm normalization factor and divide by lm+1 normalization factor
        Plm_plus1corrected = Plm_plus1 * (Plm_scaling[0:-1, :] / Plm_plus1scaling)
        Plm_plus1corrected[-1, :] = 0

        # Calculating Equation (8-25) - *************** Ensure np.sum takes the sum up the columns to make a [1xN]
        dudphi = np.sum((pos_ratio**Degree) * (Plm_plus1corrected - (mtanphi * Plm_bar[0:Degree+1, :])) * (C_calc + S_calc))

        # Calculating Equation (8-25) - *************** Ensure np.sum takes the sum up the columns to make a [1xN]
        dudlambda = np.sum((pos_ratio**Degree) * (Order[0:Degree+1] * Plm_bar[0:Degree+1, :]) * (C_calc + S_calc))

        # Summation Calculations for All Degree and Order
        dudr_sum = dudr_sum + dudr
        dudphi_sum = dudphi_sum + dudphi
        dudlambda_sum = dudlambda_sum + dudlambda

    # Final Calculations for Equation (8-25)
    dUdR = (-mu / pos_norm**2) * dudr_sum
    dUdPhi = (mu / pos_norm) * dudphi_sum
    dUdLambda = (mu / pos_norm) * dudlambda_sum

    # Spherical Coordinates to Cartesian
    drdr = pos_sat / pos_norm
    dphidr1 = -pos_sat * pos_sat[2, :] / pos_norm**2
    dphidr1[2, :] = dphidr1[2, :] + 1
    xysum = pos_sat[0:1, :]**2
    dphidr = dphidr1 / (np.sqrt(np.sum(xysum)))
    dlambdadr = np.zeros(dphidr.shape)
    dlambdadr[0, :] = -pos_sat[1, :]
    dlambdadr[1, :] = pos_sat[0, :]
    dlambdadr = dlambdadr / (np.sum(xysum))

    # Acceleration components in te x, y, z directions
    gx = np.linalg.norm(dUdR * drdr)
    gy = np.linalg.norm(dUdPhi * dphidr)
    gz = np.linalg.norm(dUdLambda * dlambdadr)

    # Acceleration due to Spherical Harmonics
    a_spherharm = np.array([[gx], [gy], [gz]])

    # Velocity and Acceleration State
    accel_state = vel_sat.tolist() + a_spherharm.tolist()

    return accel_state


def Unperturbed_Orbit(pos_sat: np.array, vel_sat: np.array, a: 'float', mu: 'float'= E.mu, num_orbit =1) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    state_sat_init = np.concatenate((pos_sat, vel_sat)) # np.array(6,)
    period_sat = 2 * np.pi * (np.sqrt(a**3 / mu))
    num_time_int = 1000
    time_vec = np.linspace(0, period_sat * num_orbit, num_time_int)
    states, positions = OrbitProp(time_vec, state_sat_init, mu)
    return states, positions


def rkutta4(f, t, y, h):

    # Calculate one Runge-Kutta4 step
    k1 = f(t, y)
    k2 = f(t + 0.5 * h, y + 0.5 * k1 * h)
    k3 = f(t + 0.5 * h, y + 0.5 * k2 * h)
    k4 = f(t + h, y + k3 * h)

    sol = y + h / 6 * (k2 + 2 * k2 + 2 * k3 + k4)

    return sol