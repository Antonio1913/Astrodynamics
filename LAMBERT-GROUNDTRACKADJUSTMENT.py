

# INPUTS
#   t1                          - Departing Time, sec
#   a_star                      - Semi-Major Axis Before Maneuver, km
#   ecc_star                    - Eccentricity Before Maneuver
#   incl_star                   - Inclination Before Maneuver, rad
#   ascending_node_star         - Right Ascension of the Ascending Node Before Maneuver, rad
#   arg_perigee_star            - Argument of Perigee Before Maneuver km
#   f_star                      - True Anomaly Before Maneuver, rad
#   L_G                         - Geocentric Longitude of Ground Target or Intersection Point, rad
#   phi_G                       - Geocentric Latitude of Ground Target or Intersection Point, rad
#   r_min                       - Minimum Position Vector, km
#   r_max                       - Maximum Position Vector, km
#   D                           - Parameter (integer) of Time of Flight
#   incl                        -
#   H_2                         - Altitude at Time 2, km
#   AA                          -
#   kappa                       - Transfer Direction
#   K                           - Mode Parameter

# OUTPUTS


import numpy as np
from BODY_VALUES.MEAN_PLANETARY_CONSTANTS import Earth as E


def ExtendedLambert(t1, a_star, ecc_star, incl_star, ascending_node_star, arg_perigee_star, true_anomaly_star, L_G, phi_G,
                     r_min, r_max, D, incl, H_2, kappa, K, *arg):


#   Argument of Latitude of the Spacecrafts Initial Orbit at T1
    u_star = arg_perigee_star + true_anomaly_star

#   Phi_1 = Phi_star - Geocentric Latitude of the Departing Point P1
    phi_star = np.asin(np.sin(incl_star) * np.sin(u_star))

#   Right Ascension of the Ascending Node at Time 1 at the time of Impulsive Maneuver
    ascending_node_1 = (ascending_node_star + (np.asin(np.tan(phi_star) / np.tan(incl_star))) -
                        (kappa * (np.asin(np.tan(phi_star) / np.tan(incl)))) + ((1 - kappa) * (np.pi / 2)))

#   Unit Vector in the Direction of the Ascending Node at Time 1
    i_vec_Omega1 = np.array([np.cos(ascending_node_1), np.sin(ascending_node_1), 0])

#   Unit Vector in the Direction of Position Vector r1 or rstar
    i_vec_r1 = np.array([(np.cos(ascending_node_star) * np.cos(u_star)) - (np.sin(ascending_node_star) * np.sin(u_star) * np.cos(incl_star)),
                         (np.sin(ascending_node_star) * np.cos(u_star)) + (np.cos(ascending_node_star) * np.sin(u_star) * np.cos(incl_star)),
                         np.sin(u_star) * np.sin(incl_star)])

#   Argument of Latitude at Time 1
    u_1 = np.acos(np.dot(i_vec_Omega1, i_vec_r1))

#   Calculating the Central Angle between Departing Point P1 and Arriving Point P2 and thre longitudinal difference
#   between target G and the intersection G', delta_l
#   Theta is Dependent on whether the input was for the ascending stage (AA) or descending stage (DA)

    # Converting time in UTC into Greenwich Mean Sidereal Time(GMST)

    Theta_gmst_t1 =

    if arg == "AA":
        theta = np.mod(np.asin(np.sin(phi_G / np.sin(incl)) - u_1), 2 * np.pi)
        delta_l = np.mod(ascending_node_1 + (np.asin(np.tan(phi_G) / np.tan(incl))) - L_G - Theta_gmst_t1, 2 * np.pi) + 2 * np.pi * D

    elif arg == "DA":
        theta = np.mod(np.pi - (np.asin(np.sin(phi_G / np.sin(incl)) - u_1)), 2 * np.pi)
        delta_l = np.mod(ascending_node_1 + (np.pi - (np.asin(np.tan(phi_G) / np.tan(incl)))) - L_G - Theta_gmst_t1, 2 * np.pi) + 2 * np.pi * D

    else:
        print(f"Incorrect Orbit Adjustment. Only valid inputs, 'DA', 'AA'.")

#   Defining the two position vectors
    r2 = E.Earth_Radius + H_2
    lambda_max = np.sqrt((r1 + r2 - abs(r1 - r2)) / (r1 + r2 + abs(r1 - r2)))

#Determining Keplerian Starters, lambda, x
# Step 1 - Setting a0 = a_star, ecc0 = ecc_star
    a0 = a_star
    ecc0 = ecc_star

# Step 2 - Solving Equations 5 and 9 to compute TOF
#   Calculating the secular drift rate of the right ascension of the ascending node caused by J2 Perturbation
    B_Omega = 1.5 * E.Earth_J2 * E.Earth_Radius**2 * np.sqrt(E.Earth_mu) * np.cos(incl)
    B_e = 1 - ecc0**2

    Omega_J2 = B_Omega * a0**(-7 / 2) * B_e**-2

#   Secular drift rate of the argument of perigee caused by J2 Perturbation
    B_omega = -1.5 * E.Earth_J2 * E.Earth_Radius**2 * np.sqrt(E.Earth_mu) * ((2.5 * np.sin(incl**2)) - 2)
    omeegadot_J2 = B_omega * a0**(-7 / 2) * B_e**-2

#   Initial Guess for TOF
    tau_0 = delta_l / (E.Earth_Rotation - Omega_J2)

#   Rotation Angle Between t1 and t2 of the Flyby Orbit
    theta_omega = omeegadot_J2 * tau_0

#   Auxiliary Angle
    theta_C_hat = theta - theta_omega

#Ensuring input is within range
    low_range_k = -theta_C_hat / (2 * np.pi)
    upper_range_k = 1 - (theta_C_hat / (2 * np.pi))
    if K < low_range_k or K > upper_range_k:
        print(f"K input is not an ideal value with other input parameters.")

#   Initial Guess for Transfer angle
    theta_C0 = np.mod(theta_C_hat, 2 * np.pi)

#   Replacing Theta_C with theta_C0 in Equations 15 21 26 to compute initial guesses for Semi perimeter and the parameter lambda
    c = np.sqrt((r1 + r1)**2 - 4 * r1 *r2 * np.cos((theta_C0 / 2)**2))
    s0 = (r1 + r2 + c) / 2
    lambda0 = (np.sqrt(r1 * r2) * np.cos(theta_C0 / 2)) / s0 + 7

















