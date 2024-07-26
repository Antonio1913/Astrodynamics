

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

#   Calculating the Central Angle between Departing Point P1 and Arriving Point P2
#   Theta is Dependent on whether the input was for the ascending stage (AA) or descending stage (DA)
    adjustment = arg

    if adjustment == "AA":
        theta = np.mod(np.asin(np.sin(phi_G / np.sin(incl)) - u_1), 2 * np.pi)

    elif adjustment == "DA":
        theta = np.mod(np.pi - (np.asin(np.sin(phi_G / np.sin(incl)) - u_1)), 2 * np.pi)
    else:
        print(f"Incorrect Orbit Adjustment. Only valid inputs, 'DA', 'AA'.")




