# this function will calculate the orbital elements of a orbiting body given a known position vector and its velocity
# at that position.
# Inputs must be a [1X3] and in km
# Outputs are in Km and degrees
import math

# INPUTS
#   r0_vec      - [1x3] initial position
#   v0_vec      - [1x3] initial velocity
#   mu                  - Gravitational Constant, km^3/s^2

# OUTPUTS
#   p                   - Semi-Parameter, km
#   ecc                 - Eccentricity
#   incl                - Inclination, rad
#   ascending_node      - Ascending Node, rad
#   arg_perigee         - Argument of Perigee, rad
#   true_anomaly        - True Anomaly, rad
#   mu                  - Gravitational Constant, km^3/s^2
#   *args
#   lambda_true         - Lambda True (Circular Equatorial Orbit)
#   arg_latitude        - Argument of Latitude (Circular Inclined Orbit)
#   arg_perigee_true    - True Argument of Perigee (Equatorial Elliptical Orbit)

import numpy as np

def RV2COE (r_ijk, v_ijk, mu):
    # Angular Momentum
    h_vec = np.cross(r_ijk, v_ijk) # km^2/s
    h_mag = np.linalg.norm(h_vec) # km^2/s

    # Node Vector
    k_vec = np.array([0, 0, 1])
    n_vec = np.cross(k_vec, h_vec) # km^2/s
    n_mag = np.linalg.norm(n_vec) # km^2/s

    # Eccentricity
    pos_mag = np.linalg.norm(r_ijk) # km
    vel_mag = np.linalg.norm(v_ijk) # km/s
    ecc_vec = (((vel_mag**2 - (mu / pos_mag)) * r_ijk) - (np.dot(r_ijk, v_ijk) * v_ijk)) / mu
    ecc = np.linalg.norm(ecc_vec)

    # Specific Orbital Energy
    energy = vel_mag**2 / 2 - (mu / pos_mag) # km^2/s^2

    # Semi-Major Axis and Semi-Minor Axis
    if ecc != 1:
        a = -(mu / (2 * energy)) # km
        p = a * (1 - ecc**2) # km

    else:
        p = h_mag**2 / mu # km
        a = float('inf') # km

    #Inclination
    incl = (np.acos(h_vec[2] / h_mag))

    #Ascending Node
    ascending_node = (np.acos(n_vec[0] / n_mag))
    if n_vec[1] < 0:
        ascending_node = (2 * math.pi) - ascending_node

    # Argument of Perigee
    arg_perigee = (np.acos(np.dot(n_vec, ecc_vec) / (n_mag * ecc)))
    if ecc_vec[2] < 0:
        arg_perigee = (2 * math.pi) - arg_perigee

    # True Anomaly
    true_anomaly = (np.acos(np.dot(ecc_vec, r_ijk) / (ecc * pos_mag)))
    if np.dot(r_ijk, v_ijk) < 0:
        true_anomaly = (2 * math.pi) - true_anomaly

    #Special Cases
    # Elliptical and Equatorial - True Argument of Perigee
    arg_perigee_true = (np.acos(ecc_vec[0] / ecc))
    if ecc_vec[1] < 0:
        arg_perigee_true = (2 * math.pi) - arg_perigee_true

    # Circular Inclined - Argument of Latitude
    arg_latitude = (np.acos(np.dot(n_vec, r_ijk) / (n_mag * pos_mag)))
    if r_ijk[2] < 0:
        arg_latitude = (2 * math.pi) - arg_latitude

    # Circular Equatorial - True Lambda
    lambda_true  = (np.acos(r_ijk[0] / pos_mag))
    if r_ijk[1] < 0:
        lambda_true = (2 * math.pi) - lambda_true

    return a, p, ecc, incl, ascending_node, arg_perigee, true_anomaly, arg_perigee_true, arg_latitude, lambda_true








# r_ijk = np.array([6524.834, 6862.875, 6448.296])  # Position vector in kilometers
# v_ijk = np.array([4.901327, 5.533756, -1.976341])   # Velocity vector in kilometers per second
# mu = 398600.4418  # Standard gravitational parameter for Earth in km^3/s^2
#
# # Call the function
# a, p, ecc, incl, ascending_node, arg_perigee, true_anomaly, arg_perigee_true, arg_latitude, lambda_true = RV2COE(r_ijk, v_ijk, mu)
#
# # Print the results
# print("Semi-major axis (a):", a, "km")
# print("Semi-latus rectum (p):", p, "km")
# print("Inclination (incl):", incl, "degrees")
# print("Right Ascension of Ascending Node (RAAN):", ascending_node, "degrees")
# print("Argument of Perigee (arg_perigee):", arg_perigee, "degrees")
# print("True Anomaly (true_anomaly):", true_anomaly, "degrees")
# print("arg_perigee_true (arg_perigee_true):", arg_perigee_true, "degrees")
# print("Argument of Latitude (arg_latitude):", arg_latitude, "degrees")
# print("True Longitude (lambda_true):", lambda_true, "degrees")
#




