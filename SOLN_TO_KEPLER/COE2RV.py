#function explanation and limitation if explained

# INPUTS
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

# OUTPUTS
#   r_vec       - position vector at the observation time, km
#   v_vec       - velocity at new position, km/s


import numpy as np
from TOOLS.FUNCTIONS import Rot1, Rot3
from BODY_VALUES.MEAN_PLANETARY_CONSTANTS import Earth as E

def COE2RV (a, ecc, incl, ascending_node, arg_perigee, true_anomaly, mu, *args):
    if len(args) > 2:
        raise ValueError(f"Too many inputs, ensure the only args inputted correspond to the type of orbit")

        #Setting Conditional Terms
        # Circular and Equatorial
    if args[0] == "lambda_true" and ecc == 0 and incl == 0:
            arg_perigee = 0
            ascending_node = 0
            lambda_true = args[1]
            true_anomaly = lambda_true

#       Circular and Inclined
    elif args[0] == "arg_latitude" and ecc == 0 and incl < 0:
            arg_perigee = 0
            arg_latitude = args[1]
            true_anomaly = arg_latitude

#       Elliptical  and Equatorial
    elif args[0] == "arg_perigee_true" and ecc < 0 and incl == 0:
            ascending_node = 0
            arg_perigee_true = args[1]
            arg_perigee = arg_perigee_true
    elif args[0] == "none":
        pass

    else:
        raise ValueError(f"Unexpected value for args {args}")

#   Vector array of the position of the body in the PQW axis
    r_PQW = np.array([[(p * (np.cos(true_anomaly))) / (1 + (ecc * np.cos(true_anomaly)))],
                      [(p * np.sin(true_anomaly)) / (1 + (ecc * np.cos(true_anomaly)))],
                          [0]])

    v_PQW = np.array([[-np.sqrt(mu / p) * np.sin(true_anomaly)],
                        [np.sqrt(mu / p) * (ecc + np.cos(true_anomaly))],
                        [0]])
#   Rotation operation in order to get the vectors in the geocentric equatorial system
    Rotations = np.dot(Rot3(-ascending_node), np.dot(Rot1(-incl), Rot3(-arg_perigee)))

#   Position vector in the IJK reference frame
    r_vec_IJK = np.dot(Rotations, r_PQW)

#   Velocity vector in the IJK reference frame
    v_vec_IJK = np.dot(Rotations, v_PQW)

    return r_vec_IJK, v_vec_IJK


# p = 11067.790
# e = 0
# i = 0
# Omega = 227.89
# omega = 53.38
# nu = 92.335
# mu = E.Earth_mu
# case = "lambda_true"
#
# r_vec_IJK, v_vec_IJK = COE2RV(p, e, i, Omega, omega, nu, mu, case, 60)
# print("position vector in the ijk reference frame", r_vec_IJK, "km")
# print("velocity vector in the ijk reference frame", v_vec_IJK, "km/s")

