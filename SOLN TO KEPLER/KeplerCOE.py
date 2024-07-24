# This function provides a solution for Kepler's problem using classical orbital elements. This method does not include
# any perturbations. Therefore, the only value to change in the two positions is the true anomaly or the special orbit
# parameters.

# This method uses the known position and velocity to find the initial position's anomaly. Then determining the
# Mean anomaly. Finally using the mean anomaly to determine the true anomaly at the new position all information is
# present to determine the new position and velocity.

# INPUTS
#   r0_vec      - [1x3] initial position
#   v0_vec      - [1x3] initial velocity
#   delta_t     - time frame of observation

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
from Rv2COE import RV2COE
from TrueAnomaly2Anomaly import nutoAnomaly
from KepEqtnElliptical import kepEqtnE
from KepEqtnParabolic import kepeqtnP
from KepEqtnHyperbolic import kepeqtnH
from Anomaly2nu import anomaly2nu
from COE2RV import COE2RV
# from BODY_VALUES.MEAN_PLANETARY_CONSTANTS import Earth as E

def KeplerCOE(r0_vec, v0_vec, delta_t, mu):

    a, p, ecc, incl, ascending_node, arg_perigee, true_anomaly, arg_perigee_true, arg_latitude, lambda_true = RV2COE(r0_vec, v0_vec, mu)

    if ecc != 0:
        anomaly0 = nutoAnomaly(ecc, true_anomaly)
    else:
        anomaly0 = arg_latitude

#   Eccentric Case
    if ecc < 1:
        M0 = anomaly0 - (ecc * np.sin(anomaly0)) # Mean Anomaly
#       Mean Motion
        n = np.sqrt(mu / a**3)
        M = M0 + (n * delta_t)
        anomaly = kepEqtnE(M, ecc) #Eccentric Anomaly
        anomaly = anomaly
        anomaly_type = "Eccentric"
        arg = (anomaly,)

#   Parabolic Case
    if ecc == 1:
        h_vec = np.linalg.cross(r0_vec, v0_vec)
        h_mag = np.linalg.norm(h_vec)
        p = h_mag**2 / mu
        #M0 = anomaly0 + (anomaly0**3 / 3)
        anomaly = kepeqtnP(delta_t, p, mu) # Parabolic anomaly
        r_mag = p / 2 * (1 + anomaly**2)
        anomaly_type = "Parabolic"
        arg = (anomaly, p, r_mag)


#   Hyperbolic Case
    if ecc > 1:
        M0 = ecc * np.sinh(anomaly0) - anomaly0
#       Mean Motion
        n = np.sqrt(mu / a ** 3)
        M = M0 + (n * delta_t)
        anomaly = kepeqtnH(M , ecc) # Hyperbolic Anomaly
        anomaly_type = "Hyperbolic"
        arg = (anomaly,)


    if ecc != 0:
        # arg is E or H for Eccentric or Hyperbolic case but B,p, r for parabolic
        true_anomaly_new = anomaly2nu(ecc, anomaly_type, arg)
        arg2 = "none"
    else:
        arg_latitude = anomaly
        arg2 = arg_latitude

    r_vec_IJK, v_vec_IJK = COE2RV(a, ecc, incl, ascending_node, arg_perigee, true_anomaly_new, mu, arg2)

    return r_vec_IJK, v_vec_IJK



r_ijk = np.array([1131.340, -2282.343, 6672.43])  # Position vector in kilometers
v_ijk = np.array([-5.64305, 4.30333, 2.42879])   # Velocity vector in kilometers per second
mu = 398600.4418  # Standard gravitational parameter for Earth in km^3/s^2
delta_t = 40 * 60 # sec

r_vec, v_vec = KeplerCOE(r_ijk, v_ijk, delta_t, mu)

print(r_vec)
print(v_vec)









