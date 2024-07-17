# This function provides a solution for Kepler's problem using classical orbital elements. This method does not include
# any perturbations. Therefore, the only value to change in the two positions is the true anomaly or the special orbit
# parameters.

# This method uses the known position and velocity to find the initial position's anomaly. Then determining the Mean anomaly.
# Finally using the mean anomaly to determine the true anomaly at the new position all information is present to
# determine the new position and velocity.

# INPUTS
#   r0_vec      - initial position
#   v0_vec      - initial velocity
#   delta_t     - time frame of observation

# OUTPUTS
#   r_vec       - position vector at the observation time
#   v_vec       - velocity at new position


import numpy as np
 from Rv2COE import RV2COE
 from TrueAnomaly2Anomaly import nutoAnomaly
 from KepEqtnElliptical import kepEqtnE
 from KepEqtnParabolic import kepeqtnP
 from KepEqtnHyperbolic import kepeqtnH
# from BODY_VALUES.MEAN_PLANETARY_CONSTANTS import Earth as E

def KeplerCOE (r0_vec, v0_vec, delta_t, mu):


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
        anomaly_type = "Eccentric"
        arg = (anomaly,)

#   Parabolic Case
    if ecc == 1:
        h_vec = np.linalg.cross(r0_vec, v0_vec)
        h_mag = np.linalg.norm(h_vec)
        p = h_mag**2 / mu
        M0 = anomaly0 + (anomaly0**3 / 3)
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
        true_anomaly_new = anomaly2nu(ecc, anomaly_type, *arg)
    else:
    arg_latitude = anomaly

      r_vec_IJK, v_vec_IJK = COE2RV(p, ecc, incl, ascending_node, arg_perigee, true_anomaly, mu, *args)













