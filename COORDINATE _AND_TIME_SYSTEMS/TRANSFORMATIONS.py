import numpy as np
from BODY_VALUES.MEAN_PLANETARY_CONSTANTS import Earth as E

# THIS FUNCTION TRANSFORMS THE INPUTTED POSITION VECTOR IN THE EARTH CENTERED EACH FIXED FRAME INTO THE
# GEODETIC HEIGHT, LONGITUDE, AND ELLIPSOIDAL HEIGHT. ALSO CONTAINS AN ADJUSTMENT FUNCTION WHEN NEAR THE POLES. Utilizes
# iterations to converge for phi_gd.

#INPUTS
#   r_vec_ecef                   - [1x3] Initial Position

# OUTPUTS
#   phi_gd                      - Geodetic Height, rad
#   longitude                   - Longitude, rad
#   h_ellp                      - Ellipsoidal Height, km

def ECEFtoLATLON(r_vec_ecef):
    r_mag_ecef = np.linalg.norm(r_vec_ecef)
    r_I = r_vec_ecef[0]
    r_J = r_vec_ecef[1]
    r_k = r_vec_ecef[2]

    r_delta_sat = np.sqrt(r_I**2 + r_J**2)
    longitude = np.asin(r_J / r_delta_sat)
    delta = np.asin(r_vec_ecef[2] / r_mag_ecef)

    #Inital Guess
    phi_gd_old = delta

    #Iteration Parameters
    max_iterations = 100
    tolerance = 1 * 10 ** -8

    for i in range(max_iterations):
        C = E.Earth_Radius / np.sqrt(1 - (E.Earth_eccentricity**2 * np.sin(phi_gd_old**2)))
        phi_gd_it = np.atan((r_k + (C* E.Earth_eccentricity**2 * np.sin(phi_gd_old))) / r_delta_sat)

        if abs(phi_gd_it - phi_gd_old) < tolerance:
            phi_gd = phi_gd_it
            break
        phi_gd_old = phi_gd_it

    h_ellp = (r_delta_sat / np.cos(phi_gd)) - C

    if phi_gd * (180 / np.pi) > 89 or phi_gd * (180 / np.pi) < -89:
        S = C * (1 - E.Earth_eccentricity**2)
        h_ellp = (r_k / np.sin(phi_gd)) - S

    return phi_gd, longitude, h_ellp



