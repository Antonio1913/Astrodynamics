import numpy as np
from .BODY_CONSTANTS import Earth as E


# THIS FUNCTION TRANSFORMS THE INPUTTED POSITION VECTOR IN THE EARTH CENTERED EACH FIXED FRAME INTO THE
# GEODETIC HEIGHT, LONGITUDE, AND ELLIPSOIDAL HEIGHT. ALSO CONTAINS AN ADJUSTMENT FUNCTION WHEN NEAR THE POLES. Utilizes
# iterations to converge for phi_gd.

# INPUTS
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

    r_delta_sat = np.sqrt(r_I ** 2 + r_J ** 2)
    longitude = np.asin(r_J / r_delta_sat)
    delta = np.asin(r_vec_ecef[2] / r_mag_ecef)

    # Inital Guess
    phi_gd_old = delta

    # Iteration Parameters
    max_iterations = 100
    tolerance = 1 * 10 ** -8

    for i in range(max_iterations):
        C = E.Radius / np.sqrt(1 - (E.eccentricity ** 2 * np.sin(phi_gd_old ** 2)))
        phi_gd_it = np.atan((r_k + (C * E.eccentricity ** 2 * np.sin(phi_gd_old))) / r_delta_sat)

        if abs(phi_gd_it - phi_gd_old) < tolerance:
            phi_gd = phi_gd_it
            break
        phi_gd_old = phi_gd_it

    h_ellp = (r_delta_sat / np.cos(phi_gd)) - C

    if phi_gd * (180 / np.pi) > 89 or phi_gd * (180 / np.pi) < -89:
        S = C * (1 - E.eccentricity ** 2)
        h_ellp = (r_k / np.sin(phi_gd)) - S

    return phi_gd, longitude, h_ellp


#
class TimeConversions:
    solarday2siderealday = 1.002737909350795
    siderealday2solarday = 0.997269566329084
    solarday2siderealtime = r"24h 3m 56.5553678s"
    siderealday2solartime = r"23h 56m 4.090524s"


# This function takes the inputted date ranging from March 1, 1900, to February 28, 2100. If the day consisted of a leap
# second use 61 seconds. The year must consist of four digits.

# INPUT
#   year                - 4 digit year ranging from 1900 to 2100
#   month               - Starting from January as 1 to December as 12
#   hour                - Using Military time or 24 hour timescale
#   min                 - 60 min timescale
#   sec                 - If day consists of leap second then use out of 61

# OUTPUT
#   JD                  - Julian Date in, days

def JULIANDATE(year, month, day, hour, minute, sec):
    JD = (367 * year) - np.trunc((7 * (year + np.trunc((month + 9) / 12))) / 4) + np.trunc(
        (275 * month) / 9) + day + 1721013.5 + (((((sec / 60) + minute) / 60) + hour) / 24)
    return JD


# Rotation Matrices
# Rot 1 finds the rotation about the X axis
def Rot1(alpha):
    matrix = np.array([[1, 0, 0],
                       [0, np.cos(alpha), np.sin(alpha)],
                       [0, -np.sin(alpha), np.cos(alpha)]])
    return matrix


# Rot 2 finds the rotation about the Y axis
def Rot2(alpha):
    matrix = np.array([[np.cos(alpha), 0, -np.sin(alpha)],
                       [0, 1, 0],
                       [np.sin(alpha), 0, np.cos(alpha)]])
    return matrix


# Rot 3 finds the rotation about the Z axis
def Rot3(alpha):
    matrix = np.array([[np.cos(alpha), np.sin(alpha), 0],
                       [-np.sin(alpha), np.cos(alpha), 0],
                       [0, 0, 1]])
    return matrix
