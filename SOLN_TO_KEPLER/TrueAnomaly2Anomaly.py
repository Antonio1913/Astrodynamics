# THIS FUNCTION WILL OUTPUT THE ECCENTRIC, PARABOLIC, AND HYPERBOLIC ANOMALY DEPENDING ON THE TYPE OF THE INPUT ORBIT.
# THIS FUNCTION NEEDS THE ECCENTRICITY AND THE TRUE ANOMALY TO GIVE THE ANOMALY NEEDED TO DERIVE THE FLIGHT PATH ANGLE.

# INPUTS
# e             - Eccentricity
# nu            - True Anomaly, rad

# OUTPUTS
# E             - Eccentric Anomaly
# P             - Parabolic Anomaly
# H             - Hyperbolic Anomaly

import numpy as np

def nutoAnomaly (e, nu):
    if e < 1.0:
        #E1 = np.asin((math.sin(nu) * np.sqrt(1 - e**2)) / (1 + (e * np.cos(nu))))
        E = np.acos((e + np.cos(nu)) / (1 + (e * np.cos(nu))))
        return E
    elif e == 1:
        B = np.tan(nu / 2)
        return  B
    else:
        #H1 = math.asinh((math.sin(nu) * math.sqrt(e**2 - 1)) / (1 + (e * math.cos(nu))))
        H = np.acosh((e + np.cos(nu)) / (1 + (e * np.cos(nu))))
        return H
