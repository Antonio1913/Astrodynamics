

import numpy as np

def FINDTOF(r0_vec, r_vec, p, mu):

    r0_mag = np.linalg.norm(r0_vec)#Finds the magnitude of the original position
    r_mag = np.linalg.norm(r_vec)#Finds the magnitude of the second position
    delta_true_anomaly = np.acos(np.dot(r0_vec, r_vec) / (r0_mag * r_mag)) #Finds the change in true anoaly at the two positions
    k = r0_mag * r_mag * (1 - np.cos(delta_true_anomaly))
    l = r0_mag + r_mag
    m = r0_mag * r_mag * (1 + np.cos(delta_true_anomaly))
    val1 = (((2 * m) - l**2) * p**2)
    val2 = (2 * k * l * p) - k**2
    a = (m * k * p) / (val1 + val2)
    f = 1 - ((r_vec / p) * (1 - np.cos(delta_true_anomaly)))
    g = (r0_mag * r_vec * np.sin(delta_true_anomaly)) / (np.sqrt(mu * p))

    #Elliptical Orbit case where a > 0
    if a > 0:






