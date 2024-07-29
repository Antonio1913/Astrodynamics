# THIS FUNCTION CONVERTS THE KNOWN ANOMALY TO OUTPUT THE TRUE ANOMALY.
# THE ORDER OF THE INPUTS IS VERY IMPORTANT

# INPUTS
#   ecc                 - Eccentricity of Orbit
#   anomaly_type        - "Eccentric" or "Parabolic" or "Hyperbolic"
#   delta_t             - time frame of observation, sec

# SPECIAL CASE INPUT
# If Parabolic *args should be "B" , "p", "r"
#   p                   - semi-minor axis, km
#   r                   - position magnitude, km

# OUTPUTS
#   nu2                 - True Anomaly, rad


import math as m

def anomaly2nu(ecc, anomaly_type, *arg):

    if anomaly_type == "Eccentric":
        if len(arg) != 1:
            raise ValueError("For Eccentric Anomaly, only one argument (E) is required or input is misspelled.")
        E = arg[0]
        # nu1 = m.asin((m.sin(E) * m.sqrt(1 - ecc**2)) / (1 - (ecc * m.cos(E))))
        nu2 = m.acos((m.cos(E[0]) - ecc) / (1 - (ecc * m.cos(E[0]))))


    elif anomaly_type == "Parabolic":
        if len(arg) != 3:
             raise ValueError("For Parabolic Anomaly, three arguments (B, p, r) are required.")
        B, p, r = arg
        # nu1 = m.asin(p * B / r)
        nu2 = m.acos((p - r) / r)

    elif anomaly_type == "Hyperbolic":
        if len(arg) != 1:
            raise ValueError("For Hyperbolic Anomaly, only one argument (H) is required.")
        H = arg[0]
        # nu1 = m.asin((-m.sinh(H) * m.sqrt(ecc ** 2 - 1)) / (1 - (ecc * m.cosh(H))))
        nu2 = m.acos((m.cosh(H[0]) - ecc) / (1 - (ecc * m.cosh(H[0]))))

    else:
        raise ValueError("Invalid anomaly_type provided.")

    return nu2




# ecc = 0.008101039975211971
# anomaly_type = "Eccentric"
# arg = 2.491107053821039
# nu2 = anomaly2nu(ecc, anomaly_type, arg)
# print(nu2)


