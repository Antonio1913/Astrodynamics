# THIS FUNCTION CONVERTS THE KNOWN ANOMALY TO OUTPUT THE TRUE ANOMALY.
# THE ORDER OF THE INPUTS IS VERY IMPORTANT
# e - Eccentricity
# anomaly_type = "Eccentric" or "Parabolic" or "Hyperbolic"
# If Eccentric *args should be "E"
# If Parabolic *args should be "B" , "p", "r"
#   p - semi-minor axis
#   r - position
# If Hyperbolic *args should be "H"
import math as m

def anomaly2nu(e, anomaly_type, *args):
    if anomaly_type == "Eccentric":
        if len(args) !=1:
            raise ValueError("For Eccentric Anomaly, only one argument (E) is required.")
        E = args[0]
        nu1 = m.asin((m.sin(E) * m.sqrt(1 - e**2)) / (1 - (e * m.cos(E))))
        nu2 = m.acos((m.cos(E) - e) / (1 - (e * m.cos(E))))
    
    elif anomaly_type == "Parabolic":  
        if len(args) != 3:
            raise ValueError("For Parabolic Anomaly, three arguments (B, p, r) are required.")
        B, p, r = args
        nu1 = m.asin(p * B / r)
        nu2 = m.acos((p - r) / r)

    elif anomaly_type == "Hyperbolic":
        if len(args) != 1:
            raise ValueError("For Hyperbolic Anomaly, only one argument (H) is required.")
        H = args[0]
        nu1 = m.asin((-m.sinh(H) * m.sqrt(e ** 2 - 1)) / (1 - (e * m.cosh(H))))
        nu2 = m.acos((m.cosh(H) - e) / (1 - (e * m.cosh(H))))

    else:
        raise ValueError("Invalid anomaly_type provided.")
    return nu1, nu2
