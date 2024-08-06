


# INPUTS
#   M           - Mean Anomaly, rad
#   ecc         - Eccentricity

# OUTPUTS
#   H           - Hyperbolic Anomaly, rad

import math
from TOOLS.FUNCTIONS import sign

def kepeqtnH(M, e):
    if e < 1.6:
        if -math.pi < M < 0 or M > math.pi:
            Hn = M - e
        else:
            Hn = M + e
    else:
        if e <3.6 and abs(M) > math.pi:
            Hn = M - sign(M) * e
        else:
            Hn = M / (e - 1)

    tolerance = 1 * 10 ** -8  # The standard tolerance for Newton_Raphson Method
    max_iterations = 100  # Maximizes the number of iterations to 20
    for i in range(max_iterations):
        Hnplus1 = Hn + ((M - e * math.sinh(Hn) + Hn) / (e * math.cosh(Hn) - 1))
        if abs(Hnplus1 - Hn) < tolerance:
            H = Hnplus1
            break
        Hn = Hnplus1
    print(f"Warning: Tolerance not met after {max_iterations} iterations")



#Test Case
M = 235.4 * math.pi / 180
e = 2.4
H = kepeqtnH(M,e)
# print(f"Hyperbolic Anomaly H = {H}")



