# Assumptions must be close enough to the true solution so there is no violation of the linear assumption
# of the Newton-Raphson method
# This function will provide the Eccentric Anomaly given the Mean Anomaly and the eccentricity.

import math

def kepEqtnE(M,e):
    if -math.pi < M < 0 or M > math.pi:
        En = M - e
    else:
        En = M + e

    tolerance = 1 * 10**-8    # The standard tolerance for Newton_Raphson Method
    max_iterations = 100     # Maximizes the number of iterations to 20

    for i in range(max_iterations):
        Enplus1 = En + ((M - En + (e * math.sin(En)))/(1 - e * math.cos(En)))
        if abs(Enplus1 - En) < tolerance:
            E = Enplus1
            return E
        En = Enplus1
    print(f"Warning: Tolerance not met after {max_iterations} iterations")



M =  4.10850506
e = 0.4
E = kepEqtnE(M,e)
print(f"Eccentric Anomaly E = {E}")






