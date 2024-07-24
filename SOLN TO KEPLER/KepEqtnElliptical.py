# Assumptions must be close enough to the true solution so there is no violation of the linear assumption
# of the Newton-Raphson method
# This function will provide the Eccentric Anomaly given the Mean Anomaly and the eccentricity.
# INPUTS
#   M           - Mean Anomaly, rad
#   ecc         - Eccentricity

# OUTPUTS
#   E           - Eccentric Anomaly, deg



import math

def kepEqtnE(M,ecc):
    if -math.pi < M < 0 or M > math.pi:
        En = M - ecc
    else:
        En = M + ecc

    tolerance = 1 * 10**-8    # The standard tolerance for Newton_Raphson Method
    max_iterations = 100     # Maximizes the number of iterations to 20

    for i in range(max_iterations):
        Enplus1 = En + ((M - En + (ecc * math.sin(En)))/(1 - ecc * math.cos(En)))
        if abs(Enplus1 - En) < tolerance:
            E = Enplus1 * (180 / math.pi)
            return E
        En = Enplus1
    print(f"Warning: Tolerance not met after {max_iterations} iterations")




# Test Case
# M =  4.10850506
# e = 0.4
# E = kepEqtnE(M,e)
# print(f"Eccentric Anomaly E = {E}")






