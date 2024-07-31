# This function solves C2(Chi) and C3(Chi) functions. These values are used to plug into the Kepler equations in term of
# the universal-variable.
# The following statements are ordered elliptical, parabolic, hyperbolic.

# INPUTS
#   psi- Universal Variable


# OUTPUTS
#   C2              - Universal Constant
#   C3              - Universal Constant

import math

def findc2c3(psi):
    if psi> 1e-6:
        C2 = (1 - math.cos(math.sqrt(psi))) / psi

        C3 = (math.sqrt(psi) - math.sin(math.sqrt(psi))) / math.sqrt(psi**3)
    elif psi< -1e-6:
        C2 = (1 - math.cosh(math.sqrt(-psi))) / psi

        C3 =(math.sinh(math.sqrt(-psi)) - math.isqrt(-psi)) / (math.sqrt(-psi**3))
    else:
        C2 = 1/2
        C3 = 1/6

    return C2, C3


# ---------------------------------------------------------------------------
psi_value = 0.0001
C2, C3 = findc2c3(psi_value)
print(f"C2 = {C2}, C3 = {C3}")