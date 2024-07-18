# This function solves C2(Chi) and C3(Chi) functions. These values are used to plug into the Kepler equations in term of
# the universal-variable.
# The following statements are ordered elliptical, parabolic, hyperbolic.
import math

def findc2c3(chi):
    if chi > 1e-6:
        C2 = (1 - math.cos(math.sqrt(chi))) / chi
        C3 = (math.sqrt(chi) - math.sin(math.sqrt(chi))) / math.sqrt(chi**3)
    elif chi < -1e-6:
        C2 = (1 - math.cosh(math.sqrt(-chi))) / chi
        C3 =(math.sinh(math.sqrt(-chi)) - math.isqrt(-chi)) / (math.sqrt(-chi**3))
    else:
        C2 = 1/2
        C3 = 1/6

    return C2, C3


# ---------------------------------------------------------------------------
Chi_value = 0.0001
C2, C3 = findc2c3(Chi_value)
# print(f"C2 = {C2}, C3 = {C3}")