#gg

import math

def findC2C3(Chi):
    if Chi > 1*10**-6:
        C2 = (1-math.cos(math.isqrt(Chi)))/Chi
        C3 = (math.isqrt(Chi) - math.sin(math.isqrt(Chi))) / math.isqrt(Chi**3)
    elif Chi < -1*10**-6:
        C2 = (1 - math.cosh(math.isqrt(-Chi))
        C3 =(math.sinh(math.isqrt(-Chi)) - math.isqrt(-Chi)) / (math.isqrt(-Chi**3))
    else:
        C2 = 1/2
        C3 = 1/6
print("C2 = {C2}, C3 = {C3}")