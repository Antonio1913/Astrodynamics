# This function will solve for the Parabolic Anomaly, B using trigonometric substitutions. The following method is
# derived using Thomas Barkers method or Barker's solution.
# delta_t  = change in time , sec
# p = semi-minor axis, km
# mass = mass of the object that is being orbited, kg
# G = gravitational constant
# B = Parabolic Anomaly

import math
from FUNCTIONS.ARCCOTANGENT import arccot


def kepeqtnP(delta_t, p, mass):
    G = 6.673*10E-20 # km^3/kg*s^2
    mu = G * mass
    np = 2 * math.sqrt(mu / p**3) # mean motion of the parabolic orbit
    s = 2 * arccot(3/2 * np * delta_t)
    w = math.atan(math.tan(s)**(1/3))
    B = 2 * (1 / math.tan(2 * w))
    print(mu)
    print(np)
    return math.radians(B)

delta_t = 53.7874 * 60
p = 25512 # km
mass = 5.9742 * 10E24
B = kepeqtnP(delta_t, p, mass)
print(f"Parabolic Anomaly B = {B}")

