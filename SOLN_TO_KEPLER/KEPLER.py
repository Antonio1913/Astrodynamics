# This function provides a solution for Kepler's problem using universal variables and one formulation for all orbit
# types. For this reason this method has an advantage over KeplerCOE other than having many more operations.
# Newton-Raphson iteration is used to determine when the universal variable converges.
#
# This method does not include any perturbations. Therefore, the only value to change in the two positions is the
# true anomaly or the special orbit parameters.

# INPUTS
#   r0_vec      - initial position
#   v0_vec      - initial velocity
#   delta_t     - time frame of observation
#   mu          - Gravitational Constant, km^3/s^2

# OUTPUTS
#   r_vec       - position vector at the observation time
#   v_vec       - velocity at new position


import numpy as np
from FUNCTIONS.FUNCTIONS import arccot, sign
from FIndC2C3 import findc2c3

def KEPLER(r0_vec, v0_vec, delta_t, mu):

#   Calculating the magnitude of the position and velocity
    r0_mag = np.linalg.norm(r0_vec)
    v0_mag = np.linalg.norm(v0_vec)

#   Saved value to decrease the number of calculations
    sqrtmu = np.sqrt(mu)
    val1 = np.dot(r0_vec, v0_vec) / sqrtmu
    val2 = sqrtmu * delta_t

#   orbits specific mechanical energy
    energy = (v0_mag**2 / 2) - (mu / r0_mag)
    a = -(mu / (2 * energy))
    alpha = 1 / a

#   Checking for orbital types

#   circular or elliptical
    if alpha > 0.000001:
        chi_n = sqrtmu * delta_t * alpha

#   parabolic
    elif abs(alpha) < 0.000001:
#   Calculating the vector and magnitude of the angular momentum
        h_vec = np.linalg.cross(r0_vec, v0_vec)
        h_mag = np.linalg.norm(h_vec)
        p = h_mag**2 / mu # Semi-Minor Axis
        s = arccot(3 * (np.sqrt(mu / p**3) * delta_t)) / 2
        w = (np.atan(np.tan(s)**(1/3)))
        chi_n = np.sqrt(p) * 2 * (1 / np.tan(2 * w))

#   hyperbolic
    elif alpha < -0.000001:
        a = 1 / alpha
        chi_n = sign(delta_t) * np.sqrt(-a) * np.log((-2 * mu * alpha * delta_t) / (np.dot(r0_vec, v0_vec) + (sign(delta_t) * np.sqrt(-mu * a) * (1 - (r0_mag * alpha)))))

#   Iteration settings
    tolerance = 1 * 10**-6
    max_iterations = 100

#   Newton-Raphson iteration
    for i in range(max_iterations):
        psi = chi_n**2 * alpha
        C2, C3 = findc2c3(psi)
        r = (chi_n**2 * C2) + (val1 * chi_n * (1 - (psi * C3))) + (r0_mag * (1 - (psi * C2)))
        chi_n_plus1 = chi_n + ((val2 - (chi_n**3 * C3) - (val1 * chi_n**2 * C2) - (r0_mag * chi_n * (1 - (psi * C3)))) / r)

        if abs(chi_n - chi_n_plus1) < tolerance:
            chi_n = chi_n_plus1
            break
        chi_n = chi_n_plus1

#   f and g function calculations
    f = 1 - ((chi_n**2 / r0_mag) * C2)
    g = delta_t - ((chi_n**3 / sqrtmu) * C3)
    g_dot = 1 - ((chi_n**2 / r) * C2)
    f_dot = (sqrtmu / (r * r0_mag)) * chi_n * ((psi * C3) - 1)

#   New Position Vector
    r_vec = (f * r0_vec) + (g * v0_vec)

#   New Velocity Vector
    v_vec = (f_dot * r0_vec) + (g_dot * v0_vec)

    check = (f * g_dot) - (f_dot * g)
    if check > 1 - tolerance and check < 1 + tolerance:
        return r_vec, v_vec
    else:
        raise ValueError("f and g function check did not equal 1, Check input value, units, and format")





#test case
r_ijk = np.array([1131.340, -2282.343, 6672.423])  # Position vector in kilometers
v_ijk = np.array([-5.64305, 4.30333, 2.42879])   # Velocity vector in kilometers per second
mu = 398600.4415  # Standard gravitational parameter for Earth in km^3/s^2
delta_t = 40 * 60 # sec

r_vec, v_vec = KEPLER(r_ijk, v_ijk, delta_t, mu)

print(r_vec)
print(v_vec)




