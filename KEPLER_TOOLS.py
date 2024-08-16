# THIS FILE STORES THE CODE THAT PROVIDES THE BASES TO PROPAGATE ORBITS.

import numpy as np
from TOOLS.FUNCTIONS import Rot1, Rot3, sign, arccot
from TOOLS.STRUCTURE_TOOLS import ensure_numpy_array

#Universly Defined Values
tolerance = 1 * 10**-8    # The standard tolerance for Newton_Raphson Method
max_iterations = 100     # Maximizes the number of iterations to 20

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

def anomaly2nu(ecc, anomaly_type, *arg):

        if anomaly_type == "Eccentric":
            if len(arg) != 1:
                raise ValueError("For Eccentric Anomaly, only one argument (E) is required or input is misspelled.")
            E = arg[0]
            # nu1 = m.asin((m.sin(E) * m.sqrt(1 - ecc**2)) / (1 - (ecc * m.cos(E))))
            nu2 = np.acos((np.cos(E[0]) - ecc) / (1 - (ecc * np.cos(E[0]))))


        elif anomaly_type == "Parabolic":
            if len(arg) != 3:
                raise ValueError("For Parabolic Anomaly, three arguments (B, p, r) are required.")
            B, p, r = arg
            # nu1 = np.asin(p * B / r)
            nu2 = np.acos((p - r) / r)

        elif anomaly_type == "Hyperbolic":
            if len(arg) != 1:
                raise ValueError("For Hyperbolic Anomaly, only one argument (H) is required.")
            H = arg[0]
            # nu1 = np.asin((-np.sinh(H) * np.sqrt(ecc ** 2 - 1)) / (1 - (ecc * np.cosh(H))))
            nu2 = np.acos((np.cosh(H[0]) - ecc) / (1 - (ecc * np.cosh(H[0]))))

        else:
            raise ValueError("Invalid anomaly_type provided.")

        return nu2

################################################################################################
# INPUTS
#   p                   - Semi-Parameter, km
#   ecc                 - Eccentricity
#   incl                - Inclination, rad
#   ascending_node      - Ascending Node, rad
#   arg_perigee         - Argument of Perigee, rad
#   true_anomaly        - True Anomaly, rad
#   mu                  - Gravitational Constant, km^3/s^2
#   *args
#   lambda_true         - Lambda True (Circular Equatorial Orbit)
#   arg_latitude        - Argument of Latitude (Circular Inclined Orbit)
#   arg_perigee_true    - True Argument of Perigee (Equatorial Elliptical Orbit)

# OUTPUTS
#   r_vec       - position vector at the observation time, km
#   v_vec       - velocity at new position, km/s

def COE2RV (a, ecc, incl, ascending_node, arg_perigee, true_anomaly, mu, *args):
    if len(args) > 2:
        raise ValueError(f"Too many inputs, ensure the only args inputted correspond to the type of orbit")

        #Setting Conditional Terms
        # Circular and Equatorial
    if args[0] == "lambda_true" and ecc == 0 and incl == 0:
            arg_perigee = 0
            ascending_node = 0
            lambda_true = args[1]
            true_anomaly = lambda_true

#       Circular and Inclined
    elif args[0] == "arg_latitude" and ecc == 0 and incl < 0:
            arg_perigee = 0
            arg_latitude = args[1]
            true_anomaly = arg_latitude

#       Elliptical  and Equatorial
    elif args[0] == "arg_perigee_true" and ecc < 0 and incl == 0:
            ascending_node = 0
            arg_perigee_true = args[1]
            arg_perigee = arg_perigee_true
    elif args[0] == "none":
        pass

    else:
        raise ValueError(f"Unexpected value for args {args}")

#   Vector array of the position of the body in the PQW axis
    r_PQW = np.array([[(a * (np.cos(true_anomaly))) / (1 + (ecc * np.cos(true_anomaly)))],
                      [(a * np.sin(true_anomaly)) / (1 + (ecc * np.cos(true_anomaly)))],
                          [0]])

    v_PQW = np.array([[-np.sqrt(mu / a) * np.sin(true_anomaly)],
                        [np.sqrt(mu / a) * (ecc + np.cos(true_anomaly))],
                        [0]])
#   Rotation operation in order to get the vectors in the geocentric equatorial system
    Rotations = np.dot(Rot3(-ascending_node), np.dot(Rot1(-incl), Rot3(-arg_perigee)))

#   Position vector in the IJK reference frame
    r_vec_IJK = np.dot(Rotations, r_PQW)

#   Velocity vector in the IJK reference frame
    v_vec_IJK = np.dot(Rotations, v_PQW)

    return r_vec_IJK, v_vec_IJK

################################################################################################
# This function solves C2(Chi) and C3(Chi) functions. These values are used to plug into the Kepler equations in term of
# the universal-variable.
# The following statements are ordered elliptical, parabolic, hyperbolic.

# INPUTS
#   psi- Universal Variable


# OUTPUTS
#   C2              - Universal Constant
#   C3              - Universal Constant
def findc2c3(psi):
    if psi> 1e-6:
        C2 = (1 - np.cos(np.sqrt(psi))) / psi

        C3 = (np.sqrt(psi) - np.sin(np.sqrt(psi))) / np.sqrt(psi**3)
    elif psi< -1e-6:
        C2 = (1 - np.cosh(np.sqrt(-psi))) / psi

        C3 =(np.sinh(np.sqrt(-psi)) - np.sqrt(-psi)) / (np.sqrt(-psi**3))
    else:
        C2 = 1/2
        C3 = 1/6

    return C2, C3

################################################################################################
# INPUTS
#   r0_vec      - [1x3] initial position
#   v0_vec      - [1x3] initial velocity
#   p           - time frame of observation
#   mu          - Gravitational Constant, km^3/s^2

# OUTPUTS
#   TOF         - Time of Flight, sec

def FINDTOF(r0_vec, r_vec, p, mu):

    r0_mag = np.linalg.norm(r0_vec)#Finds the magnitude of the original position
    r_mag = np.linalg.norm(r_vec)#Finds the magnitude of the second position
    delta_true_anomaly = np.acos(np.dot(r0_vec, r_vec) / (r0_mag * r_mag)) #Finds the change in true anoaly at the two positions
    k = r0_mag * r_mag * (1 - np.cos(delta_true_anomaly))
    l = r0_mag + r_mag
    m = r0_mag * r_mag * (1 + np.cos(delta_true_anomaly))
    val1 = (((2 * m) - l**2) * p**2)
    val2 = (2 * k * l * p) - k**2
    a = (m * k * p) / (val1 + val2)
    f = 1 - ((r_vec / p) * (1 - np.cos(delta_true_anomaly)))
    g = (r0_mag * r_vec * np.sin(delta_true_anomaly)) / (np.sqrt(mu * p))

    #Elliptical Orbit case where a > 0
    if a > 0:
        f_dot = np.sqrt(mu/p) * np.tan(delta_true_anomaly / 2) * (((1 - np.cos(delta_true_anomaly)) / p) - (1 / r0_mag) - (1 / r_mag))
        delta_E = np.acos(1 - ((r0_mag / a) * (1 - f)))
        sin_delta_E = -r0_mag * r_mag * f_dot / (np.sqrt(mu * a))
        TOF = g + (np.sqrt(a**3 / mu) * (delta_E - sin_delta_E))

    elif a == float('inf'):
        c = np.sqrt(r0_mag**2 + r_mag**2 - (2 * r0_mag * r_mag * np.cos(delta_true_anomaly)))
        s = (r0_mag + r_mag + c) / 2
        TOF = 2/3 * np.sqrt(s**3 / 2 * mu) * (1 - ((s - c) / s)**(3/2))


    elif a < 0:
        delta_H = np.acosh(1 + ((f - 1)) * (r0_mag / a))
        TOF = g + (np.sqrt((-a)**3 / mu) * (np.sinh(delta_H) - delta_H))

    return TOF

################################################################################################
# Assumptions must be close enough to the true solution so there is no violation of the linear assumption
# of the Newton-Raphson method
# This function will provide the Eccentric Anomaly given the Mean Anomaly and the eccentricity.

# INPUTS
#   M           - Mean Anomaly, rad
#   ecc         - Eccentricity

# OUTPUTS
#   E           - Eccentric Anomaly, rad

def kepEqtnE(M,ecc):
    if -np.pi < M < 0 or M > np.pi:
        En = M - ecc
    else:
        En = M + ecc

    for i in range(max_iterations):
        Enplus1 = En + ((M - En + (ecc * np.sin(En)))/(1 - ecc * np.cos(En)))
        if abs(Enplus1 - En) < tolerance:
            E = Enplus1
            return E
        En = Enplus1
    print(f"Warning: Tolerance not met after {max_iterations} iterations")

################################################################################################
# INPUTS
#   M           - Mean Anomaly, rad
#   ecc         - Eccentricity

# OUTPUTS
#   H           - Hyperbolic Anomaly, rad

def kepeqtnH(M, e):
    if e < 1.6:
        if -np.pi < M < 0 or M > np.pi:
            Hn = M - e
        else:
            Hn = M + e
    else:
        if e <3.6 and abs(M) > np.pi:
            Hn = M - sign(M) * e
        else:
            Hn = M / (e - 1)

    for i in range(max_iterations):
        Hnplus1 = Hn + ((M - e * np.sinh(Hn) + Hn) / (e * np.cosh(Hn) - 1))
        if abs(Hnplus1 - Hn) < tolerance:
            H = Hnplus1
            break
        Hn = Hnplus1
    print(f"Warning: Tolerance not met after {max_iterations} iterations")

################################################################################################
# This function will solve for the Parabolic Anomaly, B using trigonometric substitutions. The following method is
# derived using Thomas Barkers method or Barker's solution.

# INPUTS
#   delta_t      - change in time , sec
#   p            - semi-minor axis, km
#   mu           - gravitational constant body mass, km^3 / s^2

# OUTPUTS
#   B            - Parabolic Anomaly , rad

def kepeqtnP(delta_t, p, mu):
    n = 2 * np.sqrt(mu / p**3) # mean motion of the parabolic orbit
    s = (1/2 * arccot(3/2 * n * delta_t))
    w = (np.atan(np.tan(s)**(1/3)))
    B = (2 * (1 / np.tan(2 * w)))

    return B

################################################################################################
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

def KEPLER(r0_vec, v0_vec, delta_t, mu):

    # Ensures Correct Variable Structure
    r0_vec = ensure_numpy_array(r0_vec)
    v0_vec = ensure_numpy_array(v0_vec)


#   Calculating the magnitude of the position and velocity
    r0_mag = np.linalg.norm(r0_vec)
    v0_mag = np.linalg.norm(v0_vec)

#   Saved value to decrease the number of calculations
    sqrtmu = np.sqrt(mu)
    val1 = np.dot(r0_vec.T, v0_vec) / sqrtmu
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

################################################################################################
# This function provides a solution for Kepler's problem using classical orbital elements. This method does not include
# any perturbations. Therefore, the only value to change in the two positions is the true anomaly or the special orbit
# parameters.

# This method uses the known position and velocity to find the initial position's anomaly. Then determining the
# Mean anomaly. Finally using the mean anomaly to determine the true anomaly at the new position all information is
# present to determine the new position and velocity.

# INPUTS
#   r0_vec      - [1x3] initial position
#   v0_vec      - [1x3] initial velocity
#   delta_t     - time frame of observation

# OUTPUTS
#   p                   - Semi-Parameter, km
#   ecc                 - Eccentricity
#   incl                - Inclination, rad
#   ascending_node      - Ascending Node, rad
#   arg_perigee         - Argument of Perigee, rad
#   true_anomaly        - True Anomaly, rad
#   mu                  - Gravitational Constant, km^3/s^2
#   *args
#   lambda_true         - Lambda True (Circular Equatorial Orbit)
#   arg_latitude        - Argument of Latitude (Circular Inclined Orbit)
#   arg_perigee_true    - True Argument of Perigee (Equatorial Elliptical Orbit)

def KeplerCOE(r0_vec, v0_vec, delta_t, mu):

# All angles are in radians
    a, p, ecc, incl, ascending_node, arg_perigee, true_anomaly, arg_perigee_true, arg_latitude, lambda_true = RV2COE(r0_vec, v0_vec, mu)

    if ecc != 0:
        anomaly0 = nutoAnomaly(ecc, true_anomaly) #radians
    else:
        anomaly0 = arg_latitude

#   Eccentric Case
    if ecc < 1:
        M0 = anomaly0 - (ecc * np.sin(anomaly0)) # Mean Anomaly
#       Mean Motion
        n = np.sqrt(mu / a**3)
        M = M0 + (n * delta_t)
        anomaly = kepEqtnE(M, ecc) #Eccentric Anomaly
        anomaly = anomaly
        anomaly_type = "Eccentric"
        arg = (anomaly,)

#   Parabolic Case
    if ecc == 1:
        h_vec = np.linalg.cross(r0_vec, v0_vec)
        h_mag = np.linalg.norm(h_vec)
        p = h_mag**2 / mu
        #M0 = anomaly0 + (anomaly0**3 / 3)
        anomaly = kepeqtnP(delta_t, p, mu) # Parabolic anomaly
        r_mag = p / 2 * (1 + anomaly**2)
        anomaly_type = "Parabolic"
        arg = (anomaly, p, r_mag)


#   Hyperbolic Case
    if ecc > 1:
        M0 = ecc * np.sinh(anomaly0) - anomaly0
#       Mean Motion
        n = np.sqrt(mu / a ** 3)
        M = M0 + (n * delta_t)
        anomaly = kepeqtnH(M , ecc) # Hyperbolic Anomaly
        anomaly_type = "Hyperbolic"
        arg = (anomaly,)


    if ecc != 0:
        # arg is E or H for Eccentric or Hyperbolic case but B,p, r for parabolic
        true_anomaly_new = anomaly2nu(ecc, anomaly_type, arg)
        arg2 = "none"
    else:
        arg_latitude = anomaly
        arg2 = arg_latitude

    r_vec_IJK, v_vec_IJK = COE2RV(a, ecc, incl, ascending_node, arg_perigee, true_anomaly_new, mu, arg2)

    return r_vec_IJK, v_vec_IJK

################################################################################################
# this function will calculate the orbital elements of a orbiting body given a known position vector and its velocity
# at that position.
# Inputs must be a [1X3] and in km
# Outputs are in Km and degrees

# INPUTS
#   r0_vec      - [1x3] initial position
#   v0_vec      - [1x3] initial velocity
#   mu                  - Gravitational Constant, km^3/s^2

# OUTPUTS
#   p                   - Semi-Parameter, km
#   ecc                 - Eccentricity
#   incl                - Inclination, rad
#   ascending_node      - Ascending Node, rad
#   arg_perigee         - Argument of Perigee, rad
#   true_anomaly        - True Anomaly, rad
#   mu                  - Gravitational Constant, km^3/s^2
#   *args
#   lambda_true         - Lambda True (Circular Equatorial Orbit)
#   arg_latitude        - Argument of Latitude (Circular Inclined Orbit)
#   arg_perigee_true    - True Argument of Perigee (Equatorial Elliptical Orbit)

def RV2COE (r_ijk, v_ijk, mu):
    # Angular Momentum
    h_vec = np.cross(r_ijk.T, v_ijk.T).T # km^2/s
    h_mag = np.linalg.norm(h_vec) # km^2/s

    # Node Vector
    k_vec = np.array([[0, 0, 1]])
    n_vec = np.cross(k_vec, h_vec.T).T # km^2/s
    n_mag = np.linalg.norm(n_vec) # km^2/s

    # Eccentricity
    pos_mag = np.linalg.norm(r_ijk) # km
    vel_mag = np.linalg.norm(v_ijk) # km/s
    ecc_vec = (((vel_mag**2 - (mu / pos_mag)) * r_ijk) - (np.dot(r_ijk.T, v_ijk) * v_ijk)) / mu
    ecc = np.linalg.norm(ecc_vec)

    # Specific Orbital Energy
    energy = vel_mag**2 / 2 - (mu / pos_mag) # km^2/s^2

    # Semi-Major Axis and Semi-Minor Axis
    if ecc != 1:
        a = -(mu / (2 * energy)) # km
        p = a * (1 - ecc**2) # km

    else:
        p = h_mag**2 / mu # km
        a = float('inf') # km

    #Inclination
    incl = (np.acos(h_vec[2] / h_mag))

    #Ascending Node
    if n_mag > 0:
        ascending_node = (np.acos(n_vec[0] / n_mag))
    else:
        ascending_node = 0

    if n_vec[1] < 0:
        ascending_node = (2 * np.pi) - ascending_node

    # Argument of Perigee
    arg_perigee = (np.acos(np.dot(n_vec.T, ecc_vec) / (n_mag * ecc)))
    if ecc_vec[2] < 0:
        arg_perigee = (2 * np.pi) - arg_perigee


    # True Anomaly
    true_anomaly = (np.acos(np.dot(ecc_vec.T, r_ijk) / (ecc * pos_mag)))
    if np.dot(r_ijk.T, v_ijk) < 0:
        true_anomaly = (2 * np.pi) - true_anomaly

    #Special Cases
    # Elliptical and Equatorial - True Argument of Perigee
    arg_perigee_true = (np.acos(ecc_vec[0] / ecc))
    if ecc_vec[1] < 0:
        arg_perigee_true = (2 * np.pi) - arg_perigee_true

    # Circular Inclined - Argument of Latitude
    arg_latitude = (np.acos(np.dot(n_vec.T, r_ijk) / (n_mag * pos_mag)))
    if r_ijk[2] < 0:
        arg_latitude = (2 * np.pi) - arg_latitude

    # Circular Equatorial - True Lambda
    lambda_true  = (np.acos(r_ijk[0] / pos_mag))
    if r_ijk[1] < 0:
        lambda_true = (2 * np.pi) - lambda_true

    return a, ecc, incl, ascending_node, arg_perigee, true_anomaly, arg_perigee_true, arg_latitude, lambda_true

################################################################################################
# THIS FUNCTION WILL OUTPUT THE ECCENTRIC, PARABOLIC, AND HYPERBOLIC ANOMALY DEPENDING ON THE TYPE OF THE INPUT ORBIT.
# THIS FUNCTION NEEDS THE ECCENTRICITY AND THE TRUE ANOMALY TO GIVE THE ANOMALY NEEDED TO DERIVE THE FLIGHT PATH ANGLE.

# INPUTS
# e             - Eccentricity
# nu            - True Anomaly, rad

# OUTPUTS
# E             - Eccentric Anomaly
# P             - Parabolic Anomaly
# H             - Hyperbolic Anomaly

def nutoAnomaly (e, nu):
    if e < 1.0:
        #E1 = np.asin((np.sin(nu) * np.sqrt(1 - e**2)) / (1 + (e * np.cos(nu))))
        E = np.acos((e + np.cos(nu)) / (1 + (e * np.cos(nu))))
        return E
    elif e == 1:
        B = np.tan(nu / 2)
        return  B
    else:
        #H1 = np.asinh((np.sin(nu) * np.sqrt(e**2 - 1)) / (1 + (e * np.cos(nu))))
        H = np.acosh((e + np.cos(nu)) / (1 + (e * np.cos(nu))))
        return H
