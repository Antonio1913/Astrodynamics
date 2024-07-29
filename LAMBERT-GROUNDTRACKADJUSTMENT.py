

# INPUTS
#   t1                          - Departing Time, sec
#   a_star                      - Semi-Major Axis Before Maneuver, km
#   ecc_star                    - Eccentricity Before Maneuver
#   incl_star                   - Inclination Before Maneuver, rad
#   ascending_node_star         - Right Ascension of the Ascending Node Before Maneuver, rad
#   arg_perigee_star            - Argument of Perigee Before Maneuver km
#   f_star                      - True Anomaly Before Maneuver, rad
#   L_G                         - Geocentric Longitude of Ground Target or Intersection Point, rad
#   phi_G                       - Geocentric Latitude of Ground Target or Intersection Point, rad
#   r_min                       - Minimum Position Vector, km
#   r_max                       - Maximum Position Vector, km
#   D                           - Parameter (integer) of Time of Flight
#   incl                        -
#   H_2                         - Altitude at Time 2, km
#   AA                          -
#   kappa                       - Transfer Direction
#   K                           - Mode Parameter

# OUTPUTS


import numpy as np
from BODY_VALUES.MEAN_PLANETARY_CONSTANTS import Earth as E
from SOLN_TO_KEPLER.COE2RV import COE2RV


def ExtendedLambert(t1, a_star, ecc_star, incl_star, ascending_node_star, arg_perigee_star, true_anomaly_star, L_G, phi_G,
                     r_min, r_max, D, incl, H_2, kappa, K, *arg):

#   Defining r1 using the initial Satellite Parameters
    r1_vec, v1_vec = COE2RV (a_star, ecc_star, incl_star, ascending_node_star, arg_perigee_star, true_anomaly_star, E.Earth_mu, "none")
    r1 = np.linalg.norm(r1_vec)

#   Argument of Latitude of the Spacecrafts Initial Orbit at T1
    u_star = arg_perigee_star + true_anomaly_star

#   Phi_1 = Phi_star - Geocentric Latitude of the Departing Point P1
    phi_star = np.asin(np.sin(incl_star) * np.sin(u_star))

#   Right Ascension of the Ascending Node at Time 1 at the time of Impulsive Maneuver
    ascending_node_1 = (ascending_node_star + (np.asin(np.tan(phi_star) / np.tan(incl_star))) -
                        (kappa * (np.asin(np.tan(phi_star) / np.tan(incl)))) + ((1 - kappa) * (np.pi / 2)))

#   Unit Vector in the Direction of the Ascending Node at Time 1
    i_vec_Omega1 = np.array([np.cos(ascending_node_1), np.sin(ascending_node_1), 0])

#   Unit Vector in the Direction of Position Vector r1 or rstar
    i_vec_r1 = np.array([(np.cos(ascending_node_star) * np.cos(u_star)) - (np.sin(ascending_node_star) * np.sin(u_star) * np.cos(incl_star)),
                         (np.sin(ascending_node_star) * np.cos(u_star)) + (np.cos(ascending_node_star) * np.sin(u_star) * np.cos(incl_star)),
                         np.sin(u_star) * np.sin(incl_star)])

#   Argument of Latitude at Time 1
    u_1 = np.acos(np.dot(i_vec_Omega1, i_vec_r1))

#   Calculating the Central Angle between Departing Point P1 and Arriving Point P2 and thre longitudinal difference
#   between target G and the intersection G', delta_l
#   Theta is Dependent on whether the input was for the ascending stage (AA) or descending stage (DA)

    # Converting time in UTC into Greenwich Mean Sidereal Time(GMST)

    Theta_gmst_t1 =

    if arg == "AA":
        theta = np.mod(np.asin(np.sin(phi_G / np.sin(incl)) - u_1), 2 * np.pi)
        delta_l = np.mod(ascending_node_1 + (np.asin(np.tan(phi_G) / np.tan(incl))) - L_G - Theta_gmst_t1, 2 * np.pi) + 2 * np.pi * D

    elif arg == "DA":
        theta = np.mod(np.pi - (np.asin(np.sin(phi_G / np.sin(incl)) - u_1)), 2 * np.pi)
        delta_l = np.mod(ascending_node_1 + (np.pi - (np.asin(np.tan(phi_G) / np.tan(incl)))) - L_G - Theta_gmst_t1, 2 * np.pi) + 2 * np.pi * D

    else:
        print(f"Incorrect Orbit Adjustment. Only valid inputs, 'DA', 'AA'.")

#   Defining the two position vectors
    r2 = E.Earth_Radius + H_2
    lambda_max = np.sqrt((r1 + r2 - abs(r1 - r2)) / (r1 + r2 + abs(r1 - r2)))

#Determining Keplerian Starters, lambda, x0
# Step 1 - Setting a0 = a_star, ecc0 = ecc_star
    a0 = a_star
    ecc0 = ecc_star

# Step 2 - Solving Equations 5 and 9 to compute TOF
#   Calculating the secular drift rate of the right ascension of the ascending node caused by J2 Perturbation
    B_Omega = 1.5 * E.Earth_J2 * E.Earth_Radius**2 * np.sqrt(E.Earth_mu) * np.cos(incl)
    B_e = 1 - ecc0**2

    Omega_J2 = B_Omega * a0**(-7 / 2) * B_e**-2

#   Secular drift rate of the argument of perigee caused by J2 Perturbation
    B_omega = -1.5 * E.Earth_J2 * E.Earth_Radius**2 * np.sqrt(E.Earth_mu) * ((2.5 * np.sin(incl**2)) - 2)
    omeegadot_J2 = B_omega * a0**(-7 / 2) * B_e**-2

#   Initial Guess for TOF
    tau_0 = delta_l / (E.Earth_Rotation - Omega_J2)

#   Rotation Angle Between t1 and t2 of the Flyby Orbit
    theta_omega = omeegadot_J2 * tau_0

#   Auxiliary Angle
    theta_C_hat = theta - theta_omega

#Ensuring input is within range
    low_range_k = -theta_C_hat / (2 * np.pi)
    upper_range_k = 1 - (theta_C_hat / (2 * np.pi))
    if K < low_range_k or K > upper_range_k:
        print(f"K input is not an ideal value with other input parameters.")

#   Initial Guess for Transfer angle
    theta_C0 = np.mod(theta_C_hat, 2 * np.pi)

#   Replacing Theta_C with theta_C0 in Equations 15 21 26 to compute initial guesses for Semi perimeter and the parameter lambda
    c = np.sqrt((r1 + r1)**2 - 4 * r1 *r2 * np.cos((theta_C0 / 2)**2))
    s0 = (r1 + r2 + c) / 2
    lambda0 = (np.sqrt(r1 * r2) * np.cos(theta_C0 / 2)) / s0

# Determining the Allowable Number of Revolutions using modified Newton's Iteration
#   Defining iteration methods
    max_iterations = 100
    tolerance = 1 * 10**-6

    for i in range(max_iterations):

        # Defining initial guess for x0
        x0 = 0

        # Calculating other values needed for the partial derivatives below
        sigma = 1 - x0**2
        y = np.sqrt(1 - lambda0**2 * (1 - x0**2))
        s = (r1 + r2) / (1 + lambda0**2)
        B_M = -1.5 * E.Earth_J2 * E.Earth_Radius**2 * (1.5 * np.sin(incl**2) - 1)
        zeta = (1 + B_M * a0**-2 * B_e**(-3/2) * 0.5 * np.sqrt(E.Earth_mu) * tau_0)
        cos_psi = (x0 * y) + (lambda0 * (1 - x0**2))
        theta_ck = theta - (B_Omega * a0**(-7/2) * B_e**-2 * tau_0) + (2 * np.pi * K)

        # Calculating Auxiliary Function, W
        W = lambda0 - (np.sqrt(r1 * r2))**-1 * np.cos(((theta_C0 - (omeegadot_J2 * tau_0)) / 2) + (np.pi * K))

        # Calculating partial derivatives
        dLambda_F_dlambda = (sigma**-1 * y) - (lambda0**2 * y**-1) - (6 * np.sqrt(2) * s**(-3/2) * lambda0 * (1 +  lambda0**2)**-1 * zeta)
        dLambda_upsilon_dlambda = -np.pi**-1 * sigma**(3/2) * dLambda_F_dlambda

        dLambda_F_da = (2 * np.sqrt(2 * E.Earth_mu) * B_M * a0**-3 * B_e**(-3/2) * s(-3/2) * tau_0) + (7 * np.sqrt(2) * B_Omega * delta_l**-1 * a0**(-9/2) * B_e**-2 * s(-3/2) * tau_0 * zeta)
        dLambda_upsilon_da = -np.pi**-1 * sigma**(3/2) * dLambda_F_da

        da_dlambda = -2 * a0 * lambda0 * (1 + lambda0**2)**-1

        dLambda_F_dB_e = (1.5 * np.sqrt(2 * E.Earth_mu) * B_M * a0**(-2) * B_e**(5/2) * s(-3/2) * tau_0) + (4 * np.sqrt(2) * B_Omega * delta_l**-1 * a0**(-7/2) * B_e**-3 * s**(-3/2) * tau_0 * zeta)
        dLambda_Upsilon_dB_e = -np.pi**-1 * sigma**(3/2) * dLambda_F_dB_e

        dB_e_dlambda = (-8 * sigma * ((s**2 * lambda0**5 * sigma) + (lambda0**3 * ((s**2 * y**2) -
                    (sigma * ((r1 * r2) - s**2)))) + (lambda0 * (y**2 * (s**2 - (2 * r1 * r2)) - ((r1 * r2 * sigma)))) -
                    (r1 * r2 * x0 * y * (1 - lambda0**2)))) / (s**2 * (1 + lambda0**2) * (y - (lambda0 * x0))**3 * y)

        dLambda_Upsilon_dpsi = -np.pi**-1

        dpsi_dlambda = (((lambda0 * x0**2) - (x0 * y) - lambda0) * sigma**(1/2) * y**-1) / cos_psi

        dLambda_Upsilon_dx = np.pi**-1 * ((-6 * np.sqrt(2) * s**(-3/2) * x0 * sigma**(1/2) * zeta) + ((1 - (lambda0**3 * x0 * y**-1)) * sigma**(1/2)) - ((x0 - (lambda0 * y)) * x0 * sigma**(-1/2)))

        da_dx = 2 * a0 * x0 * sigma**-1

        dB_e_d_x = -2* B_e * ((-lambda0 * y**-1) + (x0 * sigma**-1))

        dPsi_dx = (((x0 - (2 * x0 * y**2)) + (((2 * lambda0 * x0**2) - lambda0) * y)) * sigma**(-1/2) * y**-1) / cos_psi

        dLambda_W_dx = 0

        dLambda_W_da = (7/4) * np.sqrt(r1 * r2) * B_omega * a0**(-9/2) * B_e**(-2) * s**(-1) * tau_0(1 + (B_Omega * a0**(-7/2) * B_e**(-2) * delta_l**(-1) * tau_0)) * np.sin(theta_ck / 2)

        dLambda_W_d_B_e = np.sqrt(r1 * r2) * B_omega * a0**(-7/2) * B_e**(-3) * s**(-1) * tau_0(1 + (B_Omega * a0**(-7/2) * B_e**(-2) * delta_l**(-1) * tau_0)) * np.sin(theta_ck / 2)

        dLambda_W_dlambda = (1 - lambda0**2 + (2 * lambda0 * W)) * (1 + lambda0**2)**-1


    # Values needed to calculate dUpsilon_dx
        dPsi_Upsilon_dlambda = dLambda_upsilon_dlambda + (dLambda_upsilon_da * da_dlambda) + (dLambda_Upsilon_dB_e * dB_e_dlambda) + (dLambda_Upsilon_dpsi * dpsi_dlambda)
        dPsi_Upsilon_dx = dLambda_Upsilon_dx + (dLambda_upsilon_da * da_dx) + (dLambda_Upsilon_dB_e * dB_e_d_x) + (dLambda_Upsilon_dpsi * dPsi_dx)

    #   Values needed to calculate dlambda_dx
        dW_dx = dLambda_W_dx + (dLambda_W_da * da_dx) + (dLambda_W_d_B_e * dB_e_d_x)
        dW_dlambda = dLambda_W_dlambda + (dLambda_W_da * da_dlambda) + (dLambda_W_d_B_e * dB_e_dlambda)
        dlambda_dx = -dW_dx / dW_dlambda

        dUpsilon_dx = (dPsi_Upsilon_dlambda * dlambda_dx) + dPsi_Upsilon_dx

        jacobian = np.array([[dPsi_Upsilon_dlambda, dPsi_Upsilon_dx], [dW_dlambda, dW_dx]])
        jacobian_inv = np.linalg.inv(jacobian)
        jacobian_inv_ident = np.identity(jacobian_inv)

    #   Damping Factor
        rho_Upislonx = 0.5

     #  Equation 60
        lambdaplus1 = lambda0 - (rho_Upislonx * jacobian_inv_ident * dUpsilon_dx)
        xplus1 = x0 - (rho_Upislonx * jacobian_inv_ident * W)


        if abs(lambda0 - lambdaplus1) < tolerance and abs(x0 - xplus1) < tolerance:
            lambda_star = lambdaplus1
            x_star = xplus1
            break
        lambda0 = lambdaplus1
        x0 = xplus1


    Nmaxtype1 =
















