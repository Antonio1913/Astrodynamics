

import numpy as np

    def COE2RV (a, p, ecc, incl, ascending_node, arg_perigee, true_anomaly, arg_perigee_true, arg_latitude, lambda_true, mu):

        #Setting Conditional Terms
        # Circular and Equatorial
         if ecc == 0 and incl == 0:
            arg_perigee = 0
            ascending_node = 0
            true_anomaly = lambda_true

        # Circular and Inclined
        if ecc == 0 and incl < 0:
            arg_perigee = 0
            true_anomaly = arg_latitude

        # Circular and Equatorial
        if ecc < 0  and incl = 0
            ascending_node = 0
            arg_perigee = arg_perigee_true


    r_PQW = np.array([(p * np.cos(true_anomaly)) / (1 + (e * np.cos(true_anomaly))),
                      (p * np.sin(true_anomaly)) / (1 + (e * np.cos(true_anomaly)))
                         , 0])

    v_PQW = np.array([-np.sqrt(mu / p) * np.sin(true_anomaly),
                        np.sqrt(mu / p) * (e + np.cos(true_anomaly)),
                        0])



