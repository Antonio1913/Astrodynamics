

import numpy as np

    def COE2RV (a, p, ecc, incl, ascending_node, arg_perigee, true_anomaly, arg_perigee_true, arg_latitude, lambda_true):

        #Setting Conditional Terms
        # Circular and Equatorial
         if ecc == 0 and incl == 0
             set(arg_perigee, ascending_node) == 0