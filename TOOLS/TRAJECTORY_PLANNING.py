import numpy as np
import TOOLS as tls

# check_impact IS A FUNCTION THAT WILL DETERMINE A FLYBY SURFACE GRACE OR A DIRECT IMPACT
# USE REF 1 TO FIND MORE INFORMATION
#
# INPUTS:
#
#   target planet           - NAME OF TARGET PLANET
#   V_inf                   - THE RESIDUAL VELOCITY OF VEHICLE APPROACHING TARGET PLANET
#   phi                     - ELEVATION ANGLE OF VEHICLE APPROACHING TARGET PLANET

# noinspection PyInconsistentReturns
def check_impact(target_planet, V_inf, phi):

    # ENSURES THAT THE INPUT IS CORRECTED TO THE EXPECTED FORMAT
    class_name = target_planet.capitalize()

    # RETRIEVES THE CLASS INFORMATION FOR THE TARGET PLANET
    planet_class = getattr(tls.BODY_CONSTANTS, class_name)

    # DEFINING THE TWO VARIABLES NEEDED FROM THE TARGET PLANET CLASS
    r_SOI = planet_class.SOI
    radius = planet_class.Radius
    V_esc = planet_class.V_esc

    # CALCULATING THE APPROACH DISTANCE - EQ. 5-3-10
    d = r_SOI * np.cos(phi)

    # CALCULATING THE IMPACT PARAMETER - EQ. 5-3-9
    b = radius * np.sqrt(1 + (V_esc**2 / V_inf**2))

    # COMPARISON OF D AND B WILL DETERMINE IF THERE IS IMPACT
    if d > b:
        # RETURNING 0 MEANS NO IMPACT
        impact_param = 0
        return impact_param, "THERE WILL BE NO IMPACT."

    elif d == b:
        # RETURNING 1 MEANS THERE WILL BE A SURFACE GRAZE
        impact_param = 1
        return impact_param, "THERE WILL BE A SURFACE GRAZE."

    else:
        # RETURNING 2 MEANS THERE WILL BE AN IMPACT
        impact_param = 2
        return impact_param, "THERE WILL BE IMPACT."



