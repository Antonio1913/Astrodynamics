import numpy as np

# check_impact IS A FUNCTION THAT WILL DETERMINE A FLYBY SURFACE GRACE OR A DIRECT IMPACT
#
# INPUTS:
#
#   target planet           - NAME OF TARGET PLANET
#   V_esc                   - ESCAPE VELOCITY OF TARGET PLANET
#   V_inf                   - THE RESIDUAL VELOCITY OF VEHICLE APPROACHING TARGET PLANET
#   phi                     - ELEVATION ANGLE OF VEHCILE APPROACHING TARGET PLANET

def check_impact(target_planet, V_esc, V_inf, phi):

