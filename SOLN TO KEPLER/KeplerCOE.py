# This function provides a solution for Kepler's problem using classical orbital elements. This method does not include
# any perturbations. Therefore, the only value to change in the two positions is the true anomaly or the special orbit
# parameters.

# This method uses the known position and velocity to find the initial position's anomaly. Then determining the Mean anomaly.
# Finally using the mean anomaly to determine the true anomaly at the new position all information is present to
# determine the new position and velocity.

# INPUTS
#   r0_vec      - initial position
#   v0_vec      - initial velocity
#   delta_t     - time frame of observation

OUTPUTS
#   r_vec       - position vector at the observation time
#   v_vec       - velocity at new position


import numpy as np

