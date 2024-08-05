
import spiceypy as spice
import numpy as np
from SPICE.SPICE_TOOLS import load_kernels, get_objects, tvlist2array, ephemdata
from TOOLS.PLOTTING_TOOLS import orbitplot
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


load_kernels('solar_system_kernels.tm')
ids, names, tcs_sec, tcs_cal = get_objects(r'..//SPICE_KERNELS\de440.bsp', display=True)

dates = ['1849 DEC 26 00:12:00', '2125 Dec 31 00:00:00']

et0 = spice.str2et(dates[0])
etf = spice.str2et(dates[1])
steps = 10000

time_vec = tvlist2array(et0, etf, steps)

names = [f for f in names if 'BARYCENTER' in f]
body_data = []

for name in names:

    # Adding Ephem Data into the List
    body_data.append(ephemdata(name, time_vec, 'ECLIPJ2000', '10'))

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
#
# for bodies in body_data:
#     ax.plot(bodies[:, 0], bodies[:, 1], bodies[:, 2])
#     ax.plot([bodies[0, 0]], [bodies[0, 1]], [bodies[0, 2]], 'o')
#
# plt.show()


orbitplot(body_data, names, AU=True)



