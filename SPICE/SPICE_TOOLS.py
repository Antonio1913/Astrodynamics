import numpy as np
import spiceypy as spice
import os


# Finds the current metakernel
def load_kernels(metakernel_filename='solar_system_kernels.tm'):
    # Define the path to the metakernel file
    script_dir = os.path.dirname(__file__)
    kernel_dir = os.path.join(script_dir, '..\\SPICE_KERNELS')
    # kernel_dir = r'C:\Users\antonio.garcia\PycharmProjects\Astrodynamics\SPICE_KERNELS'
    metakernel_path = os.path.join(kernel_dir, metakernel_filename)

    # Ensure the metakernel file exists
    if not os.path.exists(metakernel_path):
        raise FileNotFoundError(f'Metakernel file not found: {metakernel_path}')

    # Load the metakernel
    spice.furnsh(metakernel_path)
    print(f'Loaded metakernel from {metakernel_path}')


# test to ensure load_kernels is operating correctly
def main():
    load_kernels()

    # Example usage of a SPICE function
    utc_time = '2004 jun 11 19:32:00'
    et = spice.str2et(utc_time)
    print(f'ET Seconds Past J2000: {et:.3f}')

    # Unload kernels after use
    spice.kclear()
    return et



# Pulls ID, names and time coverages of the objects in the spk file
def get_objects(filename, display=False):
    objects = spice.spkobj(filename)
    ids, names, tcs_sec, tcs_cal = [], [], [], []
    n = 0
    if display:
        print('\nObjects in %s: ' % filename)

    for o in objects:
        # Retrieving ID number
        ids.append(o)

        # Time Coverage Available for Objects since J2000
        tc_sec = spice.wnfetd(spice.spkcov(filename, ids[n]), n)

        # Converts Time Coverage into Readable Form
        tc_cal = [spice.timout(f, "YYYY MON DD HR:MM:SC.### (TDB) ::TDB") for f in tc_sec]

#       Append Time Coverage to Output List
        tcs_sec.append(tc_sec)
        tcs_cal.append(tc_cal)

#       Getting the Name of Bodies
        try:
            # Adding the Names to the Output List
            names.append(spice.bodc2n(o))

        except:
            # Called if Body Name Does Not Exist
            names.append('Unknown Name')

        if display:
            print('id: %i\t\tname: %s\t\t\ttc: %s --> %s' % (ids[-1], names[-1], tc_cal[0], tc_cal[1]))

    return ids, names, tcs_sec, tcs_cal


# Converts Range of Time Into An Array
def tvlist2array(et0, etf, steps):
    time_vec = np.zeros((steps, 1))
    time_vec[:, 0] = np.linspace(et0, etf, steps)
    return time_vec

# Gets Ephemeris Data Using Spice.spekzr
def ephemdata(target, time, ref, observer):
    return np.array(spice.spkezr(target, time, ref, 'NONE', observer)[0])



