import numpy as np
import spiceypy as spice
import os


# Finds the current metakernel
def load_kernels(file_name):
    # Define the path to the metakernel file
    metakernel_path = os.path.abspath(file_name)

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


# Gets Ephemeris Data Using Spice.spkezr
def ephemdata(target, time, ref, observer):
    return np.array(spice.spkezr(target, time, ref, 'NONE', observer)[0])

# LOADS THE PATH OF THE SPK FILE TO EXTRACT INFORMATION WITHIN USING "get_objects" FUNCTION
# MAKE SURE TO USE .bsp AT THE END OF THE SPK_NAME
def load_spk(kernel_base_folder_name, spk_name):

    # DEFINE BASE PATH OF COMPUTER
    base_path = os.path.expanduser("~")

    # USE INSTRUCTIONS LOCATED IN THE SPICE FOLDER IN THE SPK_NOTES.TXT FILE
    # LOADING THE PATH OF THE SPK_NAME USING PREDETERMINED PATHS
    kernel_path = f"{base_path}/ASTRO/KERNELS/{kernel_base_folder_name}/{spk_name}"

    return kernel_path


