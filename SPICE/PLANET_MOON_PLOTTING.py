import spiceypy as sp
import TOOLS as tls
import os

# LOADS THE LEAPSECOND KERNEL (LSK), BODY SIZE AND ORIENTATION KERNEL (PCK), AND THE EPHEMERIS KERNEL (SPK)
tls.load_kernels(r'planet_moon_kernels.tm')

#kerel file name
file_name = tls.load_spk("JUPITER", "jup365.bsp")

# EXTRACTING INFORMATION FROM SPK FILE
ids, names, tcs_sec, tcs_cal = tls.get_objects(file_name, display=True)

# DATES FOR ANALYSIS
dates = ['2025 JUN 1 00:12:00', '2025 JUN 20 00:00:00']
et0 = sp.str2et(dates[0])
etf = sp.str2et(dates[1])
steps = 10000

# MAKING TIME VECTOR
time_vec = tls.tvlist2array(et0, etf, steps)

body_numbers = [599, 501, 502, 503, 504, 505, 514, 515, 516]
jup_moons = ["IO", "EUROPA", "GANYMEDE", "CALLISTO", "AMALTHEA", "ADRASTEA", "METIS"]

# PRE ALLOCATING VARIABLE FOR PLANETARY DATA
body_data = []

# EXTRACTING THE POSITION DATA FOR EACH PLANET
for jup_moon in jup_moons:
    # Adding Ephem Data into the List
    body_data.append(tls.ephemdata(jup_moon, time_vec, 'ECLIPJ2000', '599'))

tls.orbitplot(body_data, jup_moons, AU=True)