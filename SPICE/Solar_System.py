import spiceypy as spice
import TOOLS as tls

# LOADING KERNEL FILES
tls.load_kernels('solar_system_kernels.tm')

# EXTRACTING INFORMATION FROM SPK FILE
ids, names, tcs_sec, tcs_cal = tls.get_objects(r'..\\SPICE\\SPICE_KERNELS\de440.bsp', display=True)

# DATES FOR ANALYSIS
dates = ['1950 DEC 26 00:12:00', '2025 Dec 31 00:00:00']
et0 = spice.str2et(dates[0])
etf = spice.str2et(dates[1])
steps = 10000

# MAKING TIME VECTOR
time_vec = tls.tvlist2array(et0, etf, steps)

# GRABBING THE NAMES FOR THE PLANETS THAT ONLY CONTAIN BARYCENTER
names = [f for f in names if 'BARYCENTER' in f]

# PRE ALLOCATING VARIABLE FOR PLANETARY DATA
body_data = []

# EXTRACTING THE POSITION DATA FOR EACH PLANET
for name in names:
    # Adding Ephem Data into the List
    body_data.append(tls.ephemdata(name, time_vec, 'ECLIPJ2000', '10'))

tls.orbitplot(body_data, names, AU=True)



