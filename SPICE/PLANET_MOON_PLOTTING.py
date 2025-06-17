import spiceypy as sp
import TOOLS as tls

tls.load_kernels(r'planet_moon_kernels.tm')

# EXTRACTING INFORMATION FROM SPK FILE
ids, names, tcs_sec, tcs_cal = tls.get_objects(r'C:\Users\antoniogarcia778e\ASTRO\KERNELS\JUPITER\jup365.bsp', display=True)

# DATES FOR ANALYSIS
dates = ['2025 JUN 1 00:12:00', '2025 JUN 20 00:00:00']
et0 = sp.str2et(dates[0])
etf = sp.str2et(dates[1])
steps = 10000

# MAKING TIME VECTOR
time_vec = tls.tvlist2array(et0, etf, steps)

body_numbers = [599, 501, 502, 503, 504, 505, 514, 515, 516]
names = ["IO", "EUROPA", "GANYMEDE", "CALLISTO", "AMALTHEA", "ADRASTEA", "METIS"]

# PRE ALLOCATING VARIABLE FOR PLANETARY DATA
body_data = []

# EXTRACTING THE POSITION DATA FOR EACH PLANET
for name in names:
    # Adding Ephem Data into the List
    body_data.append(tls.ephemdata(name, time_vec, 'ECLIPJ2000', '599'))

tls.orbitplot(body_data, names, AU=True)