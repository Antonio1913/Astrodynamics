# VALUES DISPLAYED ARE IN REFERENCE OF EPOCH J2000 AND COE FROM SEIDELANN (1992:704) AND (1992:706)

class Mercury:
    SMA = 57909083  # km
    ecc = 0.205631752
    Inc = 7.00498625  # degrees
    RAAN = 48.33089304  # degrees
    long_Perihelion = 77.45611904  # degrees
    True_Longitude = 252.25090551  # degrees
    Orbital_Period = 0.24084445  # years
    Orbital_Period_Tropical = 87.9666  # days
    Orbital_Velocity = 47.8725  # km/s
    Radius = 2439.0  # km
    mu = 2.2032 * 10**4  # km^3/s^2
    mass = 3.3022 * 10**23  # kg
    Rotation = 58.6462  # days
    Incl_equatorial = 0.0  # degrees
    J2 = 0.00006
    Density = 5.43 # gm/cm^3)
    SOI = 1.117 * 10 **5 # km WITH RESPECT TO THE SUN

class Venus:
    SOI = 6.163 * 10 **5 # km WITH RESPECT TO THE SUN

class Moon:
    SMA = 384,400  # km
    ecc = 0.05490
    Inc = 5.145396  # degrees
    Orbital_Period = 0.0748 # years
    Orbital_Period_Tropical = 27.321582  # days
    Orbital_Velocity = 1.0232  # km/s
    Radius = 1738.0  # km
    mu = 4902.799  # km^3/s^2
    mass = 7.3483 * 10 ** 22  # kg
    Rotation = 27.32166 # days
    Incl_equatorial = 6.68  # degrees
    J2 = 0.0002027
    density = 3.34 # gm/cm^3)

class Earth:
    SMA = 149598023  # km
    ecc = 0.016708617
    Inc = 0  # degrees
    RAAN = 0  # degrees
    long_Perihelion = 102.93734808  # degrees
    True_Longitude = 100.46644851  # degrees
    Orbital_Period = 0.99997862  # years
    Orbital_Period_Tropical = 365.2421897  # days
    Orbital_Velocity = 29.7859  # km/s
    Radius = 6378.1363  # km
    flattening = 0.0033528131
    eccentricity = 0.081819221456  # Shape of Earth
    mu = 3.986004415 * 10**5  # km^3/s^2
    mass = 5.9742 * 10**24  # kg
    Rotation = 0.99726968  # days
    Incl_equatorial = 23.45  # degrees
    J2 = 0.0010826269
    angularvelocity = 7.2722 * 10**-5  # rad/s
    SOI = 9.245 * 10 ** 5  # km WITH RESPECT TO THE SUN

    class Mars:
        SOI = 5.781 * 10 ** 5  # km WITH RESPECT TO THE SUN

    class Jupiter:
        SOI = 4.820 * 10 ** 7  # km WITH RESPECT TO THE SUN

    class Saturn:
        SOI = 5.627 * 10 ** 7  # km WITH RESPECT TO THE SUN

    class Uranus:
        SOI = 5.172 * 10 ** 7  # km WITH RESPECT TO THE SUN

    class Neptune:
        SOI = 8.664 * 10 ** 7  # km WITH RESPECT TO THE SUN

    class Pluto:
        SOI = 3.320 * 10 ** 7  # km WITH RESPECT TO THE SUN