# VALUES DISPLAYED ARE IN REFERENCE OF EPOCH J2000
class Earth:
    SMA = 149598023 #km
    ecc = 0.016708617
    Inc = 0 #degrees
    RAAN = 0 #degrees
    long_Perihelion = 102.93734808 #degrees
    True_Longitude = 100.46644851 #degrees
    Orbital_Period = 0.99997862 #years
    Orbital_Period_Tropical = 365.2421897#days
    Orbital_Velocity = 29.7859 #km/s
    Radius = 6378.1363 #km
    flattening = 0.0033528131
    eccentricity = 0.081819221456 # Shpe of Earth
    mu = 3.986004415 * 10**5 #km^3/s^2
    mass = 5.9742 * 10**24 #kg
    Rotation = 0.99726968 #days
    Incl_equatorial = 23.45 #degrees
    J2 = 0.0010826269
    angularvelocity = 7.2722 * 10**-5 # rad/s

class Semi_Major_Axis:
    Moon_SMA = 384400
    Mercury_SMA = 57909083
    Venus_SMA = 108208601
    Earth_SMA = 149598023
    Mars_SMA = 227939186
    Jupiter_SMA = 778298361
    Saturn_SMA = 1429394133
    Uranus_SMA = 2875038615
    Neptune_SMA = 4504449769
    Pluto_SMA = 5915799000

class Gravitational_Parameter: # km^3/s^2
    Moon_mu = 4902.799
    Mercury_mu = 2.2032 * 10e4
    Venus_mu = 3.257 * 105
    Earth_mu = 3.986004415 * 10e5
    Mars_mu = 4.305 * 10e4
    Jupiter_mu = 1.268 * 10e8
    Saturn_mu = 3.794 * 10e7
    Uranus_mu = 5.794 * 10e6
    Neptune_mu = 6.809 * 10e6
    Pluto_mu = 9.00 * 10e2