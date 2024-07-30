

from __future__ import  print_function
from builtins import input
import spiceypy as spice

def convtm():
    METAKR = r'C:\path\to\your\convtm.tm'
    SCLKID = -82

    spice.furnsh( METAKR )

    utcim = input( ' Input UTC Time:' )

    print( 'Converting UTC Time: {:s}'.format( utcim))

    et = spice.str2et( utcim )

    print( ' ET Seconds Past J2000: {:16.3f}'.format(et))

    calet = spice.etcal(et)

    print( ' Calender ET (etcal): {:s}'.format( calet ))

    calet = spice.timout(et, 'YYYY-MON-DDTHR:MN:SC :: TDB00')

    print(' Calender ET (timeout): {:s}'.format(calet))

    sclkst = spice.sce2s( SCLKID, et)

    print( ' Spacecraft Clock Time: {:s}'.format( sclkst))

    spice.unload( METAKR )

if __name__ == '__main__':
    convtm()





