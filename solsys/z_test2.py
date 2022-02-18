from core import day, sun
import numpy as np
from datetime import datetime

d2r = np.pi/180
r2d = 180/np.pi

d = day(1990, 4, 19, 0)
#d = day(2022, 2, 11, 10)
t = datetime(1990, 4, 19, 0)


def radec_to_altaz(ra, dec, obs_loc, t):
    lon, lat = obs_loc

    J2000 = datetime(2000,1,1,12)
    d = (t - J2000).total_seconds() / 86400 #day offset

    UT = t.hour + t.minute/60 + t.second/3600
    LST = (100.46 + 0.985647 * d + lon + 15*UT + 360) % 360
    ha = (LST - ra + 360) % 360
    
    x = np.cos(ha*d2r) * np.cos(dec*d2r)
    y = np.sin(ha*d2r) * np.cos(dec*d2r)
    z = np.sin(dec*d2r)
    xhor = x*np.cos((90-lat)*d2r) - z*np.sin((90-lat)*d2r)
    yhor = y
    zhor = x*np.sin((90-lat)*d2r) + z*np.cos((90-lat)*d2r)
    az = np.arctan2(yhor, xhor)*r2d + 180
    alt = np.arcsin(zhor)*r2d
    return az, alt

s = sun(d)
obs_loc = (2.2945, 48.8584)
print('ra, dec:', (s.ra, s.dec))
print('az, alt:', s._get_altaz(obs_loc))


az, alt = radec_to_altaz(s.ra, s.dec, obs_loc, t)
print('az, alt:', (az, alt))

from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
obs = EarthLocation(lon=obs_loc[0], lat=obs_loc[1])
tt = Time(t)
c = SkyCoord(ra=s.ra, dec=s.dec, unit='deg')
cc = c.transform_to(AltAz(obstime=tt, location=obs))
print('az, alt:', (cc.az.value, cc.alt.value))


