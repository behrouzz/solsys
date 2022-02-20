import numpy as np
from numpy import pi, sin, cos, tan, sqrt, arctan2, arcsin, arctan, arccos, log10
from .orbital_elements import elements, jupiter_oe, saturn_oe, uranus_oe
from .utils import *
from .transform import cartesian_to_spherical, spherical_to_cartesian, ecliptic_to_equatorial, elements_to_ecliptic, radec_to_altaz
from .correction import moon_perts, topocentric, jupiter_lon_perts, saturn_lon_perts, saturn_lat_perts, uranus_lon_perts

class sun:
    """
    Sun positional parameters
    
    Parameters
    ----------
        d (datetime): time of observation
        epoch (int): year of epoch 
        obs_loc : tuple of observer location (longtitude, latitude)

    Attributes
    ----------
        ecl_sph  : ecliptic spherical coordinates (lon, lat, r)
        ecl_car  : ecliptic cartesian coordinates (x, y, z)
        equ_car  : equatorial cartesian coordinates (x, y, z)
        equ_sph  : equatorial spherical coordinates (ra, dec, r)
        alt      : azimuth
        az       : altitude
    """
    def __init__(self, t, obs_loc=None, epoch=None):
        self.name = 'sun'
        d = datetime_to_day(t)
        ecl = obl_ecl(d)
        self.d = d
        N,i,w,a,e,M = elements(self.name, d)
        self.elements = {'N':N, 'i':i, 'w':w, 'a':a, 'e':e, 'M':M}
        self.L = rev(w+M)
        self.ecl_car = elements_to_ecliptic('sun', N,i,w,a,e,M)
        self.ecl_sph = cartesian_to_spherical(self.ecl_car)
        if epoch is not None:
            self.ecl_sph[0] = self.ecl_sph[0] + 3.82394E-5 * (365.2422 * (epoch-2000) - d)
            self.ecl_car = spherical_to_cartesian(self.ecl_sph)
        self.equ_car = ecliptic_to_equatorial(self.ecl_car, d)
        self.equ_sph = cartesian_to_spherical(self.equ_car)
        self.ra, self.dec, self.r = self.equ_sph # just for ease of use

        if obs_loc is None:
            self.az, self.alt = None, None
        else:
            self.az, self.alt = radec_to_altaz(self.ra, self.dec, obs_loc, t)
        

class moon:
    """
    Moon positional parameters

    Parameters
    ----------
        d (datetime): time of observation
        epoch (int): year of epoch 
        obs_loc : tuple of observer location (longtitude, latitude)

    Attributes
    ----------
        elements : dictionary of orbital elements

        geo_ecl_car  : Geocentric ecliptic cartesian coordinates (x, y, z)
        geo_ecl_sph  : Geocentric ecliptic spherical coordinates (lon, lat, r)
        geo_equ_car  : Geocentric equatorial cartesian coordinates (x, y, z)
        geo_equ_sph  : Geocentric equatorial spherical coordinates (ra, dec, r)
        ra           : Right Ascension (GCRS)
        dec          : Declination (GCRS)
        r            : distance to earth (in Earth Radii)
        alt          : azimuth
        az           : altitude
        elongation   : elongation
        FV           : phase angle (0:full, 90:half, 180:new)
    """
        
    def __init__(self, t, obs_loc=None, epoch=None):
        self.name = 'moon'
        d = datetime_to_day(t)
        ecl = obl_ecl(d)
        #self.obs_loc = obs_loc
        self._sun = sun(t=t, obs_loc=obs_loc, epoch=epoch)
        N,i,w,a,e,M = elements(self.name, d)
        self.elements = {'N':N, 'i':i, 'w':w, 'a':a, 'e':e, 'M':M}
        geo_ecl_car = elements_to_ecliptic(self.name, N,i,w,a,e,M) # OK
        self.geo_ecl_sph = cartesian_to_spherical(geo_ecl_car) # OK
        
        self.L = rev(N+w+M) # Moon's mean longitude

        # Correcting perturbations
        self.geo_ecl_sph = moon_perts(self.geo_ecl_sph, self._sun.elements, self.elements)
        
        self.geo_ecl_car = spherical_to_cartesian(self.geo_ecl_sph)
        self.geo_equ_car = ecliptic_to_equatorial(self.geo_ecl_car, d)
        self.geo_equ_sph = cartesian_to_spherical(self.geo_equ_car)
        
        
        if obs_loc is None:
            self.az, self.alt = None, None
        else:
            self.geo_equ_sph = topocentric(obs_loc, self._sun.L, self.geo_equ_sph, d)
            self.az, self.alt = radec_to_altaz(self.geo_equ_sph[0], self.geo_equ_sph[1], obs_loc, t)

        self.ra, self.dec, self.r = self.geo_equ_sph # For ease of use
        
        self.elongation = arccos( cos((self._sun.ecl_sph[0]-self.geo_ecl_sph[0])*rd) * cos(self.geo_ecl_sph[1]*rd) )*(180/pi)
        self.FV = 180 - self.elongation
        

class planet:
    """
    Planets positional parameters
    
    Parameters
    ----------
        name (str)      : name of the planet
        t (datetime)    : time of observation
        obs_loc (tuple) : observer location (longtitude, latitude)
        epoch (int)     : year of epoch

    Attributes
    ----------
        hel_ecl_sph  : Heliocentric ecliptic spherical coordinates (lon, lat, r)
        hel_ecl_car  : Heliocentric ecliptic cartesian coordinates (x, y, z)
        geo_ecl_car  : Geocentric ecliptic cartesian coordinates (x, y, z)
        geo_ecl_sph  : Geocentric ecliptic spherical coordinates (lon, lat, r)
        geo_equ_car  : Geocentric equatorial cartesian coordinates (x, y, z)
        geo_equ_sph  : Geocentric equatorial spherical coordinates (ra, dec, r)
        ra           : Right Ascension (GCRS)
        dec          : Declination (GCRS)
        r            : distance to earth (in AU)
        alt          : azimuth
        az           : altitude
        elongation : elongation
        FV         : phase angle
        mag        : Apparent magnitude
        diameter   : Apparent diameter
    """
    def __init__(self, name, t, obs_loc=None, epoch=None):
        self.name = name.lower()
        #self.obs_loc = obs_loc
        d = datetime_to_day(t)
        ecl = obl_ecl(d)
        self._sun = sun(t=t, obs_loc=obs_loc, epoch=epoch)
        N,i,w,a,e,M = elements(self.name, d)
        self.elements = {'N':N, 'i':i, 'w':w, 'a':a, 'e':e, 'M':M}
        hel_ecl_car = elements_to_ecliptic(self.name, N,i,w,a,e,M)
        lon, lat, r = cartesian_to_spherical(hel_ecl_car)

        # Correcting perturbations of Jupiter, Saturn and Uranus
        if self.name in ['jupiter', 'saturn', 'uranus']:
            Mj = jupiter_oe(d)[-1]
            Ms = saturn_oe(d)[-1]
            Mu = uranus_oe(d)[-1]
            if self.name=='jupiter':
                lon = lon + jupiter_lon_perts(Mj, Ms, Mu)
            elif self.name=='saturn':
                lon = lon + saturn_lon_perts(Mj, Ms, Mu)
                lat = lat + saturn_lat_perts(Mj, Ms, Mu)
            elif self.name=='uranus':
                lon = lon + uranus_lon_perts(Mj, Ms, Mu)
        
        # Precession
        if epoch is not None:
            lon = lon + 3.82394E-5 * (365.2422 * (epoch-2000) - d)

        # heliocentric
        self.hel_ecl_sph = np.array([lon, lat, r])
        self.hel_ecl_car = spherical_to_cartesian(self.hel_ecl_sph)

        # To geocentric
        self.geo_ecl_car = self._sun.ecl_car + self.hel_ecl_car # sun check shavad
        self.geo_ecl_sph = cartesian_to_spherical(self.geo_ecl_car)
        self.geo_equ_car = ecliptic_to_equatorial(self.geo_ecl_car, d)
        self.geo_equ_sph = cartesian_to_spherical(self.geo_equ_car)
        self.ra, self.dec, self.r = self.geo_equ_sph # just for ease of use

        if obs_loc is None:
            self.az, self.alt = None, None
        else:
            self.az, self.alt = radec_to_altaz(self.ra, self.dec, obs_loc, t)
        #=====================================================================
        

        # Phase angle and the elongation
        R = self.geo_ecl_sph[-1] # ehtemalan
        r = self.hel_ecl_sph[-1]
        s = self._sun.r

        self.elongation = arccos((s**2 + R**2 - r**2)/(2*s*R))*(180/pi)
        FV    = arccos((r**2 + R**2 - s**2)/(2*r*R))*(180/pi)
        self.FV = FV
        #self.phase =  (1 + cos(self.FV*rd))/2

        # Magnitude
        if self.name=='mercury':
            d0 = 6.74
            mag = -0.36 + 5*log10(r*R) + 0.027 * FV + 2.2E-13 * FV**6
        elif self.name=='venus':
            d0 = 16.92
            mag = -4.34 + 5*log10(r*R) + 0.013 * FV + 4.2E-7  * FV**3
        elif self.name=='mars':
            d0 = 9.32
            mag = -1.51 + 5*log10(r*R) + 0.016 * FV
        elif self.name=='jupiter':
            d0 = 191.01
            mag = -9.25 + 5*log10(r*R) + 0.014 * FV
        elif self.name=='saturn':
            d0 = 158.2
            ir = 28.06 # tilt rings to ecliptic
            Nr = 169.51 + 3.82E-5 * d # ascending node of plane of rings
            los = self.geo_ecl_sph[0] # Saturn's geocentric ecliptic longitude
            las = self.geo_ecl_sph[1] # Saturn's geocentric ecliptic latitude
            # B : tilt of Saturn's rings
            B = arcsin(sin(las*rd) * cos(ir*rd) - cos(las*rd) * sin(ir*rd) * sin((los-Nr)*rd))*(180/pi)
            ring_magn = -2.6 * sin(abs(B)*rd) + 1.2 * (sin(B*rd))**2
            mag = -9.0  + 5*log10(r*R) + 0.044 * FV + ring_magn
        elif self.name=='uranus':
            d0 = 63.95
            mag = -7.15 + 5*log10(r*R) + 0.001 * FV
        elif self.name=='neptune':
            d0 = 61.55
            mag = -6.90 + 5*log10(r*R) + 0.001 * FV
        else:
            mag = None
        self.mag = round(mag,2)
        self.diameter = d0 / self.r
        
