import numpy as np
from numpy import pi, sin, cos, tan, sqrt, arctan2, arcsin, arctan, arccos, log10
from orbital_elements import elements, jupiter_oe, saturn_oe, uranus_oe
from utils import *
from transform import cartesian_to_spherical, spherical_to_cartesian, ecliptic_to_equatorial, elements_to_ecliptic, radec_to_altaz
from correction import moon_perts, jupiter_lon_perts, saturn_lon_perts, saturn_lat_perts, uranus_lon_perts

class sun:
    """
    Sun positional parameters
    
    Parameters
    ----------
        d (datetime): time of observation
        obs_loc : tuple of observer location (longtitude, latitude)

    Attributes
    ----------
        ecl_sph  : ecliptic spherical coordinates (lon, lat, r)
        ecl_car  : ecliptic cartesian coordinates (x, y, z)
        equ_car  : equatorial cartesian coordinates (x, y, z)
        equ_sph  : equatorial spherical coordinates (ra, dec, r)
        alt      : azimuth
        az       : altitude
        L   : Sun's mean longitude
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
        obs_loc : tuple of observer location (longtitude, latitude)

    Attributes
    ----------
        elements : dictionary of orbital elements
        v       : true anomaly
        E       : eccentric anomaly
        L       : mean longitude
        ra      : Right Ascension (GCRS or topocentric if obs_loc is provided)
        dec     : Declination (GCRS or topocentric if obs_loc is provided)
        r       : distance to Earth
        ecl_lon : ecliptic longitude (GCRS)
        ecl_lat : ecliptic latitude (GCRS)
        x_ecl   : ecliptic x (GCRS)
        z_ecl   : ecliptic y (GCRS)
        y_ecl   : ecliptic z (GCRS)
        x_equ   : equatorial x (GCRS)
        y_equ   : equatorial y (GCRS)
        z_equ   : equatorial z (GCRS)
        elongation : elongation
        FV         : phase angle
    """
        
    def __init__(self, t, obs_loc=None, epoch=None):
        self.name = 'moon'
        d = datetime_to_day(t)
        ecl = obl_ecl(d)
        #self.obs_loc = obs_loc
        self._sun = sun(t)
        N, i, w, a, e, M = elements(self.name, d)
        self.elements = {'N':N, 'i':i, 'w':w, 'a':a, 'e':e, 'M':M}
        #x_ecl, y_ecl, z_ecl = elements_to_ecliptic('moon', N,i,w,a,e,M)
        x, y, z = elements_to_ecliptic('moon', N,i,w,a,e,M)
        #ecl_lon, ecl_lat, ecl_r = cartesian_to_spherical(x_ecl, y_ecl, z_ecl)
        ecl_lon, ecl_lat, ecl_r = cartesian_to_spherical(x, y, z)

        self.L = rev(N+w+M) # Moon's mean longitude

        # CONSIDERING Perturbations
        # =========================
        Ls = self._sun.L
        Ms = self._sun.elements['M']
        
        Lm = self.L
        Mm = self.elements['M']
        Nm = self.elements['N']

        D = Lm - Ls # Moon's mean elongation
        F = Lm - Nm # Moon's argument of latitude

        pert_lon, pert_lat, pert_r = moon_perts(Ls, Ms, Lm, Mm, Nm, D, F)

        # Add this to the ecliptic positions we earlier computed:
        self.ecl_lon = ecl_lon + pert_lon
        self.ecl_lat = ecl_lat + pert_lat
        self.r       = ecl_r   + pert_r

        r_ = 1

        x_ecl = r_ * cos(self.ecl_lon*rd) * cos(self.ecl_lat*rd)
        y_ecl = r_ * sin(self.ecl_lon*rd) * cos(self.ecl_lat*rd)
        z_ecl = r_ * sin(self.ecl_lat*rd)

        x_equ, y_equ, z_equ = ecliptic_to_equatorial(x_ecl, y_ecl, z_ecl, d)

        self.ra, self.dec, _ = cartesian_to_spherical(x_equ, y_equ, z_equ)

        # kh: in qesmato khodam anjam dadam (motmaen nistam)
        # badan check shavad
        self.x_ecl = x_ecl * self.r
        self.y_ecl = y_ecl * self.r
        self.z_ecl = z_ecl * self.r

        self.x_equ, self.y_equ, self.z_equ = ecliptic_to_equatorial(self.x_ecl, self.y_ecl, self.z_ecl, d)


        if obs_loc is not None:
            self.ra, self.dec = self.topocentric_correction(obs_loc)
        
        self.elongation = arccos( cos((self._sun.ecl_sph[0]-self.ecl_lon)*rd) * cos(self.ecl_lat*rd) )*(180/pi)
        self.FV = 180 - self.elongation

    def topocentric_correction(self, obs_loc):
        LON, LAT = obs_loc
        mpar = arcsin(1/self.r)*(180/pi) # moon parallax
        gclat = LAT - 0.1924 * sin(2*LAT*rd)
        rho   = 0.99833 + 0.00167 * cos(2*LAT*rd)

        UT = getUT(d)

        GMST0 = rev(self._sun.L + 180) / 15
        LST = GMST0 + UT + LON/15 # in hours
        LST_deg = LST * 15
        HA = rev(LST_deg - self.ra)

        # auxiliary angle
        g = arctan( tan(gclat*rd) / cos(HA*rd) )*(180/pi)

        topRA  = self.ra  - mpar * rho * cos(gclat*rd) * sin(HA*rd) / cos(self.dec*rd)
        topDEC = self.dec - mpar * rho * sin(gclat*rd) * sin((g-self.dec)*rd) / sin(g*rd)
        return topRA, topDEC
        

class planet:
    """
    Planets positional parameters
    
    Parameters
    ----------
        d (datetime): time of observation
        name (str) : name of the planet
        obs_loc : tuple of observer location (longtitude, latitude)

    Attributes
    ----------
        
        hel_ecl_sph  : Heliocentric ecliptic spherical coordinates (lon, lat, r)
        hel_ecl_car  : Heliocentric ecliptic cartesian coordinates (x, y, z)
        geo_ecl_car  : Geocentric ecliptic cartesian coordinates (x, y, z)
        geo_ecl_sph  : Geocentric ecliptic spherical coordinates (lon, lat, r)
        geo_equ_car  : Geocentric equatorial cartesian coordinates (x, y, z)
        geo_equ_sph  : Geocentric equatorial spherical coordinates (ra, dec, r)
        
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
        
