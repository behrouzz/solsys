from numpy import pi, sin, cos, tan, sqrt, arctan2, arcsin, arctan, arccos, log10
from orbital_elements import elements
from utils import *
from transform import cartesian_to_spherical, spherical_to_cartesian, ecliptic_to_equatorial, elements_to_ecliptic
from correction import moon_perts, jupiter_lon_perts, saturn_lon_perts, saturn_lat_perts, uranus_lon_perts

class sun:
    """
    L   : Sun's mean longitude
    lon : Sun's true longitude
    """
    def __init__(self, d, obs_loc=None):
        self.name = 'sun'
        ecl = obl_ecl(d)
        self.d = d
        self.N, self.i, self.w, self.a, self.e, self.M = elements(self.name, d)
        self.L = rev(self.w + self.M)
        self.x_ecl, self.y_ecl, self.z_ecl, self.lon = elements_to_ecliptic('sun', self.N, self.i, self.w, self.a, self.e, self.M)
        self.x_equ, self.y_equ, self.z_equ = ecliptic_to_equatorial(self.x_ecl, self.y_ecl, self.z_ecl, d)
        self.ra, self.dec, self.r = cartesian_to_spherical(self.x_equ, self.y_equ, self.z_equ)

        if obs_loc is None:
            self.az, self.alt = None, None
        else:
            self.az, self.alt = self.altaz(obs_loc)
        

    def altaz(self, obs_loc):
        LON, LAT = obs_loc
        UT = getUT(self.d)
        
        GMST0 = rev(self.L + 180) / 15
        SIDTIME = GMST0 + UT + LON/15 # in hours
        RA_hour = self.ra/15 # convert RA from degrees to hours
        HA_hour = SIDTIME - RA_hour
        HA = HA_hour *15 # converts HA from hours to degrees

        x = cos(HA *rd) * cos(self.dec*rd)
        y = sin(HA *rd) * cos(self.dec*rd)
        z = sin(self.dec*rd)

        xhor = x * sin(LAT*rd) - z * cos(LAT*rd)
        yhor = y
        zhor = x * cos(LAT*rd) + z * sin(LAT*rd)

        az  = arctan2(yhor, xhor)*(180/pi) + 180
        alt = arcsin(zhor)*(180/pi)
        return az, alt
        

class moon:
    """
    Moon positional parameters
    
    Parameters
    ----------
        d (datetime): time of observation
        obs_loc : tuple of observer location (longtitude, latitude)

    Attributes
    ----------
        N       : longitude of the ascending node
        i       : inclination to the ecliptic
        w       : argument of perihelion
        a       : semi-major axis, or mean distance from Sun
        e       : eccentricity
        M       : mean anomaly
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
        
    def __init__(self, d, obs_loc=None):
        self.name = 'moon'
        ecl = obl_ecl(d)
        #self.obs_loc = obs_loc
        self._sun = sun(d)
        self.N, self.i, self.w, self.a, self.e, self.M = elements(self.name, d)
        x_ecl, y_ecl, z_ecl = elements_to_ecliptic('moon', self.N, self.i, self.w, self.a, self.e, self.M)
        ecl_lon, ecl_lat, ecl_r = cartesian_to_spherical(x_ecl, y_ecl, z_ecl)

        # CONSIDERING Perturbations
        # =========================

        self.L = rev(self.N + self.w + self.M) # Moon's mean longitude

        Ls = self._sun.L
        Ms = self._sun.M
        
        Lm = self.L
        Mm = self.M
        Nm = self.N

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

        self.elongation = arccos( cos((self._sun.lon-self.ecl_lon)*rd) * cos(self.ecl_lat*rd) )*(180/pi)
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
        
        lon_h_ecl  : Heliocentric ecliptic longitude
        lat_h_ecl  : Heliocentric ecliptic latitude
        r_h_ecl    : Heliocentric ecliptic distance
        xh_ecl     : Heliocentric ecliptic x
        yh_ecl     : Heliocentric ecliptic y
        zh_ecl     : Heliocentric ecliptic z
        lon_g_ecl  : Geocentric ecliptic longitude
        lat_g_ecl  : Geocentric ecliptic latitude
        r_g_ecl    : Geocentric ecliptic distance
        xg_ecl     : Geocentric ecliptic x
        yg_ecl     : Geocentric ecliptic y
        zg_ecl     : Geocentric ecliptic z
        ra         : Right Ascension (GCRS)
        dec        : Declination (GCRS)
        r          : Distance to Earth
        x          : equatorial x (GCRS)
        y          : equatorial y (GCRS)
        z          : equatorial z (GCRS)
        elongation : elongation
        FV         : phase angle
        mag        : Apparent magnitude
        diameter   : Apparent diameter
    """
    def __init__(self, name, d, obs_loc=None, epoch=None):
        self.name = name.lower()
        #self.obs_loc = obs_loc
        ecl = obl_ecl(d)
        self._sun = sun(d)
        N,i,w,a,e,M = elements(self.name, d)
        
        xh_ecl, yh_ecl, zh_ecl = elements_to_ecliptic(self.name, N,i,w,a,e,M)

        # spherical heliocentric ecliptic
        self.lon_h_ecl, self.lat_h_ecl, self.r_h_ecl = \
                        cartesian_to_spherical(xh_ecl, yh_ecl, zh_ecl)

        # Correcting perturbations of Jupiter, Saturn and Uranus
        if self.name in ['jupiter', 'saturn', 'uranus']:
            Mj = jupiter_oe(d)[-1]
            Ms = saturn_oe(d)[-1]
            Mu = uranus_oe(d)[-1]
            if self.name=='jupiter':
                self.lon_h_ecl = self.lon_h_ecl + jupiter_lon_perts(Mj, Ms, Mu)                
            elif self.name=='saturn':
                self.lon_h_ecl = self.lon_h_ecl + saturn_lon_perts(Mj, Ms, Mu)
                self.lat_h_ecl = self.lat_h_ecl + saturn_lat_perts(Mj, Ms, Mu)
            elif self.name=='uranus':
                self.lon_h_ecl = self.lon_h_ecl + uranus_lon_perts(Mj, Ms, Mu)
                
        # Precession
        if epoch is not None:
            self.lon_h_ecl = self.lon_h_ecl + 3.82394E-5 * (365.2422 * (epoch-2000) - d)

        # kh: hypatie
        self.xh_ecl, self.yh_ecl, self.zh_ecl = \
                     spherical_to_cartesian(self.lon_h_ecl, self.lat_h_ecl, self.r_h_ecl)

        self.xg_ecl = self._sun.x_ecl + self.xh_ecl
        self.yg_ecl = self._sun.y_ecl + self.yh_ecl
        self.zg_ecl = self._sun.z_ecl + self.zh_ecl

        self.lon_g_ecl, self.lat_g_ecl, self.r_g_ecl = \
                        cartesian_to_spherical(self.xg_ecl, self.yg_ecl, self.zg_ecl)

        # equatorial rectangular geocentric
        self.x, self.y, self.z = ecliptic_to_equatorial(self.xg_ecl, self.yg_ecl, self.zg_ecl, d)

        # converet to spherical coordinates (radec)
        self.ra, self.dec, self.r = cartesian_to_spherical(self.x, self.y, self.z)

        # Phase angle and the elongation
        R = self.r
        r = self.r_h_ecl
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
            los = self.lon_g_ecl # Saturn's geocentric ecliptic longitude
            las = self.lat_g_ecl # Saturn's geocentric ecliptic latitude
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
        self.mag = mag
        self.diameter = d0 / self.r
