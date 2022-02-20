from numpy import array, pi, sin, cos, tan, arcsin, arctan
from .utils import rev, getUT

rd = pi/180
dg = 180/pi






def moon_perts(geo_ecl_sph, sun_elem, moon_elem):
    """Perturbations in Moon's coordinates"""
    Ls = rev(sun_elem['w']+sun_elem['M'])# Sun's mean longitude
    Ms = sun_elem['M']
    Lm = rev(moon_elem['N']+moon_elem['w']+moon_elem['M'])# Moon's mean longitude
    Mm = moon_elem['M']
    Nm = moon_elem['N']
    D = Lm - Ls # Moon's mean elongation
    F = Lm - Nm # Moon's argument of latitude
    
    # Perturbations in longitude (degrees)
    pert_lon = - 1.274 * sin((Mm-2*D)*rd) + 0.658 * sin(2*D*rd) \
               - 0.186 * sin(Ms*rd) - 0.059 * sin((2*Mm-2*D)*rd)\
               - 0.057 * sin((Mm-2*D+Ms)*rd) + 0.053 * sin((Mm+2*D)*rd)\
               + 0.046 * sin((2*D-Ms)*rd) + 0.041 * sin((Mm-Ms)*rd)\
               - 0.035 * sin(D*rd) - 0.031 * sin((Mm+Ms)*rd)\
               - 0.015 * sin((2*F-2*D)*rd) + 0.011 * sin((Mm-4*D)*rd)

    # Perturbations in latitude (degrees):
    pert_lat = - 0.173 * sin((F-2*D)*rd) - 0.055 * sin((Mm-F-2*D)*rd)\
               - 0.046 * sin((Mm+F-2*D)*rd) + 0.033 * sin((F+2*D)*rd)\
               + 0.017 * sin((2*Mm+F)*rd)

    # Perturbations in lunar distance (Earth radii):
    pert_r = -0.58 * cos((Mm-2*D)*rd) - 0.46 * cos(2*D*rd)

    # Add to the ecliptic positions computed earlier:
    geo_ecl_sph[0] = geo_ecl_sph[0] + pert_lon
    geo_ecl_sph[1] = geo_ecl_sph[1] + pert_lat
    geo_ecl_sph[2] = geo_ecl_sph[2] + pert_r
    
    return geo_ecl_sph



def topocentric(obs_loc, Lsun, geo_equ_sph, d):
    LON, LAT = obs_loc
    ra, dec, r = geo_equ_sph
    mpar = arcsin(1/r)*(180/pi) # moon parallax
    gclat = LAT - 0.1924 * sin(2*LAT*rd)
    rho   = 0.99833 + 0.00167 * cos(2*LAT*rd)

    UT = getUT(d)

    GMST0 = rev(Lsun + 180) / 15
    LST = GMST0 + UT + LON/15 # in hours
    LST_deg = LST * 15
    HA = rev(LST_deg - ra)

    # auxiliary angle
    g = arctan( tan(gclat*rd) / cos(HA*rd) )*(180/pi)

    topRA  = ra  - mpar * rho * cos(gclat*rd) * sin(HA*rd) / cos(dec*rd)
    topDEC = dec - mpar * rho * sin(gclat*rd) * sin((g-dec)*rd) / sin(g*rd)
    return array([topRA, topDEC, r])



def jupiter_lon_perts(Mj, Ms, Mu):
    pert_lon =  -0.332 * sin((2*Mj-5*Ms-67.6)*rd)\
               - 0.056 * sin((2*Mj-2*Ms+21)*rd) \
               + 0.042 * sin((3*Mj-5*Ms+21)*rd) \
               - 0.036 * sin((Mj-2*Ms)*rd) \
               + 0.022 * cos((Mj-Ms)*rd) \
               + 0.023 * sin((2*Mj-3*Ms+52)*rd) \
               - 0.016 * sin((Mj-5*Ms-69)*rd)
    return pert_lon


def saturn_lon_perts(Mj, Ms, Mu):
    pert_lon = +0.812 * sin((2*Mj-5*Ms-67.6)*rd) \
               - 0.229 * cos((2*Mj-4*Ms-2)*rd) \
               + 0.119 * sin((Mj-2*Ms-3)*rd) \
               + 0.046 * sin((2*Mj-6*Ms-69)*rd) \
               + 0.014 * sin((Mj-3*Ms+32)*rd)
    return pert_lon


def saturn_lat_perts(Mj, Ms, Mu):
    pert_lat =  -0.020 * cos((2*Mj-4*Ms-2)*rd) \
               + 0.018 * sin((2*Mj-6*Ms-49)*rd)
    return pert_lat

def uranus_lon_perts(Mj, Ms, Mu):
    pert_lon =  +0.040 * sin((Ms-2*Mu+6)*rd) \
               + 0.035 * sin((Ms-3*Mu+33)*rd) \
               - 0.015 * sin((Mj-Mu+20)*rd)
    return pert_lon
