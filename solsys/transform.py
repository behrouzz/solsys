import numpy as np
from numpy import sin, cos, sqrt, arctan2, pi
from .utils import rev, obl_ecl, getE, rd, dg
from datetime import datetime

def mag(x):
    return np.linalg.norm(np.array(x))


def elements_to_ecliptic(name, N,i,w,a,e,M):
    E = getE(e, M, dp=15)
    xv = a * (cos(E*rd) - e)
    yv = a * (sqrt(1 - e**2) * sin(E*rd))
    r = sqrt(xv**2 + yv**2) # distance
    v = rev(arctan2(yv, xv)*dg) # true anomaly
    if name=='sun':
        lon = rev(w+v)
        x_ecl = r * cos(lon*rd)
        y_ecl = r * sin(lon*rd)
        z_ecl = 0.0
        return np.array([x_ecl, y_ecl, z_ecl])#, lon
    else:
        x_ecl = r * ( cos(N*rd) * cos((v+w)*rd) - sin(N*rd) * sin((v+w)*rd) * cos(i*rd) )
        y_ecl = r * ( sin(N*rd) * cos((v+w)*rd) + cos(N*rd) * sin((v+w)*rd) * cos(i*rd) )
        z_ecl = r * ( sin((v+w)*rd) * sin(i*rd) )
        return np.array([x_ecl, y_ecl, z_ecl])


def cartesian_to_spherical(xyz):
    x,y,z = xyz
    lon = rev(np.arctan2(y, x)*dg)
    lat = np.arctan2(z, np.sqrt(x**2 + y**2))*dg
    r = np.sqrt(x**2 + y**2 + z**2)
    return np.array([lon, lat, r])

def ecliptic_to_equatorial(ecl_xyz, d):
    x_ecl, y_ecl, z_ecl = ecl_xyz
    ecl = obl_ecl(d)
    x_equ = x_ecl
    y_equ = y_ecl * np.cos(ecl*rd) - z_ecl * np.sin(ecl*rd)
    z_equ = y_ecl * np.sin(ecl*rd) + z_ecl * np.cos(ecl*rd)
    return np.array([x_equ, y_equ, z_equ])

def equatorial_to_ecliptic(spherical):
    ra, dec, r = spherical
    lon, lat = equ2ecl(ra, dec)
    return np.array([lon, lat, r])

def equcar_to_eclcar(xyz):
    equ_sph = cartesian_to_spherical(xyz)
    ecl_sph = equatorial_to_ecliptic(equ_sph)
    ecl_car = spherical_to_cartesian(ecl_sph)
    return ecl_car

def equ2ecl(ra, dec):
    ra, dec = ra*(pi/180), dec*(pi/180)
    epsilon = (23 + 26/60 + 21.448/3600)*pi/180
    sb = sin(dec)*cos(epsilon) - cos(dec)*sin(epsilon)*sin(ra)
    cbcl = cos(dec)*cos(ra)
    cbsl = sin(dec)*sin(epsilon) + cos(dec)*cos(epsilon)*sin(ra)
    lon = rev(arctan2(cbsl,cbcl)*(180/pi))
    r = sqrt(cbsl**2 + cbcl**2)
    lat = arctan2(sb,r)*(180/pi)
    return lon, lat

def spherical_to_cartesian(spherical):
    lon, lat, r = spherical
    x = r * np.cos(lat*rd) * np.cos(lon*rd)
    y = r * np.cos(lat*rd) * np.sin(lon*rd)
    z = r * np.sin(lat*rd)
    return np.array([x, y, z])


def radec_to_altaz(ra, dec, obs_loc, t):
    lon, lat = obs_loc

    J2000 = datetime(2000,1,1,12)
    d = (t - J2000).total_seconds() / 86400 #day offset

    UT = t.hour + t.minute/60 + t.second/3600
    LST = (100.46 + 0.985647 * d + lon + 15*UT + 360) % 360
    ha = (LST - ra + 360) % 360
    
    x = np.cos(ha*rd) * np.cos(dec*rd)
    y = np.sin(ha*rd) * np.cos(dec*rd)
    z = np.sin(dec*rd)
    xhor = x*np.cos((90-lat)*rd) - z*np.sin((90-lat)*rd)
    yhor = y
    zhor = x*np.sin((90-lat)*rd) + z*np.cos((90-lat)*rd)
    az = np.arctan2(yhor, xhor)*dg + 180
    alt = np.arcsin(zhor)*dg
    return az, alt


def state_to_element(r, v, mu=1.32712440018e+20):
    """
    N (rad) : longitude of the ascending node
    i (rad) : inclination to the ecliptic
    w (rad) : argument of perihelion
    a (m)   : semi-major axis, or mean distance from Sun
    e       : eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M (rad) : mean anomaly
    """
    h = np.cross(r, v)

    # eccentricity vector
    ev = (np.cross(v,h)/mu) - (r/mag(r))
    e = mag(ev) # orbit eccentricity

    # vector pointing towards the ascending node
    n = np.cross(np.array([0,0,1]), h)

    # true anomaly (in radian)
    tmp = np.inner(ev,r)/(e*mag(r))

    if np.inner(r,v) >= 0:
        vv = np.arccos(tmp)
    else:
        vv = 2*np.pi - np.arccos(tmp)

    # inclination (in radian)
    i = np.arccos(h[-1]/mag(h))

    # eccentric anomaly
    E = 2 * np.arctan(np.tan(vv/2)/np.sqrt((1+e)/(1-e)))

    # longitude of the ascending node (in radian)
    if n[1] >= 0:
        N = np.arccos(n[0]/mag(n))
    else:
        N = 2*pi - np.arccos(n[0]/mag(n))

    # argument of perihelion
    tmp = np.inner(n, ev) / (mag(n)*e)
    if ev[-1] >= 0:
        w = np.arccos(tmp)
    else:
        w = 2*np.pi - np.arccos(tmp)

    # mean anomaly
    M = E - e*np.sin(E)

    # semi-major axis
    a = 1 / ((2/mag(r)) - (mag(v)**2 / mu))

    return N,i,w,a,e,M


def state_to_kepler(r_, v_, mu=0.00029591220819207774):
    """
    (All distances in AU, times in d, angles in radian)

    N   : longitude of the ascending node
    i   : inclination to the ecliptic
    w   : argument of perihelion
    a   : semi-major axis
    e   : eccentricity
    M   : mean anomaly
    TA  : true anomaly
    L   : mean longitude
    T   : period
    tsp : time since periapse passage
    EA  : eccentric anomaly
    n   : Mean angular motion
    """
    
    # specific angular momentum
    h_ = np.cross(r_, v_) # = L_/m = (r_Xp_)/m = (r_Xm*v_)/m = r_ X v_
    h = mag(h_)

    # radius and velocity
    r = mag(r_)
    v = mag(v_)

    # specific energy
    E = ((v**2)/2) - (mu/r)

    # semi-major axis
    a = -mu/(2*E)

    # eccentricity
    e = np.sqrt( 1 - (h**2/(a*mu)) )

    # inclination (0 to 180)
    i = np.arccos(h_[-1]/h)

    # right ascension of ascending node (0 to 360)
    N = np.arctan2(h_[0], -h_[1])

    # argument of latitude (w+v) (0 to 360)
    if i!=0:
        w_v = np.arctan2( (r_[-1]/np.sin(i)) , (r_[0]*np.cos(N)+r_[1]*np.sin(N)) )
    else:
        w_v = np.arctan2( np.inf , (r_[0]*np.cos(N)+r_[1]*np.sin(N)) )

    # true anomaly (0 to 360)
    p = a*(1-e**2)
    TA = np.arctan2( np.sqrt(p/mu)*np.dot(v_,r_) , p-r )

    # argument of periapse
    w = w_v - TA

    # eccentric anomaly (0 to 360)
    EA = 2 * np.arctan( np.sqrt((1-e)/(1+e)) * np.tan(TA/2) )

    # mean anomaly
    M = EA - e*np.sin(EA)

    # Mean angular motion
    n = np.sqrt(mu/a**3)

    # mean longitude
    L = M + w

    # period
    T = 2*np.pi * np.sqrt((a**3)/mu)

    # time since periapse passage
    tsp = np.sqrt((a**3)/mu) * M
    dc = {'N':N, 'i':i, 'w':w, 'a':a, 'e':e, 'M':M,
          'TA':TA, 'L':L, 'T':T, 'tsp':tsp, 'EA':EA, 'n':n}
    return dc
