from numpy import pi, sin, cos, tan, sqrt, arctan2, arcsin, arctan, arccos, log10
from datetime import datetime
import math

planets = ['mercury', 'venus', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune']

rd = pi/180
dg = 180/pi

def datetime_to_day(t):
    y = t.year
    m = t.month
    d = t.day
    UT = t.hour + t.minute/60 + t.second/3600
    return day(y,m,d,UT)

def day(y, m, D, UT):
    d = 367*y - 7 * ( y + (m+9)//12 ) // 4 + 275*m//9 + D - 730530
    d = d + UT/24.0
    return d

def getUT(d):
    if d < 0:
        UT = (d - int(d)+1) * 24
    else:
        UT = (d - int(d)) * 24
    return UT

def getUT_kh(d):
    return (d - int(d)) * 24

def floor(x):
    if x<0:
        return int(x)-1
    else:
        return int(x)
    
def rev(x):
    if x >= 360:
        while x >= 360:
            x = x - 360
    elif x < 0:
        while x < 0:
            x = 360 + x
    return x

def obl_ecl(d):
    """obliquity of the ecliptic (tilt of the Earth's axis of rotation)"""
    ecl = 23.4393 - 3.563E-7 * d
    return ecl

def obl_ecl_kh(d):
    # My own fit with JPL's JDTDB
    #ecl = 24.312212901355764 - 3.5606446971883423e-07 * d
    incl = 24.312603556024428 - 3.5623324732981026e-07 * d
    return incl

def getE(ec, m, dp=5):
    # http://www.jgiesen.de/kepler/kepler.html
    K = pi/180
    maxIter=30
    i=0
    delta = 10**-dp
    m = m/360.0
    m = 2 * pi * (m-floor(m))
    if ec<0.8:
        E=m
    else:
        E=pi
    F = E - ec*sin(m) - m
    while ((abs(F)>delta) and (i<maxIter)):
        E = E - F/(1-ec*cos(E))
        F = E - ec*sin(E) - m
        i = i + 1
    E = E/K
    return round(E*(10**dp)) / (10**dp)

def datetime_to_jd(t):
    year = t.year
    month = t.month
    t_d = t.day
    t_H = t.hour
    t_M = t.minute
    t_S = t.second
    t_MS = t.microsecond
    day = t_d + t_H/24 + t_M/(24*60) + t_S/(24*60*60) + t_MS/(24*60*60*1000000)
    
    if month == 1 or month == 2:
        yearp = year - 1
        monthp = month + 12
    else:
        yearp = year
        monthp = month
    
    if ((year < 1582) or
        (year == 1582 and month < 10) or
        (year == 1582 and month == 10 and day < 15)):
        # before start of Gregorian calendar
        B = 0
    else:
        # after start of Gregorian calendar
        A = math.trunc(yearp / 100.)
        B = 2 - A + math.trunc(A / 4.)
        
    if yearp < 0:
        C = math.trunc((365.25 * yearp) - 0.75)
    else:
        C = math.trunc(365.25 * yearp)
        
    D = math.trunc(30.6001 * (monthp + 1))
    jd = B + C + D + day + 1720994.5
    
    return jd

def jd_to_datetime(jd):

    jd = jd + 0.5
    F, I = math.modf(jd)
    I = int(I)
    A = math.trunc((I - 1867216.25)/36524.25)
    
    if I > 2299160:
        B = I + 1 + A - math.trunc(A / 4.)
    else:
        B = I
        
    C = B + 1524
    D = math.trunc((C - 122.1) / 365.25)
    E = math.trunc(365.25 * D)
    G = math.trunc((C - E) / 30.6001)
    day = C - E + F - math.trunc(30.6001 * G)
    
    if G < 13.5:
        month = G - 1
    else:
        month = G - 13
        
    if month > 2.5:
        year = D - 4716
    else:
        year = D - 4715

    d_ = day
    d = int(d_)
    h_ = (d_-d)*24
    h = int(h_)
    m_ = (h_-h)*60
    m = int(m_)
    s_ = (m_-m)*60
    s = int(s_)
    ms_ = (s_-s)*1000000
    ms = int(ms_)
    dt = datetime(year, month, d, h, m, s, ms)
        
    return dt
