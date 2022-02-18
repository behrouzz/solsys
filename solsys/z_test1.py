from core import day, sun, planet
from datetime import datetime
        
t = datetime(1990, 4, 19, 0)
#d = day(1990, 4, 19, 0)
#d = day(2022, 2, 11, 10)
obs_loc = (2.2945, 48.8584)

s = sun(t, obs_loc)
p = planet('venus', t, obs_loc)#, epoch=2000)

print(s.ra, s.dec)
print(p.ra, p.dec)
print(s.az, s.alt)
print(p.az, p.alt)
