from core import sun, planet, moon
from datetime import datetime
        
#t = datetime(1990, 4, 19, 0)
t = datetime(2022, 2, 25, 0)
#t = datetime.utcnow()
obs_loc = (15, 60)

m = moon(t, obs_loc)


print('geo_ecl_car:', m.geo_ecl_car)
print('geo_ecl_sph:', m.geo_ecl_sph)
print('geo_equ_car:', m.geo_equ_car)
print('geo_equ_sph:', m.geo_equ_sph)
print('ra, dec, r :', (m.ra, m.dec, m.r))
print('elongation :', m.elongation)
print('FV         :', m.FV)
