from core import sun, planet
from datetime import datetime
from utils import planets
        
#t = datetime(1990, 4, 19, 0)
t = datetime.utcnow()
obs_loc = (15, 60)

for i in planets:
    p = planet(i, t)#, obs_loc)
    print(i, ':', p.mag)
