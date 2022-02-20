**Author:** [Behrouz Safari](https://behrouzz.github.io/)<br/>
**License:** [MIT](https://opensource.org/licenses/MIT)<br/>

# solsys
*A python package for calculating positions of Solar System objects*


## Installation

Install the latest version of *solsys* from [PyPI](https://pypi.org/project/solsys/):

    pip install solsys

The only Requirement is *numpy*.

## Quick start

In the current version, there are three classes: *sun*, *moon* and *planet*. Let's get positions of the Sun, Moon and Venus for an observer in Strasbourg at this moment.

```python
from datetime import datetime
from solsys import sun, moon, planet

t = datetime.utcnow()
obs_loc = (7.7441, 48.5831)

s = sun(t, obs_loc)
m = moon(t, obs_loc)
p = planet('venus', t, obs_loc)

# ra and dec (GCRS) of Sun
print(s.ra, s.dec)

# azimuth and altitude of moon
print(m.az, m.alt)

# Geocentric equatorial cartesian coordinates of Venus
print(p.geo_equ_car)

# Apparent magnitude of Venus
print(p.mag)
```

See more at [astrodatascience.net](https://astrodatascience.net/)
