from core import day, sun, planet
        

#d = day(1990, 4, 19, 0)
d = day(2022, 2, 11, 10)

s = sun(d)
p = planet('venus', d)#, epoch=2000)

print(s.name, ':', (s.ra, s.dec))
print(p.name, ':', (p.ra, p.dec))
print(p.r)
