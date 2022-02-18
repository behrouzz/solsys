from numpy import sin, cos, pi

rd = pi/180
dg = 180/pi


def moon_perts(Ls, Ms, Lm, Mm, Nm, D, F):
    """Perturbations in Moon's coordinates"""
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
    return pert_lon, pert_lat, pert_r


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
