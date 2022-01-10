"""
Module plot
================
This module supplies several functions for plotting reasons.
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from collections.abc import Iterable
from .transform import radec_to_altaz, altaz_to_radec, angular_separation


def _equalize_scale(X,Y,Z, ax):
    """
    Equalize scales of x, y and z axis

    Arguments
    ---------
        X (np.array): x values
        Y (np.array): y values
        Z (np.array): z values
        ax (axes): axes object

    Returns
    -------
        matplotlib axes object
    """
    max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
    for xb, yb, zb in zip(Xb, Yb, Zb):
       ax.plot([xb], [yb], [zb], 'w')
    return ax

def plot_xyz(x, y, z, color, size, au=False):
    """
    Draw 3d cartesian plot of a body

    Arguments
    ---------
        x (np.array): x values
        y (np.array): y values
        z (np.array): z values
        color (str): color
        size (int): size
        au (bool): whether or not in AU unit; default False.

    Returns
    -------
        matplotlib axes object
    """
    fig = plt.figure(figsize=plt.figaspect(0.5)*1.2)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)
    ax.tick_params(axis='z', labelsize=8)
    au = 149597870700.0 if au else 1
    ax.scatter(x/au, y/au, z/au, c=color, s=size)
    ax = _equalize_scale(x/au, y/au, z/au, ax)
    return ax

def plot_altaz(az, alt, mag=None, size=None, color='k', alpha=1, marker='o', ax=None):
    """
    Plot positions of bodies based on altitude/azimuth coordinates

    Arguments
    ---------
        az (iter of float): azimuth values
        alt (iter of float): altitude values
        mag (iter of float): apparent magnitudes; default None.
        size (int): size; default None.
        color (str): color; default 'k'.
        alpha (float): alpha value (transparency), between 0 and 1; default 1.
        marker (str): marker shape; default 'o'.
        ax (axes): axes object; default None.

    Returns
    -------
        matplotlib axes object
    """

    if isinstance(az, Iterable):
        az = np.array(az)
        alt = np.array(alt)
        mag = np.array(mag) if mag is not None else None
    else:
        az = np.array([az])
        alt = np.array([alt])
        mag = None
        
    az  = az*(np.pi/180)

    if size is None:
        if mag is None:
            size = [5] * len(az)
        else:
            size = (0.1 + max(mag)-mag)**1.8
    
    if ax is None:
        fig = plt.figure()
        ax = fig.add_axes([0.1,0.1,0.8,0.8], polar=True)
        ax.set_theta_zero_location('N')
        ax.set_rlim(90, 0, 1)
        ax.set_yticks(np.arange(0,91,30))
        ax.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi, 5*np.pi/4, 6*np.pi/4, 7*np.pi/4])
        ax.set_xticklabels(['N','NE','E','SE','S','SW','W','NW'])
        ax.tick_params(axis=u'both', which=u'both',length=0)
    if matplotlib.__version__ < '3.0.0':
        alt = [90-i for i in alt]
    ax.scatter(az, alt, c=color, s=size, alpha=alpha, marker=marker)
    ax.grid(True, alpha=0.7)
    return ax

def plot_radec(ra, dec, mag=None, size=None):
    """
    Plot positions of bodies based on RA/Dec coordinates

    Arguments
    ---------
        ra (iter of float): Right Ascension values
        dec (iter of float): Declination values
        mag (iter of float): apparent magnitudes; default None.
        size (int): size; default None.

    Returns
    -------
        matplotlib axes object
    """
        
    if isinstance(ra, Iterable):
        ra = np.array(ra)
        dec = np.array(dec)
        mag = np.array(mag) if mag is not None else None
    else:
        ra = np.array([ra])
        dec = np.array([dec])
        mag = None

    ra = (np.mod(ra-180, 360)-180) * (np.pi/180)
    dec = dec * (np.pi/180)

    if size is None:
        if mag is None:
            size = [5] * len(ra)
        else:
            size = (0.1 + max(mag)-mag)**1.8

    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(111, projection="aitoff")
    ax.grid(True)
    ax.scatter(ra, dec, s=size, c='k')
    plt.subplots_adjust(top=0.95,bottom=0.0)
    plt.xlabel('RA')
    plt.ylabel('Dec')
    return ax



