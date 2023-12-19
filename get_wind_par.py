# Core packages
import os
import numpy as np

# Astro packages
from astropy.constants import R_sun
from astropy.constants import L_sun
from astropy.constants import sigma_sb

print('\nUnit assumptions:')
print('R_sun =', R_sun.cgs.value, 'in', R_sun.cgs.unit) # M.Lorenzo R_sun = 6.6599e10 cm
print('L_sun =', L_sun.cgs.value, 'in', L_sun.cgs.unit) # M.Lorenzo L_sun = 3.839e33 erg/s
print('Stefan-Boltzmann constant =', sigma_sb.cgs.value, 'in', sigma_sb.cgs.unit)
print('Bolometric magnitude of the Sun =', 4.74)
print('\n')

def get_Mdot_R_vinf_Miguel(teff, lgf, logg, logq):
    """
    Calculate the radius of a star from its effective temperature, surface
    gravity, luminosity and mass loss rate.
    
    Parameters
    ----------
    teff : float
        Effective temperature of the star in K.
    
    lgf : float
        Logarithm of the surface gravity of the star.
    
    logg : float
        Logarithm of the surface gravity of the star.
    
    logq : float
        Logarithm of the wind-strength parameter.
    
    Returns
    -------
    Mdot : float
        Mass loss rate of the star in ???
    
    R : float
        Radius of the star in ???
    
    v_inf : float
        Terminal velocity of the wind of the star in ???
    """
    
    # Bolometric magnitude of the Sun from Cox 2000
    Mbol_sun = 4.74

    # Calculate R (from Urbaneja+2017)
    a = 3.20
    a_low = 8.34
    b = -7.90
    lgf_break = 1.3

    if lgf >= lgf_break:
        Mbol = a * (lgf - 1.5) + b
    else:
        b_break = a * (lgf_break - 1.5) + b
        Mbol = a_low * (lgf - lgf_break) + b_break

    logL = -(1 / 2.5) * (Mbol - Mbol_sun)
    L = 10 ** logL * L_sun.cgs.value
    R = np.sqrt(L / (4 * np.pi * sigma_sb.cgs.value * teff ** 4))
    R_rsun = R / R_sun.cgs.value

    # Calculate wind parameters:
    v_esc = np.sqrt(2 * 10 ** logg * (R)) * 1e-5
    t = teff * 1e-4
    if t <= 1.7:
        v_inf = 1.25 * v_esc

    elif t <= 2.3:
        v_inf = (0.2917 * (teff * 10 ** -3) - 3.7083) * v_esc
    else:
        v_inf = 3. * v_esc

    # Q = Mdot/R*vinf)**1.5
    Mdot = 10 ** logq * (R_rsun * v_inf) ** 1.5

    return Mdot, R_rsun, v_inf