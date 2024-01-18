# Core packages
import os
import numpy as np

# Astro packages
from astropy.constants import c
from astropy.constants import R_sun
from astropy.constants import L_sun
from astropy.constants import sigma_sb

print('\nUnit assumptions:')
print('c =', c.cgs.value, 'in', c.cgs.unit) # M.A.Urbaneja uses c = 2.99e10 cm/s
print('R_sun =', R_sun.cgs.value, 'in', R_sun.cgs.unit) # M.A.Urbaneja uses R_sun = 6.6599e10 cm
print('L_sun =', L_sun.cgs.value, 'in', L_sun.cgs.unit)
print('Stefan-Boltzmann constant =', sigma_sb.cgs.value, 'in', sigma_sb.cgs.unit)
print('Bolometric magnitude of the Sun =', 4.74)
print('\n')

def get_Mdot_R_vinf(teff, logg, logq, prescription='Urbaneja', sgs=True, z=1.00, yhe=0.1):
    """
    Calculate the radius of a star from its effective temperature, surface
    gravity, luminosity and mass loss rate.
    
    Parameters
    ----------
    teff : float
        Effective temperature of the star in K.
    
    logg : float
        Logarithm of the surface gravity of the star.
    
    logq : float
        Logarithm of the wind-strength parameter.
    
    prescription : str
        Prescription to calculate the mass loss rate, radius,
        terminal velocity and escape velocity of the star.
        'Urbaneja' (default), 'Urbaneja2017' or 'Sergio'.
    
    sgs : bool, optional
        [Only for prescription='Urbaneja']
        If True, use the prescription for supergiants (default).
        If False, use the prescription for non-supergiants.
        Else, use the prescription from Kudritzki et al. 2008.
    
    z : float
        Metallicity of the star. Default is 1.00.
    
    yhe : float
        Helium abundance of the star. Default is 0.1.
    
    Returns
    -------
    Mdot : float
        Mass loss rate of the star in solar masses per year.
    
    R : float
        Radius of the star in solar radii.
    
    v_inf : float
        Terminal velocity of the wind of the star in km/s.
        
    v_esc : float
        Escape velocity of the star in km/s.
    """
    
    # From: 'Nominal values for selected solar and planetary quantities: IAU 2015 resolution B3'
    Mbol_sun = 4.74 # Bolometric magnitude of the Sun (M.A.Urbaneja uses 4.75)
    Teff_sun = 5772 # Effective temperature of the Sun in K (M.A.Urbaneja uses 5779.1389)

    # Calculate loggf from teff and logg:
    loggf = logg - 4 * np.log10(teff * 1e-4)

    if prescription == 'Urbaneja2017':

        # Calculate Mbol from loggf (Urbaneja et al. 2017):
        a = 3.20
        a_low = 8.34
        b = -7.90
        loggf_break = 1.3

        if loggf >= loggf_break:
            Mbol = a * (loggf - 1.5) + b
        else:
            b_break = a * (loggf_break - 1.5) + b
            Mbol = a_low * (loggf - loggf_break) + b_break

        # Calculate R:
        logLLsol = - 0.4 * (Mbol - Mbol_sun)
        L = 10 ** logLLsol * L_sun.cgs.value # Luminosity star
        R = np.sqrt(L / (4 * np.pi * sigma_sb.cgs.value * teff ** 4)) # Stefan-Boltzmann law
        R_rsun = R / R_sun.cgs.value # Radius star in solar radii

        # Calculate escape velocity (v_esc):
        v_esc = np.sqrt(2 * 10 ** logg * R) * 1e-5

        # Calculate terminal velocity (v_inf):
        t = teff * 1e-4
        if t <= 1.7:
            v_inf = 1.25 * v_esc
        elif t <= 2.3:
            v_inf = (0.2917 * (teff * 1e-3) - 3.7083) * v_esc
        else:
            v_inf = 3 * v_esc

        # Q = Mdot/(R*vinf)**1.5
        Mdot = 10 ** logq * (R_rsun * v_inf) ** 1.5

        return Mdot, R_rsun, v_inf, v_esc


    elif prescription == 'Urbaneja':

        # Calculate Mbol from loggf from FGLR with different prescriptions:
        # Kudritzki, Urbaneja & Rix, 2020, ApJ, 890, 28 
        if sgs == True and loggf >= 1.29:
            Mbol = 3.33 * (loggf - 1.5) - 7.89
        elif sgs == True and loggf < 1.29:
            Mbol = 8.07 * (loggf - 1.29) + 3.33 * (1.29 - 1.5) - 7.89
        elif sgs == False:
            Mbol = 3.76 * (loggf - 3.0) - 2.98 # new non-evolved eFGLR
        else:
            # Kudritzki et al. 2008
            Mbol = 3.41 * (loggf - 1.5) - 8.02

        # Calculate R:
        logLLsol = - 0.4 * (Mbol - Mbol_sun)
        R_rsun = 10 ** (0.5 * (logLLsol - 4 * np.log10(teff / Teff_sun))) # Stefan-Boltzmann law
        # NOTE: Here instead of using the L_sun, it uses the Teff_sun

        # Calculate escape velocity (v_esc):
        solar_metallicity = 0.02
        x = (1 - z * solar_metallicity) / (1 + 4 * yhe) # H mass fraction
        gamma = (sigma_sb.cgs.value / c.cgs.value) * 0.20 * (1 + x) * teff ** 4 / 10 ** logg
        v_esc = np.sqrt(2 * 10 ** logg * (R_rsun * R_sun.cgs.value) * (1 - gamma)) * 1e-5

        # Calculate terminal velocity (v_inf) using Achim's vinfty/vesc prescription:
        t = teff * 1e-4
        if t <= 1.7:
            v_inf = 1.25 * v_esc
        elif t <= 2.3:
            v_inf = (0.2917 * (teff * 1e-3) - 3.7083) * v_esc
        else:
            v_inf = 3 * v_esc

        # Q = Mdot/(R*vinf)**1.5
        Mdot = 10 ** logq * (R_rsun * v_inf) ** 1.5

        return Mdot, R_rsun, v_inf, v_esc


    elif prescription == 'Sergio':
        
        print('To be implemented')
        print('Check file: crea_chainfil_allgrid_hhex.pro')
        
        return None, None, None, None