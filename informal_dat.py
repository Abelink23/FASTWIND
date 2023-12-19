# Core packages
import os
import numpy as np

from init import *
from get_wind_par import *

def write_indat(dir, teff, lgf, logq, beta, he, c, n, o, si, mg, z, micro, micro_fw, heion):
    '''
    Write the INDAT.DAT files for FASTWIND from a set of input parameters.
    
    NOTE: This function is intended to be use alongside MAUI in a sense that grids are defined using lgf and not logg directly.
    
    Parameters
    ----------
    dir : str
        Path to the directory where the INDAT.DAT file will be saved.
    
    teff : float
        Effective temperature in K.
    
    lgf : float
        Logarithm of the surface gravity.
    
    logq : float
        Logarithm of the wind-strength parameter.
    
    beta : float
        The exponent of the wind-velocity law.
    
    he : float
        Helium abundance.
    
    c : float / bool
        Carbon abundance.
    
    n : float / bool
        Nitrogen abundance.
    
    o : float / bool
        Oxygen abundance.
    
    si : float / bool
        Silicon abundance.
    
    mg : float / bool
        Magnesium abundance.
    
    z : float
        Metallicity.
    
    micro : float
        Microturbulent velocity.
    
    micro_fw : float
        Microturbulent velocity for FASTWIND. Default is 10.
            
    heion : int
        Helium ionization stage. Default is 1.
    
    Returns
    -------
    Nothing but the INDAT.DAT files is saved in the specified directory.
    '''

    # Check if necessary directories exist, if not create them:
    if not os.path.exists(dir):
        print('The directory %s does not exist.' % dir)
        return
    elif not os.path.exists(dir+'INDAT/'):
        os.makedirs(dir+'INDAT/')
        print('The directory %s was created.' % dir+'INDAT/')

    # Check that none of the key parameters invalid:
    for var in ['teff', 'lgf', 'logq', 'beta', 'he', 'z', 'micro']:
        if eval(var) == False or eval(var) == None:
            print('Error: %s has not a valid input, please correct.' % var)
            return

    # Calculate logg from lgf:
    logg = lgf + 4 * np.log10(teff / 10000)

    # Define model name (IT CANNOT HAVE MORE THAN 30 CHARACTERS)
    # NOTE: IF YOU MODIFY THE FORMAT HERE, CHANGE IT ALSO IN input_condor.py
    model_name = 't' + '{:05}'.format(int(round(teff))) + \
                 'g' + '{:03}'.format(int(round(logg*100))) + \
                 'q' + '{:03}'.format(int(round(logq*(-10)))) + \
                 'b' + '{:03}'.format(int(round(beta*100))) + \
                 'x' + '{:02}'.format(int(round(micro))) + \
                'he' + '{:03}'.format(int(round(he*100)))

    for extra_elem,extra_elem_name in zip([si, c, n, o, mg], ['si','c','n','o','mg']):
        if extra_elem != False and len(model_name) <= 27:
            model_name = model_name + extra_elem_name + '{:03}'.format(int(round(extra_elem*100)))

    Mdot, R_rsun, v_inf = get_Mdot_R_vinf_Miguel(teff, lgf, logg, logq)

    #Create INDAT.DAT file
    f = open(dir + 'INDAT/' + model_name + '_indat.dat', "w")
    f.write('\'' + model_name + '\'' + '\n') # name of the model
    f.write(' T T           0         100' + '\n') # <True> <True> <Calculate model from 0> <max. number of iterations for convergence>
    f.write('  0.000000000000000E+000' + '\n')  # Internal FASTWIND parameter
    f.write('  ' + str(int(round(teff))) + '        ' + "{:,.2f}".format(logg) + '        ' + "{:,.2f}".format(R_rsun) + '\n')  # teff + logg + R (at tau=2/3)
    f.write('   120.000000000000       0.600000000000000' + '\n')  # <max stellar radius where the atmosphere is calculated>, <Times teff as initial value>
    f.write("{:.6E}".format(Mdot) + '    0.1000 ' + str(int(round(v_inf))) + '    ' + "{:,.2f}".format(beta) + '    0.1000' + '\n')  # <Mdot> <vmin> <v_inf> <beta>
    f.write('  ' + "{:,.2f}".format(he) + '        ' + str(int(round(heion))) + '\n')  # <Yhe> number of electrons that each He atom contributes 
    f.write(' F T F T T' + '\n')  # Internal FASTWIND parameters
    f.write('   ' + str(int(round(micro_fw))) + '        ' + "{:,.2f}".format(z) + '      T T'+ '\n') # <micro> <Z> <True> <True>
    f.write(' T F           1           2' + '\n')  # <True> <False> <teff corrections>
    f.write('1.0000      .10000      .20000' + '\n')  # Input parameters for the clumping law
    f.write('SI ' + "{:,.2f}".format(si) + '\n')  # Si abundance
    if mg != False:
        f.write('MG  ' + "{:,.2f}".format(mg) + '\n')  # Mg abundance
    if c != False:
        f.write('C  ' + "{:,.2f}".format(c) + '\n')  # C abundance
    if n != False:
        f.write('N  ' + "{:,.2f}".format(n) + '\n')  # N abundance
    if o != False:
        f.write('O  ' + "{:,.2f}".format(o) + '\n')  # O abundance

    f.close()

    return model_name


def write_formal(dir, model_name, micro, escattering=1):
    '''
    Write the FORMAL.DAT files for FASTWIND from a set of input parameters.

    Parameters
    ----------
    dir : str
        Path to the directory where the FORMAL.DAT file will be saved.

    name : str
        Name of the model.

    micro : float
        Microturbulent velocity.

    escattering : int
        Flag to include electron scattering. Default is 1.

    Returns
    -------
    Nothing but the FORMAL.DAT files is saved in the specified directory.
    '''

    # Check if necessary directories exist, if not create them:
    if not os.path.exists(dir):
        print('The directory %s does not exist.' % dir)
        return
    elif not os.path.exists(dir+'FORMAL/'):
        os.makedirs(dir+'FORMAL/')
        print('The directory %s was created.' % dir+'FORMAL/')

    f = open(dir + 'FORMAL/' + model_name + '_formal.dat', "w")
    f.write(model_name + '\n') # name
    f.write("{:,.2f}".format(micro) + '\n') # micro
    f.write('{:01}'.format(escattering))  # micro
    f.close()

