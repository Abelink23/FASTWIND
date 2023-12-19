import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import special
import AA_LIBRARY_basic as mybasic
import seaborn as sns
from scipy.io.idl import readsav
import LowZCal_Params as mycal
from astropy.convolution import convolve, Box1DKernel
from astropy.stats import sigma_clip
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mycolorpy import colorlist as mcp
import glob
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import lhsmdu
mpl.rcParams['figure.dpi'] = 120
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
import warnings
warnings.filterwarnings('ignore', category=PendingDeprecationWarning)

# --------------------------------------------------- FUNCTIONS ------------------------------------------------------ #


def write_indat(dir, teff, lgf, logQ, beta, z, he, micro, si, mg = False, c= False, n= False, o= False, micro_fw = 10, heion = 1):

    logg = lgf + 4 * np.log10(teff / 10000)

    # Define name:
    name = 't' + '{:05}'.format(int(teff)) + 'g' + '{:03}'.format(int(logg*100)) + 'q' + '{:02}'.format(int(logQ*(-1))) \
           + 'b' + '{:02}'.format(int(beta)) + 'he' + '{:03}'.format(int(he*100)) + 'z' + '{:03}'.format(int(z*100)) + 'x' + '{:02}'.format(int(micro))
    '''
           + '_si' + '{:03}'.format(int(si*100))
    if mg != False:
        name = name + 'mg' + '{:03}'.format(int(mg*100))
    if c != False:
        name = name + 'c' + '{:03}'.format(int(c*100))
    if n != False:
        name = name + 'n' + '{:03}'.format(int(n*100))
    if o != False:
        name = name + 'o' + '{:03}'.format(int(o*100))

    name = name + '_xi' + '{:04}'.format(int(micro*100))
    '''

    # Calculate R and wind parameters:
    Mbol_sun = 4.74  # Cox 2000
    SB_constant = 5.6704e-5  # Stefan-Boltzmann constant [g s-3 K-4]
    Lsun = 3.839e33  # Solar luminosity [erg s-1]
    Rsun = 6.6599E10  # Solar radius [cm]

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
    L = 10 ** logL * Lsun
    R = np.sqrt(L / (4 * np.pi * SB_constant * teff ** 4))
    R_rsun = R / Rsun

    # calculate wind parameters:
    vsc = np.sqrt(2 * 10 ** logg * (R)) * 1.E-5
    t = teff * 1.E-4
    if t <= 1.7:
        vinfty = 1.25 * vsc

    elif t <= 2.3:
        vinfty = (0.2917 * (teff * 10 ** -3) - 3.7083) * vsc
    else:
        vinfty = 3. * vsc

    # Q = Mdot/R*vinf)**1.5
    Mdot = 10 ** logQ * (R_rsun * vinfty) ** 1.5


    #Create INDAT.DAT file
    f = open(dir + name + '_indat.dat', "w")
    f.write('\'' + name + '\'' + '\n') # name
    f.write(' T T           0         100' + '\n') # fw
    f.write('  0.000000000000000E+000' + '\n')  # fw
    f.write('  ' + str(int(teff)) + '        ' + "{:,.2f}".format(logg) + '        ' + "{:,.2f}".format(R_rsun) + '\n')  # teff + logg + R
    f.write('   120.000000000000       0.600000000000000' + '\n')  # fw: Rmax, Rmin
    f.write("{:.6E}".format(Mdot) + '    0.1000 ' + str(int(vinfty)) + '    ' + "{:,.2f}".format(beta) + '    0.1000' + '\n')  # Mdot vinfty beta
    f.write('  ' + "{:,.2f}".format(he) + '        ' + str(int(heion)) + '\n')  # yhe
    f.write(' F T F T T' + '\n')  # fw
    f.write('   ' + str(int(micro_fw)) + '        ' + "{:,.2f}".format(z) + '      T T'+ '\n')
    f.write(' T F           1           2' + '\n')  # fw
    f.write('1.0000      .10000      .20000' + '\n')  # fw
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

    return name

def write_formal(dir, name, micro, escattering = 1):
    # Create FORMAL.DAT file
    f = open(dir + name + '_formal.dat', "w")
    f.write(name + '\n') # name
    f.write("{:,.2f}".format(micro) + '\n') # micro
    f.write('{:01}'.format(escattering))  # micro
    f.close()

# ---------------------------------------------------   PATHS   ------------------------------------------------------ #

project_dir = '/home/mlorenzo/Documents/Work/Thesis/Projects/SexA_Bstars/'
data_dir = project_dir + 'Data/'
grid_dir = project_dir + 'Grid/'
test_dir = '/home/mlorenzo/Documents/Work/Quick/'


# ---------------------------------------------------   RANGES   ----------------------------------------------------- #

grid_name = 'BSG_SEXA_CNOSiMg_Q14'

# Constants used:


# Constants in the grid:
z = 0.1
beta = 1
logQ = -14

# Variables in the grid:
teff_range = (15000, 35000)
loggf_range = (1.00, 1.80) #mlg test
micro_range = (0, 30)
he_range = (0.07, 0.3)
c_range = (6.00, 8.00)
n_range = (6.00,8.00)
o_range = (7.05, 8.25)
mg_range = (6.06, 7.06)
si_range = (6.00, 7.00)

var_list = [teff_range, loggf_range, micro_range, he_range, c_range, n_range, o_range, mg_range, si_range]
var_dict = {'teff':0, 'lgf':1, 'micro':2, 'he':3, 'c':4, 'n':5, 'o':6, 'mg':7, 'si':8}

# ---------------------------------------------------   STEPS   ------------------------------------------------------ #

calc_hypercube = False
write_ALL_files = True
test_ONE_files = False

# ---------------------------------------------------   CALLS   ------------------------------------------------------ #
print('Hello world!')
c = 3*10**5

if calc_hypercube == True:
    # Define seed:
    lhsmdu.setRandomSeed(123)

    # Calculate hypercube of 800 models
    lhs = lhsmdu.sample(len(var_list), 800)

    # Tranform to arrays:
    lhs_array = np.array(lhs)

    # Scale to the variable ranges:
    scaled_lhs_array = lhs_array
    i = 0
    for var in var_list:
        scaled_lhs_array[i] =  scaled_lhs_array[i] * (var[1] - var[0]) + var[0]
        i += 1

    # Plot diagram:
    plt.title('Kiel diagram')
    plt.scatter(scaled_lhs_array[0], scaled_lhs_array[1])
    plt.show()

    # Save scaled lhs:
    np.save(grid_dir + 'LHS_' + grid_name + '.npy', scaled_lhs_array)

if write_ALL_files == True:
    #Load scaled lhs:
    slhs = np.load(grid_dir + 'LHS_' + grid_name + '.npy')

    # Create the indat.dat's and formal.dat's of all models and a list of their names:
    n_models = len(slhs[0])

    f = open(grid_dir + 'ModelList_' + grid_name + '.txt', "w")
    fj = open(grid_dir + 'jobs.list_' + grid_name + '_1', "w")
    j = 2
    for i in range(n_models):
        k = i + 1

        # Cretae indat.dat's:
        model_name = write_indat(grid_dir + 'INDAT.DAT/', teff=slhs[var_dict['teff']][i], lgf=slhs[var_dict['lgf']][i],
                                 logQ = logQ, beta = beta, z=z, he=slhs[var_dict['he']][i],
                                 si=slhs[var_dict['si']][i], mg = slhs[var_dict['mg']][i],
                                 c=slhs[var_dict['c']][i], n=slhs[var_dict['n']][i], o=slhs[var_dict['o']][i],
                                 micro=slhs[var_dict['micro']][i], micro_fw=10)
        # Cretae formal.dat's:
        write_formal(grid_dir + 'FORMAL.DAT/', model_name, slhs[var_dict['micro']][i], escattering=1)

        # Create list of model names:
        f.write( model_name + '\n') # name

        fj.write(model_name + '\n')  # name
        if k % 100 == 0:
            fj.write('END')
            fj.close()
            fj = open(grid_dir + 'jobs.list_' + grid_name + '_' + str(j), "w")
            j+=1


    f.close()
    fj.write('END')
    fj.close()

if test_ONE_files == True:
    model_name = write_indat(test_dir, teff=28718, lgf=1.35, logQ=-12.88, beta=3.52, z=z,
                he=0.07,  si=6.34, c=7.32, n=6.93, o=7.71, micro = 28.8, micro_fw = 10)
    write_formal(test_dir, model_name, 28.8, escattering=1)
