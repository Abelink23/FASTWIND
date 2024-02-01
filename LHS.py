# Core packages
import os
import numpy as np

# Plotting packages
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
plt.rc('xtick', direction='in', top='on')
plt.rc('ytick', direction='in', right='on')

# Latin Hypercube Sampling
import lhsmdu

#import warnings
#warnings.filterwarnings('ignore', category=PendingDeprecationWarning)

from init import *
from informal_dat import *

def read_input_grid_lhs(file):
    """
    Read the input grid file and returns the grid input parameters.

    Parameters
    ----------
    file : str
        Path to the input grid file inside LHS/ folder.

    Returns
    -------
    Dictionary with the grid parameters.
    """

    names_all_parameters = ['teff','lgf','logg','logq','beta','he','si','mg','c','n','o','z','micro','micro_fw','heion']

    grid = {}
    with open('LHS_grids/'+file, 'r') as file:
        for line in file:

            if line.startswith('#') or line.startswith('\n'):
                continue

            key, value = line.strip().split('=')
            key = key.strip().lower()
            value = value.split(',')

            if key in grid.keys():
                print('Warning: %s is already defined in the input grid file, it will be overwritten.' % key)

            if key == 'grid_name':
                grid['grid_name'] = value[0].strip()
            else:
                if len(value) == 1:
                    grid[key] = float(value[0])
                elif len(value) == 2:
                    grid[key] = (float(value[0]), float(value[1]))

    # Check if grid_name is defined:
    if 'grid_name' not in grid.keys():
        print('Error: grid_name not defined in the input grid file.')
        return

    # Check that not both lgf and logg are defined:
    if 'lgf' in grid.keys() and 'logg' in grid.keys():
        print('Error: both lgf and logg are defined in the input grid file, please correct.')
        return

    # Check if all the fundamental variables are defined:
    for var in names_all_parameters:
        if (var in ['teff','logq','beta','he','z','micro'] and not var in grid.keys()) or \
            (var in ['lgf','logg'] and not ('lgf' in grid.keys() or 'logg' in grid.keys())):
            print('Error: %s not defined in the input grid file, please correct.' % var)
            return
        elif var not in grid.keys() and var not in ['lgf','logg']:
            if var == 'micro_fw':
                print('Warning: micro_fw not defined in the input grid file.')
                print('Typical values are 10 for SGs and 5 for Gs.')
                grid['micro_fw'] = int(input('Enter an integer value for micro_fw: '))
            elif var == 'heion':
                print('Warning: heion not defined in the input grid file, it will be set to 1.')
                grid['heion'] = 1
            else:
                print('Warning: %s not defined in the input grid file, it will be set to False.' % var)
                grid[var] = False

    return grid


def run_hypercube(grid, n_models, prescription, sgs, make_input_files=True, show_kiel=False):
    """
    Creates  for a given grid of input parameters.

    Parameters
    ----------
    grid : str
        Name of the txt file containing the input information of the grid.
        It must be located inside LHS_grids/ folder.
        An example of the format of the file can be found in LHS_grids/test.txt.

    n_models : int
        The number of models to calculate in the hypercube.

    prescription and sgs : str and bool
        See get_Mdot_R_vinf() for more information.

    make_input_files : bool, optional
        Flag indicating whether to create input files for each model. 
        Default is True.

    show_kiel : bool, optional
        Flag indicating whether to display a Kiel diagram. 
        Default is False.

    Returns:
        None
    """
    
    print('Reading grid file...')
    try:
        grid = read_input_grid_lhs(grid)
    except:
        print('Error reading grid file, please check the format.')
        return

    grid_name = grid['grid_name']

    print('\nDefining a seed for the random number generator...\n')
    lhsmdu.setRandomSeed(123)

    # Find number of variables in the grid with a range:
    n_ranges = sum(1 for key in grid.keys() if isinstance(grid[key], tuple) and len(grid[key]) == 2)
    name_ranges = [key for key in grid.keys() if isinstance(grid[key], tuple) and len(grid[key]) == 2]
    print('\nNumber of variables in the grid that cover a range of values: %i\n' % n_ranges)

    # Calculate hypercube of n_models
    print('\nCalculating hypercube of %i models...' % n_models)
    lhs = lhsmdu.sample(n_ranges, n_models)
    print('Hypercube calculated.\n')

    # Transform to arrays:
    lhs_array = np.array(lhs)

    # Scale to the variable ranges:
    print('\nScaling hypercube to the variable ranges...')
    for i in range(n_ranges):
        lhs_array[i] = lhs_array[i] * (grid[name_ranges[i]][1] - grid[name_ranges[i]][0]) + grid[name_ranges[i]][0]

    # Replace variable lhs_array with slhs (scaled lhs):
    slhs = lhs_array
    print('Hypercube has been scaled.\n')
    del lhs_array
    
    if show_kiel == True:
        plt.title('Kiel diagram')
        plt.scatter(slhs[0], slhs[1])
        plt.show(block=False)

    # Save scaled lhs:
    print('\nSaving scaled hypercube...')
    # Check if there is a grid with the same name:
    if os.path.exists(main_dir + 'LHS_grids/%s.npy' % grid_name):
        print('Warning: a grid with the same name already exists, it will be overwritten.')
        overwrite = input('Do you want to overwrite? (y/n): ')
        if overwrite.lower() == 'y':
            np.save(main_dir + 'LHS_grids/%s.npy' % grid_name, slhs)
            print('Scaled hypercube overwritten.')
        else:
            print('Grid not saved.')
            pass
    else:
        np.save(main_dir + 'LHS_grids/%s.npy' % grid_name, slhs)
        print('Scaled hypercube saved.\n')

    if make_input_files == True:
 
        # Make subfolder inside INPUT/ for the grid:
        if not os.path.exists(main_dir + 'INPUT/' + grid_name):
            print('Creating subfolder inside %sINPUT/ for the grid...' % main_dir)
            os.makedirs(main_dir + 'INPUT/' + grid_name)
            print('Subfolder created.')

        # Assign the arrays of the input variables for write_indat and write_formal functions:
        for param in grid.keys():
            if param in name_ranges:
                # Assign the value of the scaled lhs for the parameter:
                tmp_param = slhs[name_ranges.index(param)]
            else:
                # Create an array with the same length as the scaled lhs for the parameter:
                tmp_param = np.full(len(slhs[0]), grid[param])
            
            # Assign the array to the parameter:
            grid[param] = tmp_param

        # Make the trick in case logg is defined instead of lgf (in informal_dat.py is de-converted):
        if 'logg' in grid.keys():
            grid['lgf'] = grid['logg'] - 4 * np.log10(grid['teff'] / 10000)
            del grid['logg']

        f = open(main_dir + 'INPUT/%s/ModelList_%s.txt' % (grid_name, grid_name), "w")
        fj = open(main_dir + 'INPUT/%s/jobs.list_%s_1' % (grid_name, grid_name), "w")

        bad_models = 0
        for i in range(n_models):
            # Create the INDAT.dat file:
            model_name = write_indat(main_dir + 'INPUT/%s/' % grid_name,
                                     teff=grid['teff'][i],
                                     lgf=grid['lgf'][i],
                                     logq=grid['logq'][i],
                                     beta=grid['beta'][i],
                                     z=grid['z'][i],
                                     he=grid['he'][i],
                                     si=grid['si'][i],
                                     mg=grid['mg'][i],
                                     c=grid['c'][i],
                                     n=grid['n'][i],
                                     o=grid['o'][i],
                                     micro=grid['micro'][i],
                                     micro_fw=grid['micro_fw'][i],
                                     heion=grid['heion'][i],
                                     prescription=prescription,
                                     sgs=sgs)

            # If the model is not valid, skip the FORMAL_INPUT.dat file:
            if model_name is None:
                print('Skipping FORMAL_INPUT for model %s' % model_name)
                bad_models += 1
                continue

            # Create the FORMAL_INPUT.dat file:
            write_formal(main_dir + 'INPUT/%s/' % grid_name,
                         model_name,
                         grid['micro'][i],
                         scattering=1)

            # Create list of model names:
            f.write(model_name + '\n')
            fj.write(model_name + '\n')

            k = i + 1
            if k % 100 == 0:
                fj.write('END')
                fj.close()
                if k != n_models:
                    fj = open(main_dir + 'INPUT/%s/jobs.list_%s_%s' % (grid_name, grid_name, int(k/100+1)), "w")

        if k != n_models:
            fj.write('END')
            fj.close()
        
        f.close()
        
        print('Input files created.')
        print('Number of bad models: %i/%i' % (bad_models, n_models))
        
        return None
