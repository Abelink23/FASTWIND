# Core packages
import os
import numpy as np

from init import *
from LHS import read_input_grid_lhs
from get_wind_par import *

def input_condor_lhs(grid, lhs_grid_name, model_atom, lines, formal, prescription, sgs, e_scatter=1):
    '''
    Create the input files for Condor from an input Latin Hypercube Sampling grid.
    
    Parameters
    ----------
    grid : str
        Name of the txt file containing the input information of the grid.
        It must be located inside LHS_grids/ folder.
        An example of the format of the file can be found in LHS_grids/test.txt.
        NOTE: IT MUST BE THE SAME FILE USED TO CREATE THE LHS GRID.

    lhs_grid_name : str
        Name of the grid file inside LHS_grids/ folder.
    
    model_atom : str
        Name of the atomic data file used by FASTWIND.
    
    lines : str
        Name of the line data file used by FASTWIND.
    
    formal : str
        Name of the formal data file used by FASTWIND.
    
    prescription and sgs : str and bool
        See get_Mdot_R_vinf() for more information.
    
    e_scatter : int [0/1], optional
        Enable or disable the electron scattering option. 
        Default is enabled (e_scatter=1).
    
    Returns
    -------
    None
    '''

    # Find the grid file inside LHS_grids/ folder:
    print('\nTrying to find the grid file inside LHS_grids/ folder...')
    if not grid.endswith('.txt'):
        grid = grid + '.txt'
    if not os.path.exists(main_dir + 'LHS_grids/' + grid):
        print('Error: grid file not found.\n')
        return
    else:
        grid = read_input_grid_lhs(grid)
        print('Grid file found and read.\n')

    # Find the grid file inside LHS_grids/ folder:
    print('\nTrying to find the LHS file inside LHS_grids/ folder...')
    if not lhs_grid_name.endswith('.npy'):
        lhs_grid_name = lhs_grid_name + '.npy'
    if not os.path.exists(main_dir + 'LHS_grids/' + lhs_grid_name):
        print('Error: LHS file not found.\n')
        return
    else:
        print('Grid file found.\n')

    slhs = np.load(main_dir + 'LHS_grids/' + lhs_grid_name)

    # Create the grid folder inside INPUT/ folder:
    print('\nSearching for the grid folder inside INPUT/ folder...')
    grid_name = lhs_grid_name.replace('.npy', '')
    if not os.path.exists(main_dir + 'INPUT/' + grid_name):
        print('INPUT/ folder not found, creating it...')
        os.makedirs(main_dir + 'INPUT/' + grid_name)
    if not os.path.exists(main_dir + 'INPUT/' + grid_name + '/CONDOR_CHAINFIL/'):
        os.makedirs(main_dir + 'INPUT/' + grid_name + '/CONDOR_CHAINFIL/')

    print('Condor files will be saved inside INPUT/%s/CONDOR_CHAINFIL/ folder.\n' % grid_name)


    # Find names of variables in the grid with a range:
    name_ranges = [key for key in grid.keys() if isinstance(grid[key], tuple) and len(grid[key]) == 2]

    all_abundances = ['si','c','n','o','mg']
    idx_abundances = [i for i in range(len(name_ranges)) if name_ranges[i] in all_abundances]
    abundances = [name_ranges[i] for i in idx_abundances]

    # Assign the arrays of the input variables to the grid dictionary:
    for param in grid.keys():
        if param in name_ranges:
            # Assign the value of the scaled lhs for the parameter:
            tmp_param = slhs[name_ranges.index(param)]
        else:
            # Create an array with the same length as the scaled lhs for the parameter:
            tmp_param = np.full(len(slhs[0]), grid[param])

        # Assign the array to the parameter:
        grid[param] = tmp_param

    # Calculate logg from lgf and teff:
    if not 'logg' in grid.keys():
        grid['logg'] = grid['lgf'] + 4 * np.log10(grid['teff'] / 10000)

    # Create the CHAINFIL files:
    model_names = []; bad_models = 0
    for i in range(len(slhs[0])):

        Mdot, R_rsun, v_inf, v_esc = get_Mdot_R_vinf(teff=grid['teff'][i],
                                                logg=grid['logg'][i],
                                                logq=grid['logq'][i],
                                                prescription=prescription,
                                                sgs=sgs)

        if Mdot == None or R_rsun == None or v_inf == None or v_esc == None:
            print('Skipping CHAINFIL...')
            continue

        chainfil = open(main_dir + 'INPUT/%s/CONDOR_CHAINFIL/CHAINFIL_%i' % (grid_name, i), "w")

        # Write the CHAINFIL file:
        chainfil.write('ATOM ' + model_atom + ' thom_new.dat ' + lines + '\n')

        chainfil.write('# COMPILE' + '\n')
        
        # Create the ABUN string:
        abun_string = 'ABUN'
        for j in range(len(abundances)):
            abun_string += ' ' + abundances[j].upper() + ' ' + "{:,.2f}".format(grid[abundances[j]][i]-12)
        chainfil.write(abun_string + '\n')
        
        chainfil.write('NAME' + '\n')

        # Create the MODEL string (THE MODEL NAME CANNOT HAVE MORE THAN 30 CHARACTERS)
        # NOTE: IF YOU MODIFY THE FORMAT HERE, CHANGE IT ALSO IN informal_dat.py
        model_string = 't' + '{:05}'.format(int(round(grid['teff'][i]))) + \
                       'g' + '{:03}'.format(int(round(grid['logg'][i]*100))) + \
                       'q' + '{:03}'.format(int(round(grid['logq'][i]*(-10)))) + \
                       'b' + '{:03}'.format(int(round(grid['beta'][i]*100))) + \
                       'x' + '{:02}'.format(int(round(grid['micro'][i]))) + \
                      'he' + '{:03}'.format(int(round(grid['he'][i]*100)))

        for extra_elem in ['si','c','n','o','mg']:
            if extra_elem != False and len(model_string) <= 32:
                model_string = model_string + extra_elem + '{:03}'.format(int(round(grid[extra_elem][i]*100)))

        # Check that the model does not exist already:
        while model_string in model_names:
            print('Warning: model %s already exists.' % model_string)
            model_string = input('Please, enter a new name for the model: ')

        model_names.append(model_string)

        # FORMAT [IDL Sergio]            A35,2(I3,1x),   I5,F4.2, F5.1,    E8.2.,F5.0,F4.2,F4.2, F5.1,F4.2, F5.2,F4.2,F4.2)
        # FORMAT (what I use)            A35,     A12,   I5,F4.2, F5.1,    E8.2, F5.0,F4.2,F4.2, F5.1,F4.2,           A15)'
        # t15133g204q134b196x26he019si782     000 100 40000 4.20 007.2 1.28E-07 3511. 1.00 0.10 009.9 1.00 01.00 0.10 0.20
        # model_name                                   teff grav     R 10.^mdot  vinf beta  yhe micro metal   ??   ??   ??
        model_string = model_string + ' 000 100 ' + \
                                      '{:05} '.format(int(round(grid['teff'][i]))) + \
                                      '{:4.2f} '.format(grid['logg'][i]) + \
                                      '{:05.1f} '.format(R_rsun) + \
                                      '{:1.2E} '.format(Mdot) + \
                                      '{:4.0f}. '.format(v_inf) + \
                                      '{:4.2f} '.format(grid['beta'][i]) + \
                                      '{:4.2f} '.format(grid['he'][i]) + \
                                      '{:05.1f} '.format(grid['micro_fw'][i]) + \
                                      '{:4.2f} '.format(grid['z'][i]) + \
                                      '01.00 0.10 0.20\n'

        chainfil.write(model_string)

        chainfil.write('FORMAL %i %s %i' % (e_scatter, formal, int(grid['micro'][i])) + '\n')
        chainfil.write('CLEAN' + '\n' + 'COMPRESS' + '\n' + 'END' + '\n' + '# -----------------------')

        chainfil.close()

    print('Condor files created.')

    return