
from LHS import *

grid = read_input_grid_lhs('test.txt')
grid_name=grid['grid_name'].strip()

# Test LHS.py:
run_hypercube('BSGs_HHeSiCNO.txt', 100, make_input_files=True, show_kiel=False)

from input_condor import *
# Test input_condor.py:
input_condor_lhs(grid='BSGs_HHeSiCNO.txt', lhs_grid_name='BSGs_HHeSiCNO', file_atom='A10HHeCnewNOSi_MAUI.dat', file_lines='LINES_CNOSi_new_coll_MAUI.dat')


grid='BSGs_HHeSiCNO.txt'; lhs_grid_name='BSGs_HHeSiCNO.npy'
slhs = np.load(main_dir + 'LHS_grids/' + lhs_grid_name)
'teff',slhs[0][0]
'logg',slhs[1][0] + 4 * np.log10(slhs[0][0] / 10000)
'micro',slhs[2][0]
'beta',slhs[3][0]
'logq',slhs[4][0]
'he',slhs[5][0]
'c',slhs[6][0]
'n',slhs[7][0]
'o',slhs[8][0]
'si',slhs[9][0]

from get_wind_par import *
Mdot, R, v_inf = get_Mdot_R_vinf_Miguel(teff=slhs[0][0], lgf=slhs[1][0], logg=slhs[1][0] + 4 * np.log10(slhs[0][0] / 10000), logq=slhs[4][0])

'teff'=15133
'logg'=2.04
'micro'=26
'beta'=1.96
'logq'=-13.38
'he'=0.188
'c'=7.35
'n'=7.55
'o'=8.36
'si'=7.82

R = 9.2984519
print('{:05.1f}'.format(R))

Mdot = 3.56611091791121e-08
print('{:05.1E}'.format(Mdot))

v_inf = 540.81979
print('{:05.1f}'.format(v_inf))