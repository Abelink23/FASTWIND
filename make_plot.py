# Core packages
import os
import numpy as np
import re

# Plotting packages
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
plt.rc('xtick', direction='in', top='on')
plt.rc('ytick', direction='in', right='on')


dict_param = {'t':'Teff [K]', 'g':'logg', 'q':'logQ', 'b':'beta', 'x':r'\xi', 'lgf':'log(gf)',
              'he':r'$Y_{\rm He}$', 'Si':r'$\epsilon_{\rm Si}$', 'C':r'$\epsilon_{\rm C}$',
                'N':r'$\epsilon_{\rm N}$', 'O':r'$\epsilon_{\rm O}$'}

def sHR(path_list, y, x):
    """
    Function that makes an sHR diagram for a list of models.
    
    Parameters
    ----------
    path_list : str
        Path to the txt file containing the list of models.
        If you want to separate the good and bad models, you can use the
        extension '_bad' at the end of the model name.
    
    Returns
    -------
    Nothing, but it saves a png file with the plot.
    
    """
    
    with open(path_list, 'r') as f:
        files = f.read().splitlines()
        g_model = [f for f in files if not f.strip().endswith('_bad') and not f == '']
        b_model = [f.split('_bad')[0] for f in files if f.strip().endswith('_bad') and not f == '']
    
    # use the first file to find the parameters (the strings) between the values and make a dictionary of lists
    strings = [re.findall('[a-zA-Z]+', g_model[0])[i] for i in range(len(re.findall('[a-zA-Z]+', g_model[0])))]
    g_params = {s:[] for s in strings}
    b_params = {s:[] for s in strings}
    
    # loop through all the files and append the values to the dictionary
    for f in g_model:
        for s in strings:
            g_params[s].append(int(re.findall(s+'(\d+)', f)[0]))
    
    for f in b_model:
        for s in strings:
            b_params[s].append(int(re.findall(s+'(\d+)', f)[0]))
    
    # make the plot
    fig, ax = plt.subplots(figsize=(8,6))
    #ax.scatter(np.asarray(params['t']),np.asarray(params['t'])**4/np.asarray(params['g'])/100, s=1)
    # loggf = logg - 4 * np.log10(teff * 1e-4)
    g_params['lgf'] = [g_params['g'][i]/100 - 4 * np.log10(g_params['t'][i]/1e4) for i in range(len(g_params['t']))]
    b_params['lgf'] = [b_params['g'][i]/100 - 4 * np.log10(b_params['t'][i]/1e4) for i in range(len(b_params['t']))]
    ax.scatter(np.asarray(g_params[x]),np.asarray(g_params[y]), s=5, c='b', ec='None', alpha=0.6, label='good')
    ax.scatter(np.asarray(b_params[x]),np.asarray(b_params[y]), s=5, c='r', ec='None', label='bad')
    ax.set_xlabel(dict_param[x])
    ax.set_ylabel(dict_param[y])
    
    if x == 't':
        ax.invert_xaxis()
    if y == 'lgf':
        ax.invert_yaxis()
    
    # save the plot
    plt.savefig(path_list.replace('.txt', '.png')+'_%s_%s.png' % (y, x), dpi=300, bbox_inches='tight')


