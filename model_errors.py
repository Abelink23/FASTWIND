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

def plot_bad_models(path_list, y, x, path_chainfils=None):
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
        g_models = [f for f in files if not f.strip().endswith('_bad') and not f == '']
        b_models = [f.split('_bad')[0] for f in files if f.strip().endswith('_bad') and not f == '']
    
    # use the first file to find the parameters (the strings) between the values and make a dictionary of lists
    strings = [re.findall('[a-zA-Z]+', g_models[0])[i] for i in range(len(re.findall('[a-zA-Z]+', g_models[0])))]
    g_params = {s:[] for s in strings}
    b_params = {s:[] for s in strings}
    
    # loop through all the files and append the values to the dictionary
    for m in g_models:
        for s in strings:
            g_params[s].append(int(re.findall(s+'(\d+)', m)[0]))
    
    for m in b_models:
        for s in strings:
            b_params[s].append(int(re.findall(s+'(\d+)', m)[0]))
    
    # load the CHAINFIL files and make a dictionary with the model name
    if path_chainfils is not None:
        dict_chainfils = {}
        for c in os.listdir(path_chainfils):
            if c.startswith('CHAINFIL_'):
                with open(path_chainfils + c, 'r') as f:
                    content = f.read().splitlines()
                    model_name = content[4].split()[0]
                    dict_chainfils[model_name] = c
    
    # make the plot
    fig, ax = plt.subplots(figsize=(8,6))
    # loggf = logg - 4 * np.log10(teff * 1e-4)
    g_params['lgf'] = [g_params['g'][i]/100 - 4 * np.log10(g_params['t'][i]/1e4) for i in range(len(g_params['t']))]
    b_params['lgf'] = [b_params['g'][i]/100 - 4 * np.log10(b_params['t'][i]/1e4) for i in range(len(b_params['t']))]
    ax.scatter(np.asarray(g_params[x]),np.asarray(g_params[y]), s=5, c='b', ec='None', alpha=0.4, label='good')
    ax.scatter(np.asarray(b_params[x]),np.asarray(b_params[y]), s=15, c='r', ec='None', label='bad')
    ax.set_xlabel(dict_param[x])
    ax.set_ylabel(dict_param[y])
    
    # add the CHAINFIL number if the folder of the CHAINFIL files is given
    if path_chainfils is not None:
        # put the CHAINFIL number in the plot only for the bad models
        for m in range(len(b_models)):
            ax.text(b_params[x][m], b_params[y][m], dict_chainfils[b_models[m]].split('_')[1], fontsize=6)
    
    if x == 't':
        ax.invert_xaxis()
    if y == 'lgf':
        ax.invert_yaxis()
    
    # save the plot
    plt.savefig(path_list.replace('.txt','')+'_%s_%s.png' % (y, x), dpi=300, bbox_inches='tight')


# function to find the chainfil for a given model that has failed
# the name of the models are listed in a txt file (path_list)
# the chainfils are in a directory (path_chainfils)
# the name of the model is contained in the content of the chainfil
# the program has to make a dictionaty with CHAINFIL_<n>: <model_name> where model name is in the fourth line of the chainfil
def find_chainfil(path_list, path_chainfils):
    """
    Function that finds the CHAINFIL names for a list of models that have failed.
    
    Parameters
    ----------
    path_list : str
        Path to the txt file containing the list of models that have failed (no extension).
    
    path_chainfils : str
        Path to the directory containing the CHAINFIL files.
    
    Returns
    -------
    Nothing, but it saves a txt file with the list of CHAINFIL files.
    
    """
    
    with open(path_list, 'r') as f:
        models = f.read().splitlines()
        models = [m.strip() for m in models if not m == '']
    
    dict_chainfils = {}
    for c in os.listdir(path_chainfils):
        if c.startswith('CHAINFIL_'):
            with open(path_chainfils + c, 'r') as f:
                content = f.read().splitlines()
                model_name = content[4].split()[0]
                dict_chainfils[c] = model_name
    
    with open(path_list.split('.')[0] + '_chainfils.txt', 'w') as f:
        for m in models:
            for k, v in dict_chainfils.items():
                if m == v:
                    f.write(m + ' -> ' + k + '\n')
    
    return print('Done!')