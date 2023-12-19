# Core packages
import os

def read_paths(file_path):
    paths = {}
    with open(file_path, 'r') as file:
        for line in file:
            key, value = line.strip().split('=')
            paths[key.strip()] = value.strip()
            if not paths[key.strip()].endswith('/'):
                paths[key.strip()] = paths[key.strip()] + '/'
    return paths

paths = read_paths('paths.txt')

main_dir = paths.get('main')
models_dir = paths.get('models')
input_dir = paths.get('input')

if not os.path.exists(main_dir+'LHS_grids/'):
    os.makedirs(main_dir+'LHS_grids/')
    
if not os.path.exists(main_dir+'INPUT/'):
    os.makedirs(main_dir+'INPUT/')