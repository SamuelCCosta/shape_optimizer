from annealing import *
from square_solver import *
import sys
from orchestrator import *
import yaml
import copy

def shape_optimizer(config_file : str):
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    
    # Make the key perturbations tuples
    # (in optimization kwargs, with "perturbation" as part of the key)
    for key, val in config.items():
        if "kwargs_" in key:
            for opt_key, opt_val in val.items():
                if "perturbation" in opt_key and isinstance(opt_val, list):
                        config[key][opt_key] = tuple(opt_val)
    
    # Get combinations list
    n_runs = config['n_runs']
    match config['combination_type']:
        case 'grid':
            combinations = build_grid_combinations(config['geometric_params'], config['penalizations'], config['kwargs_optimization'])
            combinations = [copy.deepcopy(c) for _ in range(n_runs) for c in combinations]
        case 'zip':
            combinations = build_zip_combinations(config['geometric_params'], config['penalizations'], config['kwargs_optimization'])
            combinations = [copy.deepcopy(c) for _ in range(n_runs) for c in combinations]
        case _:
            raise ValueError("Invalid combination type")
        
    db_path = config['database_name']
    table_name = config['table_name']
    max_processes = config['max_processes']
    
    # Run experiments
    match config['optimization_type']:
        case 'SA':
            run_experiments(optimization_worker_SA, combinations, db_path=db_path, table_name=table_name, max_processes=max_processes)
        case 'DSA':
            run_experiments(optimization_worker_DSA, combinations, db_path=db_path, table_name=table_name, max_processes=max_processes)
        case _:
            raise ValueError("Invalid optimization type")


if __name__ == '__main__':
    config_file = sys.argv[1]
    shape_optimizer(config_file)