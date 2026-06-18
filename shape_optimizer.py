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
                    # If it's a list of lists, convert inner elements to tuples to allow sweeping
                    if len(opt_val) > 0 and isinstance(opt_val[0], list):
                        config[key][opt_key] = [tuple(p) for p in opt_val]
                    else:
                        config[key][opt_key] = tuple(opt_val)
    
    # get appropriate kwargs for the optimization
    match config['optimization_type']:
        case 'SA':
            kwargs_optimization = config['kwargs_SA']
        case 'DSA':
            kwargs_optimization = config['kwargs_DSA']
        case _:
            raise ValueError("Invalid optimization type")

    # Get combinations list
    n_runs = config['n_runs']
    match config['combination_type']:
        case 'grid':
            combinations = build_grid_combinations(config['geometric_params'], config['penalizations'], kwargs_optimization)
            combinations = [copy.deepcopy(c) for _ in range(n_runs) for c in combinations]
        case 'zip':
            combinations = build_zip_combinations(config['geometric_params'], config['penalizations'], kwargs_optimization)
            combinations = [copy.deepcopy(c) for _ in range(n_runs) for c in combinations]
        case 'lhs':
            n_samples = config.get('n_lhs_samples', 10)
            combinations = build_lhs_combinations(config['geometric_params'], config['penalizations'], kwargs_optimization, num_samples=n_samples)
            combinations = [copy.deepcopy(c) for _ in range(n_runs) for c in combinations]
        case _:
            raise ValueError("Invalid combination type")
    
    db_path = config['database_name']
    table_name = config['table_name']

    # Dynamic table naming
    if isinstance(table_name, (bool, type(None))):
        x_max, y_max = config['geometric_params']['geometric_config']['x_max'], config['geometric_params']['geometric_config']['y_max']
        n_ellipses = config['geometric_params']['num_ellipses']
        linear_penalization = config['penalizations']['linear']
        table_name = f'x{x_max}y{y_max}n{n_ellipses}lambda{linear_penalization}_{config["optimization_type"]}'

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