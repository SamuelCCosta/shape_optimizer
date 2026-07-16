from annealing import *
from square_solver import *
import numpy as np
import random
import multiprocessing as mp
from multiprocessing.pool import Pool
from multiprocessing.connection import wait
from pebble import concurrent
import time
import copy
import itertools
import os
import sqlite3
from db_utils import database_writer, get_column_names_dict, unfold_parameters, get_cpu_model


class EllipseSA(BaseSimulatedAnnealing):
    '''Classical SA algorithm, needs an initial suitable configuration to progress with the optimization'''
    def __init__(self, geometric_params : dict, penalizations : dict, **kwargs):
        self.gen_seed = kwargs.pop('gen_seed', 0)
        self.generation_attempts = 0 # Guarantee the seed works or gives a deterministic suitable seed

        self.num_ellipses = geometric_params['num_ellipses']
        perturbation = list(kwargs.pop('perturbation'))
        self.perturbation = perturbation * self.num_ellipses
        
        # stopping criteria
        self.no_improvement = 0
        self.improvement_threshold = kwargs.pop('improvement_threshold', 50)

        super().__init__(**kwargs)
        self.geometric_params = geometric_params
        self.initial_state = []

        #penalization related
        self.full_area = geometric_params['geometric_config']['x_max'] * geometric_params['geometric_config']['y_max']
        self.linear_penalization = penalizations['linear']

        #keyword parameters for solver / bundle constructors
        self.sqs_params = {k: v for k, v in geometric_params.items() if k not in ("num_ellipses")}
        self.ellipse_bundle_params = {k: geometric_params[k] for k in ["geometric_config", "h", "num_ellipses"]}

    def raw_cost_function(self, state) -> tuple[float, bool]:
        '''Get full cost, including (linear) penalization.'''
        ellipses = EllipseBundle(**self.ellipse_bundle_params)
        n_param = 5 #ellipse parameter number
        try:
            for i in range(self.num_ellipses):
                idx = i * n_param
                ellipses.add(Ellipse(state[idx], state[idx + 1], state[idx + 2],
                                    state[idx + 3], state[idx + 4]))
        except Exception: # fails at conception (invalid configuration)
            return float('inf'), False
            
        try:
            sqs = SquareSolver(**self.sqs_params)
            area_percent = (self.full_area - ellipses.area()) / self.full_area
            return sqs.solve(ellipses) + self.linear_penalization * area_percent, True
        except Exception: # fails in solver
            return float('inf'), True
    
    def get_neighbour(self, state):
        noise = [random.normalvariate(0.0, 1.0) * p for p in self.perturbation]
        #noise = [random.uniform(-p, p) for p in self.perturbation] 
        return [s + n for s, n in zip(state, noise)]
    
    def raw_generation_function(self, seed=0):
        '''Generates a random valid configuration, evaluates it, and returns cost & state.'''
        ellipses = EllipseBundle(**self.ellipse_bundle_params)
        
        ellipses.generate_random(seed=seed)
        
        if len(ellipses.bundle) < self.num_ellipses:
            return None, float('inf')
            
        state = []
        for e in ellipses.bundle:
            state.extend([float(e.center[0]), float(e.center[1]), float(e.quadratic_form[0,0]), float(e.quadratic_form[0,1]), float(e.quadratic_form[1,1])])
            
        try:
            sqs = SquareSolver(**self.sqs_params)
            area_percent = (self.full_area - ellipses.area()) / self.full_area
            cost = sqs.solve(ellipses) + self.linear_penalization * area_percent #OBJECTIVE FUNCTION
            return state, cost
        except Exception:
            return None, float('inf')
    
    def get_initial_state(self):
        while self.initial_state == []:
            try:
                current_seed = self.gen_seed + self.generation_attempts if self.gen_seed != 0 else 0
                self.generation_attempts += 1
                with ProcessPool(max_workers=1, context=spawn_ctx) as pool:
                    future = pool.schedule(self.raw_generation_function, args=[current_seed], timeout=60.0)
                    state, cost = future.result() # type: ignore
                    if cost != float('inf'):
                        self.initial_state = state
                        self.initial_cost = cost
            except Exception:
                pass # Ignore timeouts and crashes, try again
    
    def stopping_criteria(self):
        '''Continue if temp is above min_temp. If temp is low, also check for recent improvements.'''
        if self.temp <= self.min_temp:
            return False # Stop if temperature is below the absolute minimum
        if self.temp < 10 * self.min_temp:
            return self.no_improvement <= self.improvement_threshold # If temp is low, stop if there are no recent improvements
        return True # Otherwise, continue
    
    def run_best(self):
        self.no_improvement += 1


class EllipseDSA(DirectSimulatedAnnealing):
    '''Direct SA algorithm'''
    def __init__(self, geometric_params : dict, penalizations : dict, **kwargs):
        self.gen_seed = kwargs.pop('gen_seed', 0)
        self.generation_attempts = 0
        super().__init__(**kwargs)
        self.geometric_params = geometric_params
        self.num_ellipses = geometric_params['num_ellipses']
        self.perturbation = self.perturbation * self.num_ellipses
        self.min_perturbation = self.min_perturbation * self.num_ellipses

        #penalization related
        self.full_area = geometric_params['geometric_config']['x_max'] * geometric_params['geometric_config']['y_max']
        self.linear_penalization = penalizations['linear']

        #keyword parameters for solver / bundle constructors
        self.sqs_params = {k: v for k, v in geometric_params.items() if k not in ("num_ellipses")}
        self.ellipse_bundle_params = {k: geometric_params[k] for k in ["geometric_config", "h", "num_ellipses"]}
    
    def raw_cost_function(self, state) -> tuple[float, bool]:
        '''Get full cost, including (linear) penalization.'''
        ellipses = EllipseBundle(**self.ellipse_bundle_params)
        n_param = 5 #ellipse parameter number
        try:
            for i in range(self.num_ellipses):
                idx = i * n_param
                ellipses.add(Ellipse(state[idx], state[idx + 1], state[idx + 2],
                                    state[idx + 3], state[idx + 4]))
        except Exception: # fails at conception (invalid configuration)
            return float('inf'), False
            
        try:
            sqs = SquareSolver(**self.sqs_params)
            area_percent = (self.full_area - ellipses.area()) / self.full_area
            return sqs.solve(ellipses) + self.linear_penalization * area_percent, True #OBJECTIVE FUNCTION
        except Exception: # fails in solver
            return float('inf'), True
    
    def get_neighbour(self, state):
        noise = [random.normalvariate(0.0, 1.0) * scale for scale in self.perturbation]
        return [s + n for s, n in zip(state, noise)]
    
    def raw_generation_function(self, seed=0):
        '''Generates a random valid configuration, evaluates it, and returns cost & state.'''
        ellipses = EllipseBundle(**self.ellipse_bundle_params)
        
        ellipses.generate_random(seed=seed)
        
        if len(ellipses.bundle) < self.num_ellipses:
            return None, float('inf')
            
        state = []
        for e in ellipses.bundle:
            state.extend([float(e.center[0]), float(e.center[1]), float(e.quadratic_form[0,0]), float(e.quadratic_form[0,1]), float(e.quadratic_form[1,1])])
            
        try:
            sqs = SquareSolver(**self.sqs_params)
            area_percent = (self.full_area - ellipses.area()) / self.full_area
            cost = sqs.solve(ellipses) + self.linear_penalization * area_percent #OBJECTIVE FUNCTION
            return state, cost
        except Exception:
            return None, float('inf')

    def get_starting_configs(self):
        self.initial_configs = []
        print(f"Generating {self.num_configs} starting configurations...")
        while len(self.initial_configs) < self.num_configs:
            try:
                current_seed = self.gen_seed + self.generation_attempts if self.gen_seed != 0 else 0
                self.generation_attempts += 1
                with ProcessPool(max_workers=1, context=spawn_ctx) as pool:
                    future = pool.schedule(self.raw_generation_function, args=[current_seed], timeout=60.0)
                    state, cost = future.result() # type: ignore
                    if cost != float('inf') and state is not None:
                        self.initial_configs.append((state, cost))
            except Exception:
                pass # Ignore timeouts and crashes, try again
        self.initial_configs.sort(key=lambda x: x[1])

    def test_generation(self):
        self.get_starting_configs()
        print(len(self.initial_configs), self.initial_configs)


def optimization_worker_SA(queue : mp.Queue, run_id: int, geometric_params : dict, 
                        penalizations : dict, kwargs_SA : dict, initial_params : list = []):
    '''Defines the routine of an optimization worker'''
    solver = EllipseSA(geometric_params, penalizations, **kwargs_SA)
    
    if initial_params == []:
        solver.get_initial_state()
        actual_initial_params = solver.initial_state
        initial_cost = getattr(solver, 'initial_cost', None)
    else:
        actual_initial_params = initial_params
        initial_cost = None
    print(actual_initial_params)
    start_time = time.time()
    best_param, best_cost = solver.run(actual_initial_params, initial_cost=initial_cost)
    runtime = time.time() - start_time

    # Get all the values into the SQL queue
    update_parameters = {
        'run_id': run_id,
        'runtime': runtime,
        'best_param': best_param,
        'best_cost': best_cost,
        'initial_params': actual_initial_params,
        'cpu_model': get_cpu_model()
    }
    #Convert all non-number values into JSON strings
    update_parameters = unfold_parameters(update_parameters)

    queue.put(update_parameters)

def optimization_worker_DSA(queue : mp.Queue, run_id: int, geometric_params : dict, 
                        penalizations : dict, kwargs_DSA : dict):
    '''Defines the routine of a DSA optimization worker'''
    solver = EllipseDSA(geometric_params, penalizations, **kwargs_DSA)
    start_time = time.time()
    best_param, best_cost = solver.run() # DSA returns (state, cost)
    runtime = time.time() - start_time

    # Get all the values into the SQL queue
    update_parameters = {
        'run_id': run_id,
        'runtime': runtime,
        'best_param': best_param,
        'best_cost': best_cost,
        'initial_params': solver.initial_configs, #store all initializers
        'cpu_model': get_cpu_model()
    }
    #Convert all non-number values into JSON strings
    update_parameters = unfold_parameters(update_parameters)

    queue.put(update_parameters)

def _traverse_and_find_lists(config_bundle):
    paths = []
    sweep_lists = []
    
    def traverse(obj, path):
        if isinstance(obj, dict):
            for k, v in obj.items():
                traverse(v, path + [k])
        elif isinstance(obj, list):
            paths.append(path)
            sweep_lists.append(obj)
            
    traverse(config_bundle, [])
    return paths, sweep_lists

def build_grid_combinations(geometric_params, penalizations, kwargs_optimization):
    config_bundle = {
        'geometric_params': geometric_params,
        'penalizations': penalizations,
        'kwargs_optimization': kwargs_optimization
    }
    paths, sweep_lists = _traverse_and_find_lists(config_bundle)
    
    if not sweep_lists:
        return [config_bundle]
        
    combinations = []
    for combo in itertools.product(*sweep_lists):
        new_bundle = copy.deepcopy(config_bundle)
        for path, val in zip(paths, combo):
            target = new_bundle
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = val
        combinations.append(new_bundle)
        
    return combinations

def build_zip_combinations(geometric_params, penalizations, kwargs_optimization):
    config_bundle = {
        'geometric_params': geometric_params,
        'penalizations': penalizations,
        'kwargs_optimization': kwargs_optimization
    }
    paths, sweep_lists = _traverse_and_find_lists(config_bundle)
    
    if not sweep_lists:
        return [config_bundle]
        
    length = len(sweep_lists[0])
    if not all(len(lst) == length for lst in sweep_lists):
        raise ValueError("For zipped combinations, all parameter lists must have the same length.")
        
    combinations = []
    for combo in zip(*sweep_lists):
        new_bundle = copy.deepcopy(config_bundle)
        for path, val in zip(paths, combo):
            target = new_bundle
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = val
        combinations.append(new_bundle)
        
    return combinations

def build_lhs_combinations(geometric_params, penalizations, kwargs_optimization, num_samples):
    config_bundle = {
        'geometric_params': geometric_params,
        'penalizations': penalizations,
        'kwargs_optimization': kwargs_optimization
    }
    paths, sweep_lists = _traverse_and_find_lists(config_bundle)
    
    if not sweep_lists:
        return [config_bundle]
        
    lhs_lists = []
    for lst in sweep_lists:
        # Create a stratified sample of indices covering the list evenly
        indices = [int(i * len(lst) / num_samples) for i in range(num_samples)]
        random.shuffle(indices)
        lhs_lists.append([lst[i] for i in indices])
        
    combinations = []
    for combo in zip(*lhs_lists):
        new_bundle = copy.deepcopy(config_bundle)
        for path, val in zip(paths, combo):
            target = new_bundle
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = val
        combinations.append(new_bundle)
        
    return combinations

def run_experiments(worker_target, combinations, extra_worker_args=(),
                    db_path='experiments.db', table_name='results', max_processes=10):
        
    if not combinations:
        print("No combinations to run.")
        return

    first_bundle = combinations[0]
    column_names = get_column_names_dict(
        geometric_params=first_bundle['geometric_params'],
        penalizations=first_bundle['penalizations'],
        kwargs_optimization=first_bundle['kwargs_optimization'],
        paramtype=list
    )
    
    # --- SETUP DATABASE & INITIAL INSERTS ---
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute(f"CREATE TABLE IF NOT EXISTS {table_name} (run_id INTEGER PRIMARY KEY AUTOINCREMENT)")
    cursor.execute(f"PRAGMA table_info({table_name})")
    existing_columns = {row[1] for row in cursor.fetchall()}
    for name, sql_type in column_names.items():
        if name not in existing_columns:
            if sql_type not in ('REAL', 'INTEGER', 'TEXT'):
                sql_type = 'REAL'
            cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN {name} {sql_type}")
            
    db_name = os.path.splitext(os.path.basename(db_path))[0]
    
    for combo_bundle in combinations:
        run_parameters = combo_bundle['kwargs_optimization'].copy()
        for pen_type in combo_bundle['penalizations']:
            run_parameters[pen_type + '_penalization'] = combo_bundle['penalizations'][pen_type]
        run_parameters |= combo_bundle['geometric_params']
        run_parameters = unfold_parameters(run_parameters)
        
        columns = list(run_parameters.keys())
        columns_sql = ', '.join(columns)
        placeholders_sql = ', '.join([f':{col}' for col in columns])
        insert_sql = f'INSERT INTO {table_name} ({columns_sql}) VALUES ({placeholders_sql})'
        cursor.execute(insert_sql, run_parameters)
        run_id = cursor.lastrowid
        combo_bundle['run_id'] = run_id
        
        track_file_name = combo_bundle['kwargs_optimization'].get('track_file_name')
        if track_file_name is not None and not isinstance(track_file_name, str):
            csv_dir = f'track_csv/{db_name}/{table_name}'
            os.makedirs(csv_dir, exist_ok=True)
            path = f'{csv_dir}/{run_id}.csv'
            combo_bundle['kwargs_optimization']['track_file_name'] = path
            cursor.execute(f"UPDATE {table_name} SET track_file_name = ? WHERE run_id = ?", (path, run_id))
            
    conn.commit()
    conn.close()
    
    spawn_ctx = mp.get_context('spawn')
    queue = spawn_ctx.Queue()
    db_writer = spawn_ctx.Process(target=database_writer, args=(queue, db_path, table_name))
    db_writer.start()

    active_workers = []

    for combo_bundle in combinations:
        if len(active_workers) >= max_processes:
            sentinels = [w.sentinel for w in active_workers]
            wait(sentinels)
            active_workers = [w for w in active_workers if w.is_alive()]

        w = spawn_ctx.Process(target=worker_target, 
                              args=(queue, combo_bundle['run_id'], combo_bundle['geometric_params'], 
                              combo_bundle['penalizations'], combo_bundle['kwargs_optimization'], *extra_worker_args))
        w.start()
        active_workers.append(w)
    
    for w in active_workers:
        w.join()
    
    queue.put(None)
    db_writer.join()
    print('Done')


    '''
    Geometric info (constant in every optimization) : geometric_params
    Penalization info : penalizations
    Any SA/DSA info : best_params, best_cost, runtime
    SA specific info : kwargs_SA, initial_params
    DSA specific info : kwargs_DSA, ... (WIP)
    '''





    
