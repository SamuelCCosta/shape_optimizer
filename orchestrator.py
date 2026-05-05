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
from db_utils import database_writer, get_column_names_dict, unfold_parameters


class EllipseSA(BaseSimulatedAnnealing):
    '''Classical SA algorithm, needs an initial suitable configuration to progress with the optimization'''
    def __init__(self, scales, geometric_params : dict, penalizations : dict, **kwargs):
        super().__init__(**kwargs)
        self.geometric_params = geometric_params
        self.num_ellipses = geometric_params['num_ellipses']
        self.scales = scales * self.num_ellipses
        self.initial_state = []

        #penalization related
        self.full_area = geometric_params['geometric_config']['x_max'] * geometric_params['geometric_config']['y_max']
        self.linear_penalization = penalizations['linear']

        #keyword parameters for solver / bundle constructors
        self.sqs_params = {k: v for k, v in geometric_params.items() if k not in ("num_ellipses")}
        self.ellipse_bundle_params = {k: geometric_params[k] for k in ["geometric_config", "h", "num_ellipses"]}

    def raw_cost_function(self, state) -> float:
        '''Get full cost, including (linear) penalization.'''
        sqs = SquareSolver(**self.sqs_params)
        ellipses = EllipseBundle(**self.ellipse_bundle_params)
        n_param = 5 #ellipse parameter number
        try:
            for i in range(self.num_ellipses):
                idx = i * n_param
                ellipses.add(Ellipse(state[idx], state[idx + 1], state[idx + 2],
                                    state[idx + 3], state[idx + 4]))
            area_percent = (self.full_area - ellipses.area()) / self.full_area
            return sqs.solve(ellipses) + self.linear_penalization * area_percent
        except: #invalid configuration
            return float('inf')
    
    def get_neighbour(self, state):
        noise = [random.normalvariate(0.0, 1.0) * scale for scale in self.scales]
        return [s + n for s, n in zip(state, noise)]
    
    def raw_generation_function(self):
        '''Generates a random valid configuration, evaluates it, and returns cost & state.'''
        sqs = SquareSolver(**self.sqs_params)
        ellipses = EllipseBundle(**self.ellipse_bundle_params)
        
        ellipses.generate_random()
        
        if len(ellipses.bundle) < self.num_ellipses:
            return float('inf'), None
            
        state = []
        for e in ellipses.bundle:
            state.extend([e.center[0], e.center[1], e.quadratic_form[0,0], e.quadratic_form[0,1], e.quadratic_form[1,1]])
            
        try:
            area_percent = (self.full_area - ellipses.area()) / self.full_area
            cost = sqs.solve(ellipses) + self.linear_penalization * area_percent #OBJECTIVE FUNCTION
            return cost, state
        except Exception:
            return float('inf'), None
    
    def get_initial_state(self):
        while self.initial_state == []:
            try:
                with ProcessPool(max_workers=1, context=spawn_ctx) as pool:
                    future = pool.schedule(self.raw_generation_function, timeout=2.0)
                    cost, state = future.result() # type: ignore
                    if cost != float('inf'):
                        self.initial_state = state
            except Exception:
                pass # Ignore timeouts and crashes, try again


class EllipseDSA(DirectSimulatedAnnealing):
    '''Direct SA algorithm'''
    def __init__(self, scales, geometric_params : dict, penalizations : dict, **kwargs):
        super().__init__(**kwargs)
        self.geometric_params = geometric_params
        self.num_ellipses = geometric_params['num_ellipses']
        self.scales = scales * self.num_ellipses

        #penalization related
        self.full_area = geometric_params['geometric_config']['x_max'] * geometric_params['geometric_config']['y_max']
        self.linear_penalization = penalizations['linear']

        #keyword parameters for solver / bundle constructors
        self.sqs_params = {k: v for k, v in geometric_params.items() if k not in ("num_ellipses")}
        self.ellipse_bundle_params = {k: geometric_params[k] for k in ["geometric_config", "h", "num_ellipses"]}
    
    def raw_cost_function(self, state) -> float:
        '''Get full cost, including (linear) penalization.'''
        sqs = SquareSolver(**self.sqs_params)
        ellipses = EllipseBundle(**self.ellipse_bundle_params)
        n_param = 5 #ellipse parameter number
        try:
            for i in range(self.num_ellipses):
                idx = i * n_param
                ellipses.add(Ellipse(state[idx], state[idx + 1], state[idx + 2],
                                    state[idx + 3], state[idx + 4]))
            area_percent = (self.full_area - ellipses.area()) / self.full_area
            return sqs.solve(ellipses) + self.linear_penalization * area_percent #OBJECTIVE FUNCTION
        except: #invalid configuration
            return float('inf')
    
    def get_neighbour(self, state):
        noise = [random.normalvariate(0.0, 1.0) * scale for scale in self.scales]
        return [s + n for s, n in zip(state, noise)]
    
    def raw_generation_function(self):
        '''Generates a random valid configuration, evaluates it, and returns cost & state.'''
        sqs = SquareSolver(**self.sqs_params)
        ellipses = EllipseBundle(**self.ellipse_bundle_params)
        
        ellipses.generate_random()
        
        if len(ellipses.bundle) < self.num_ellipses:
            return float('inf'), None
            
        state = []
        for e in ellipses.bundle:
            state.extend([e.center[0], e.center[1], e.quadratic_form[0,0], e.quadratic_form[0,1], e.quadratic_form[1,1]])
            
        try:
            area_percent = (self.full_area - ellipses.area()) / self.full_area
            cost = sqs.solve(ellipses) + self.linear_penalization * area_percent #OBJECTIVE FUNCTION
            return cost, state
        except Exception:
            return float('inf'), None

    def get_starting_configs(self):
        self.configurations = []
        print(f"Generating {self.num_configs} starting configurations...")
        while len(self.configurations) < self.num_configs:
            try:
                with ProcessPool(max_workers=1, context=spawn_ctx) as pool:
                    future = pool.schedule(self.raw_generation_function, timeout=2.0)
                    cost, state = future.result() # type: ignore
                    if cost != float('inf') and state is not None:
                        self.configurations.append((cost, state))
            except Exception:
                pass # Ignore timeouts and crashes, try again
        self.configurations.sort()

    def test_generation(self):
        self.get_starting_configs()
        print(len(self.configurations), self.configurations)


def optimization_worker_SA(queue : mp.Queue, scales : list, geometric_params : dict, 
                        penalizations : dict, initial_params : list, kwargs_SA : dict):
    '''Defines the routine of an optimization worker'''
    solver = EllipseSA(scales, geometric_params, penalizations, **kwargs_SA)
    
    if initial_params == []:
        solver.get_initial_state()
        actual_initial_params = solver.initial_state
    else:
        actual_initial_params = initial_params
    print(actual_initial_params)
    start_time = time.time()
    best_param, best_cost = solver.run(actual_initial_params)
    runtime = time.time() - start_time

    # Get all the values into the SQL queue
    run_parameters = kwargs_SA.copy()
    run_parameters['runtime'] = runtime
    run_parameters['best_param'] = best_param
    run_parameters['best_cost'] = best_cost
    run_parameters['initial_params'] = actual_initial_params
    run_parameters['scales'] = scales

    for pen_type in penalizations:
        run_parameters[pen_type + '_penalization'] = penalizations[pen_type]
    run_parameters |= geometric_params

    #Convert all non-number values into JSON strings
    run_parameters = unfold_parameters(run_parameters)

    queue.put(run_parameters)

def optimization_worker_DSA(queue : mp.Queue, scales : list, geometric_params : dict, 
                        penalizations : dict, kwargs_DSA : dict):
    '''Defines the routine of a DSA optimization worker'''
    solver = EllipseDSA(scales, geometric_params, penalizations, **kwargs_DSA)
    start_time = time.time()
    best_cost, best_param = solver.run() # DSA returns (cost, state)
    runtime = time.time() - start_time

    # Get all the values into the SQL queue
    run_parameters = kwargs_DSA.copy()
    run_parameters['runtime'] = runtime
    run_parameters['best_param'] = best_param
    run_parameters['best_cost'] = best_cost
    run_parameters['initial_params'] = [state for cost, state in solver.configurations]
    run_parameters['scales'] = scales

    for pen_type in penalizations:
        run_parameters[pen_type + '_penalization'] = penalizations[pen_type]
    run_parameters |= geometric_params

    #Convert all non-number values into JSON strings
    run_parameters = unfold_parameters(run_parameters)

    queue.put(run_parameters)

def run_experiments(worker_target, geometric_params, penalizations, kwargs_optimization, scales, extra_worker_args=(), db_path='experiments.db', table_name='results', max_processes=10):
    config_bundle = {
        'geometric_params': geometric_params,
        'penalizations': penalizations,
        'kwargs_optimization': kwargs_optimization
    }
    
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
    
    combinations = []
    for combo in itertools.product(*sweep_lists):
        new_bundle = copy.deepcopy(config_bundle)
        for path, val in zip(paths, combo):
            target = new_bundle
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = val
        combinations.append(new_bundle)
        
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
    
    spawn_ctx = mp.get_context('spawn')
    queue = spawn_ctx.Queue()
    db_writer = spawn_ctx.Process(target=database_writer, args=(queue, column_names, db_path, table_name))
    db_writer.start()

    active_workers = []

    for combo_bundle in combinations:
        if len(active_workers) >= max_processes:
            sentinels = [w.sentinel for w in active_workers]
            wait(sentinels)
            active_workers = [w for w in active_workers if w.is_alive()]

        w = spawn_ctx.Process(target=worker_target, 
                              args=(queue, scales, combo_bundle['geometric_params'], 
                              combo_bundle['penalizations'], *extra_worker_args, combo_bundle['kwargs_optimization']))
        w.start()
        active_workers.append(w)
    
    for w in active_workers:
        w.join()
    
    queue.put(None)
    db_writer.join()
    print('Done')


if __name__ == '__main__':
    geometric_params = {
        "geometric_config": {'x_max': [1.0, 2.0], 'y_max': [1.0, 2.0], 'MW_x': 0.3, 'ME_x': 0.7},
        "h": 0.02,
        "heat_sources": 10.0,
        "base_temp": 0.0,
        "num_ellipses": 2
    }
    
    penalizations ={
        'linear' : 256.0
    }

    kwargs_SA = {
        'initial_temp' : 100,
        'min_temp' : 0.01,
        'cooling_rate' : 0.99
    }
    
    kwargs_DSA = {
        'initial_temp' : 1000.0, 'min_temp' : 0.0001,
        'cooling_rate_max' : 0.999, 'cooling_rate_min' : 0.99,
        'base_markov_length' : 50, 'num_configs' : 10,
        'initial_perturbation' : 0.05, 'min_perturbation' : 0.01, 
        'min_cost_gap' : 0.0001
    }

    #set up optimizations
    #initial_params = [0.5, 0.5, 25.0, 0, 25.0]
    initial_params = []
    scales = [0.1, 0.1, 5.0, 10.0, 5.0]

    run_experiments(optimization_worker_SA, geometric_params, penalizations, kwargs_SA, scales, extra_worker_args=(initial_params,), db_path = 'test.db', table_name='results', max_processes=10)
    

    '''
    Geometric info (constant in every optimization) : geometric_params
    Penalization info : penalizations
    Any SA/DSA info : best_params, best_cost, runtime
    SA specific info : kwargs_SA, initial_params
    DSA specific info : kwargs_DSA, ... (WIP)
    '''





    
