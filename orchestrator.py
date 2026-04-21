from annealing import *
from square_solver import *
import numpy as np
import random
import multiprocessing as mp
from multiprocessing.pool import Pool
from multiprocessing.connection import wait
from pebble import concurrent
import sqlite3
import json
import time


class EllipseSA(BaseSimulatedAnnealing):
    '''Classical SA algorithm, needs an initial suitable configuration to progress with the optimization'''
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
            return sqs.solve(ellipses) + self.linear_penalization * area_percent
        except: #invalid configuration
            return float('inf')
    
    def get_neighbour(self, state):
        noise = [random.normalvariate(0.0, 1.0) * scale for scale in self.scales]
        return [s + n for s, n in zip(state, noise)]


def database_writer(queue : mp.Queue, column_names : dict, db_path = "experiments.db"):
    '''Given a queue, listens indefinitely and writes to the database when the queue is non empty.
       The function terminates when a None object is inserted in the queue.'''
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    #Table creation with the columns we want
    columns_definitions = ['run_id INTEGER PRIMARY KEY AUTOINCREMENT']
    for name, sql_type in column_names.items():
        if sql_type in ('REAL','INTEGER','TEXT'):
            columns_definitions.append(f'{name} {sql_type}')
        else: #assume is real
            columns_definitions.append(f'{name} REAL')
    
    sql_make_table = f'CREATE TABLE IF NOT EXISTS results ({', '.join(columns_definitions)})'
    cursor.execute(sql_make_table)
    conn.commit()

    while True:
        record = queue.get() #dict with column name and value or None

        if record is None: #end of workers
            break
        
        columns = list(record.keys())
        columns_sql = ', '.join(columns)
        placeholders_sql = ', '.join([f':{col}' for col in columns])

        insert_sql = f'INSERT INTO results ({columns_sql}) VALUES ({placeholders_sql})'
        cursor.execute(insert_sql, record)
        conn.commit()

    conn.close()


def get_column_names_dict(geometric_params : dict, penalizations : dict, kwargs_optimization : dict, paramtype = list):
    '''Get a dict where column names map to their respective SQL type. Between optimization specific arguments and
    geometric configuration, includes the columns 'runtime', 'best_param', 'best_cost', 'initial_params', 'scales'.
    '''
    # Sort which types are mapped to SQL types, default to 'TEXT'
    sql_type = {int : 'INTEGER', float : 'REAL'}
    def get_sql_type(val):
        if isinstance(val, type):
            return sql_type.get(val, 'TEXT')
        return sql_type.get( type(val), 'TEXT')
    
    clean_kwargs_optimization = {k : get_sql_type(v) for k, v in kwargs_optimization.items()}

    param_SQL_type = get_sql_type(paramtype)
    middle = {'runtime' : 'REAL', 'best_param' : param_SQL_type, 'best_cost' : 'REAL',
              'initial_params' : param_SQL_type, 'scales' : param_SQL_type}

    clean_penalizations = {k + '_penalization' : 'REAL' for k,v in penalizations.items()}

    flat_geo = {}
    for k, v in geometric_params.items():
        if isinstance(v, dict):
            flat_geo |= v
        else:
            flat_geo[k] = v
    
    clean_geometric_params = {k : get_sql_type(v) for k, v in flat_geo.items()}
    
    return clean_kwargs_optimization | middle | clean_penalizations | clean_geometric_params


def optimization_worker(queue : mp.Queue, scales : list, geometric_params : dict, 
                        penalizations : dict, initial_params : list, kwargs_SA : dict):
    '''Defines the routine of an optimization worker'''
    solver = EllipseSA(scales, geometric_params, penalizations, **kwargs_SA)
    start_time = time.time()
    best_param, best_cost = solver.run(initial_params)
    runtime = time.time() - start_time

    # Get all the values into the SQL queue
    run_parameters = kwargs_SA.copy()
    run_parameters['runtime'] = runtime
    run_parameters['best_param'] = best_param
    run_parameters['best_cost'] = best_cost
    run_parameters['initial_params'] = initial_params
    run_parameters['scales'] = scales

    for pen_type in penalizations:
        run_parameters[pen_type + '_penalization'] = penalizations[pen_type]
    run_parameters |= geometric_params

    #Convert all non-number values into JSON strings
    run_parameters = unfold_parameters(run_parameters)

    queue.put(run_parameters)

def unfold_parameters(parameters : dict):
    '''From a dictionary, unfolds dictionary values, if the value is a list or tuple outputs a JSON string. '''
    output = {}
    for k, v in parameters.items():
        if isinstance(v, dict):
            output |= v
        else:
            output[k] = v
    
    for key, val in output.items():
        if isinstance(val, (list, tuple)):
            output[key] = json.dumps(val)
    
    return output
    
if __name__ == '__main__':
    geometric_params = {
        "geometric_config": {'x_max': 3.0, 'y_max': 2.0, 'MW_x': 0.3, 'ME_x': 0.7},
        "h": 0.02,
        "heat_sources": 10.0,
        "base_temp": 0.0,
        "num_ellipses": 1
    }
    
    penalizations = {
        "linear" : 0.0,
    } #quadratic to be implemented

    kwargs_SA = {
        'initial_temp' : 1000,
        'min_temp' : 0.0001,
        'cooling_rate' : 0.99
    }

    #set up optimizations
    initial_params = [0.5, 0.5, 25.0, 0, 25.0]
    scales = [0.1, 0.1, 5.0, 25.0, 5.0] 


    #database setup
    column_names = get_column_names_dict(geometric_params = geometric_params, penalizations=penalizations,
                                         kwargs_optimization= kwargs_SA, paramtype=list)
    spawn_ctx = mp.get_context('spawn')
    queue = spawn_ctx.Queue()
    db_writer = spawn_ctx.Process(target=database_writer, args=(queue, column_names))
    db_writer.start()

    linear_penalization_list = [128.0]
    
    max_processes = 10
    active_workers = []

    for penalization in linear_penalization_list:
        if len(active_workers) >= max_processes:
            sentinels = [w.sentinel for w in active_workers]
            wait(sentinels)
            active_workers = [w for w in active_workers if w.is_alive()]

        task_penalizations = penalizations.copy()
        task_penalizations['linear'] = penalization
        
        w = spawn_ctx.Process(target=optimization_worker, 
                              args=(queue, scales, geometric_params, task_penalizations, initial_params, kwargs_SA))
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





    
