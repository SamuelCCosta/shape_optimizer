import sqlite3
import json
import multiprocessing as mp

def database_writer(queue: mp.Queue, column_names: dict, db_path="experiments.db", table_name="results"):
    '''Given a queue, listens indefinitely and writes to the database when the queue is non empty.
       The function terminates when a None object is inserted in the queue.'''
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Create table with just the primary key if it doesn't exist
    cursor.execute(f"CREATE TABLE IF NOT EXISTS {table_name} (run_id INTEGER PRIMARY KEY AUTOINCREMENT)")
    
    # Fetch existing columns to dynamically alter table if new configurations are introduced
    cursor.execute(f"PRAGMA table_info({table_name})")
    existing_columns = {row[1] for row in cursor.fetchall()}

    # Add any missing columns dynamically
    for name, sql_type in column_names.items():
        if name not in existing_columns:
            if sql_type not in ('REAL', 'INTEGER', 'TEXT'):
                sql_type = 'REAL'
            cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN {name} {sql_type}")
            
    conn.commit()

    while True:
        record = queue.get()  # dict with column name and value or None

        if record is None:  # end of workers
            break
        
        columns = list(record.keys())
        columns_sql = ', '.join(columns)
        placeholders_sql = ', '.join([f':{col}' for col in columns])

        insert_sql = f'INSERT INTO {table_name} ({columns_sql}) VALUES ({placeholders_sql})'
        cursor.execute(insert_sql, record)
        conn.commit()

    conn.close()

def get_column_names_dict(geometric_params: dict, penalizations: dict, kwargs_optimization: dict, paramtype=list):
    '''Get a dict where column names map to their respective SQL type. Between optimization specific arguments and
    geometric configuration, includes the columns 'runtime', 'best_param', 'best_cost', 'initial_params'.
    '''
    # Sort which types are mapped to SQL types, default to 'TEXT'
    sql_type = {int: 'INTEGER', float: 'REAL'}
    def get_sql_type(val):
        if isinstance(val, type):
            return sql_type.get(val, 'TEXT')
        return sql_type.get(type(val), 'TEXT')
    
    clean_kwargs_optimization = {k: get_sql_type(v) for k, v in kwargs_optimization.items()}

    param_SQL_type = get_sql_type(paramtype)
    middle = {'runtime': 'REAL', 'best_param': param_SQL_type, 'best_cost': 'REAL',
              'initial_params': param_SQL_type}

    clean_penalizations = {k + '_penalization': 'REAL' for k, v in penalizations.items()}

    flat_geo = {}
    for k, v in geometric_params.items():
        if isinstance(v, dict):
            flat_geo |= v
        else:
            flat_geo[k] = v
    
    clean_geometric_params = {k: get_sql_type(v) for k, v in flat_geo.items()}
    
    return (clean_kwargs_optimization | middle | clean_penalizations | 
            clean_geometric_params | {'cpu_model': 'TEXT'})

def unfold_parameters(parameters: dict):
    '''From a dictionary, unfolds dictionary values, if the value is a list or tuple outputs a JSON string.'''
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

def get_cpu_model():
    '''Returns the CPU model name on Linux systems.'''
    try:
        with open("/proc/cpuinfo", "r") as f:
            for line in f:
                if "model name" in line:
                    return line.split(":")[1].strip()
    except Exception:
        pass
    return "Unknown"
