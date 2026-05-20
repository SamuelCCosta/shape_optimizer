import os
import time
import pandas as pd
from datetime import datetime
import argparse
import sqlite3

def get_active_runs(db_path, max_age_seconds=600):
    """Find active runs by querying the database and checking recent file activity."""
    active_runs = []
    if not os.path.exists(db_path):
        return active_runs
        
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = [row[0] for row in cursor.fetchall() if not row[0].startswith("sqlite_")]
        
        for table in tables:
            cursor.execute(f"PRAGMA table_info({table})")
            cols = [row[1] for row in cursor.fetchall()]
            
            if 'track_file_name' in cols and 'runtime' in cols:
                # An active run has NOT finished yet (runtime IS NULL)
                cursor.execute(f"SELECT run_id, track_file_name FROM {table} WHERE track_file_name IS NOT NULL AND runtime IS NULL")
                for run_id, track_file in cursor.fetchall():
                    if track_file and os.path.exists(track_file):
                        try:
                            mtime = os.path.getmtime(track_file)
                            # Hide crashed runs that haven't updated their file recently
                            if time.time() - mtime < max_age_seconds:
                                active_runs.append({
                                    'table': table,
                                    'run_id': run_id,
                                    'file': track_file,
                                    'mtime': mtime
                                })
                        except OSError:
                            pass
    except sqlite3.Error:
        pass
    finally:
        conn.close() #type:ignore
                    
    return active_runs

def watch(db_path, refresh_rate=2.0, max_age_minutes=10):
    """Continuously prints out the progress of active optimization runs."""
    max_age_seconds = max_age_minutes * 60
    
    while True:
        active_runs = get_active_runs(db_path, max_age_seconds)
        
        os.system('cls' if os.name == 'nt' else 'clear')
        print(f"=== Optimization Watcher ===")
        print(f"Database: {db_path} | Time: {datetime.now().strftime('%H:%M:%S')}")
        print(f"Showing runs modified in the last {max_age_minutes} minutes.")
        print("-" * 88)
        print(f"{'Table':<26} | {'Run':<5} | {'Iter':<6} | {'Temp':<10} | {'Cur Cost':<12} | {'Best Cost':<12}")
        print("-" * 88)
        
        if not active_runs:
            print("No active runs found.")
        else:
            # Sort by most recently modified
            for run in sorted(active_runs, key=lambda x: x['run_id']):
                try:
                    df = pd.read_csv(run['file'])
                    if len(df) == 0:
                        continue
                    
                    iteration = len(df)
                    last_row = df.iloc[-1]
                    
                    if 'best_cost' in df.columns:
                        # Direct Simulated Annealing tracking
                        temp = last_row.get('temp', 0.0)
                        best_cost = last_row.get('best_cost', 0.0)
                        current_cost = best_cost
                    else:
                        # Standard Simulated Annealing tracking
                        temp = last_row.get('temp', 0.0)
                        current_cost = last_row.get('cost', 0.0)
                        best_cost = df[df['is_best'] == True]['cost'].iloc[-1] if 'is_best' in df.columns else df['cost'].min()
                    
                    print(f"{run['table']:<26} | {run['run_id']:<5} | {iteration:<6} | {temp:<10.4f} | {current_cost:<12.4f} | {best_cost:<12.4f}")
                except Exception:
                    print(f"{run['table']:<26} | {run['run_id']:<5} | Error reading tracking file")
        try:
            time.sleep(refresh_rate)
        except KeyboardInterrupt:
            print("\nExiting watcher.")
            break

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Watch active optimization runs.")
    parser.add_argument("db_file", help="Path to the database file (e.g., experiments.db)")
    parser.add_argument("-r", "--refresh-rate", type=float, default=2.0, help="Refresh rate in seconds (default: 2.0)")
    parser.add_argument("-m", "--max-age", type=int, default=5, help="Max age in minutes to consider a run active (default: 5)")
    
    args = parser.parse_args()
    watch(args.db_file, refresh_rate=args.refresh_rate, max_age_minutes=args.max_age)
