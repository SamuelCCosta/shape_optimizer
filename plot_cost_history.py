import sqlite3
import os
import pandas as pd
import matplotlib.pyplot as plt

def plot_cost_history(db_path, table_name, run_ids, save_path=None, min_iter=None, max_iter=None):
    """
    Fetches the track_file_name for given run_ids from the database,
    reads the corresponding CSV files, and plots the current and best cost in a single pane.
    If save_path is provided, saves the plot as an image instead of showing it.
    """
    if not isinstance(run_ids, (list, tuple)):
        run_ids = [run_ids]

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Check if 'track_file_name' column exists to avoid SQL errors
    cursor.execute(f"PRAGMA table_info({table_name})")
    columns = [row[1] for row in cursor.fetchall()]
    
    if 'track_file_name' not in columns:
        print(f"Error: Column 'track_file_name' does not exist in table '{table_name}'.")
        conn.close()
        return
        
    plt.figure(figsize=(10, 6))
    
    for run_id in run_ids:
        query = f"SELECT track_file_name FROM {table_name} WHERE run_id = ?"
        cursor.execute(query, (run_id,))
        row = cursor.fetchone()
        
        if not row or not row[0]:
            print(f"No valid track_file_name found for run_id {run_id} in {table_name}.")
            continue
            
        track_file_name = row[0]
        if not os.path.exists(track_file_name):
            print(f"Error: CSV file not found at path: {track_file_name}")
            continue
            
        # Read the tracked history using pandas
        df = pd.read_csv(track_file_name)
        
        if 'accepted' in df.columns:
            current_costs = df['cost'].where(df['accepted']).ffill()
        else:
            current_costs = df['cost']
            
        best_costs = current_costs.cummin()
        iterations = df.index
        
        mask = pd.Series(True, index=iterations)
        if min_iter is not None:
            mask &= (iterations >= min_iter)
        if max_iter is not None:
            mask &= (iterations <= max_iter)
            
        iterations = iterations[mask]
        best_costs = best_costs[mask]
        current_costs = current_costs[mask]

        # Plot best cost and current cost with matching colors
        p = plt.plot(iterations, best_costs, linewidth=2, label=f'Best Cost (Run {run_id})')
        plt.plot(iterations, current_costs, color=p[0].get_color(), alpha=0.3, linewidth=1)
        #plt.plot(iterations, current_costs, linewidth=1, label=f'Current Cost (Run {run_id})')

    conn.close()
    
    run_ids_str = ", ".join(map(str, run_ids))
    plt.title(f"Optimization Cost History (Run IDs: {run_ids_str})")
    plt.xlabel("Iteration")
    plt.ylabel("Cost")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path)
        plt.close()
    else:
        plt.show()

def print_solution_improvements(db_path, table_name):
    """
    Reads the initial cost from the CSV and the best cost from the DB
    for all runs in the given table, and calculates the improvement.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Check if 'track_file_name' and 'best_cost' columns exist
    cursor.execute(f"PRAGMA table_info({table_name})")
    columns = [row[1] for row in cursor.fetchall()]
    
    if 'track_file_name' not in columns or 'best_cost' not in columns:
        print(f"Error: Required columns not found in table '{table_name}'.")
        conn.close()
        return

    query = f"SELECT run_id, track_file_name, best_cost FROM {table_name}"
    cursor.execute(query)
    rows = cursor.fetchall()
    conn.close()

    print(f"{'Run ID':>8} | {'Initial Cost':>15} | {'Best Cost':>15} | {'Improvement (%)':>18}")
    print("-" * 65)

    for run_id, track_file_name, best_cost in rows:
        best_str = f"{best_cost:15.4f}" if best_cost is not None else f"{'N/A':>15}"
        
        if not track_file_name or not os.path.exists(track_file_name):
            print(f"{run_id:8} | {'N/A':>15} | {best_str} | {'N/A':>18}")
            continue
            
        try:
            # Read just the first row of the CSV to get the initial cost
            df = pd.read_csv(track_file_name, nrows=1)
            initial_cost = float(df['cost'].iloc[0])
            
            improvement_pct = ((initial_cost - best_cost) / initial_cost) * 100 if initial_cost != 0 else 0.0
            
            print(f"{run_id:8} | {initial_cost:15.4f} | {best_cost:15.4f} | {improvement_pct:17.2f}%")
        except Exception:
            print(f"{run_id:8} | {'Error':>15} | {best_str} | {'N/A':>18}")

def plot_all_tables_cost_history(db_path, output_folder, min_iter=None, max_iter=None):
    """
    Iterates through all non-sqlite tables in the given database,
    applies plot_cost_history for all valid run_ids in each table,
    and saves the plots as PNGs in the specified output folder.
    """
    os.makedirs(output_folder, exist_ok=True)
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    try:
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = cursor.fetchall()
        
        for (table_name,) in tables:
            if table_name.startswith("sqlite_"):
                continue
            
            cursor.execute(f"PRAGMA table_info({table_name})")
            columns = [row[1] for row in cursor.fetchall()]
            
            if 'track_file_name' not in columns or 'run_id' not in columns:
                continue
                
            cursor.execute(f"SELECT run_id FROM {table_name} WHERE track_file_name IS NOT NULL")
            run_ids = [row[0] for row in cursor.fetchall()]
            
            if run_ids:
                save_path = os.path.join(output_folder, f"{table_name}_cost_history.png")
                print(f"Saving plot for table '{table_name}' -> {save_path}")
                plot_cost_history(db_path, table_name, run_ids, save_path=save_path, min_iter=min_iter, max_iter=max_iter)
    except sqlite3.Error as e:
        print(f"Database error occurred: {e}")
    finally:
        conn.close()

if __name__ == '__main__':
    # Example Usage
    db_path = 'jasminum-remote.db'
    table_name = 'results'
    run_ids = [i for i in range(1,10)]
    run_ids = [7, 10] #7 , 10
    run_ids = [i for i in range(1,31)]
    min_iter = None
    max_iter = None

    #print_solution_improvements(db_path, table_name)
    plot_cost_history(db_path, table_name, run_ids, min_iter=min_iter, max_iter=max_iter)
    #plot_all_tables_cost_history(db_path, 'plots/cost_history')
    
