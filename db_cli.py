import sqlite3
import json
import os
import time
import argparse
import sys
import re
from datetime import datetime
import multiprocessing as mp

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse as MplEllipse
from pebble import ProcessPool, ProcessExpired
from concurrent.futures import TimeoutError

spawn_ctx = mp.get_context('spawn')

# =============================================================================
# 1. AREA & PENALTY ANALYSIS
# =============================================================================

def analyze_results(db_path, table_name, plot=False):
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    try:
        cursor.execute(f"SELECT run_id, best_param, best_cost, linear_penalization, x_max, y_max, best_param, initial_params FROM {table_name} ORDER BY linear_penalization ASC, best_cost ASC")
        rows = cursor.fetchall()

        print(f"{'Run ID':<8} | {'Total Cost':<12} | {'Penalization':<12} | {'Raw Cost':<12} | {'Area % Remaining':<18}")
        print("-" * 75)

        best_per_penalty = {}
        n_runs = 0
        
        for row in rows:
            if row['best_param'] == row['initial_params']:
                continue

            n_runs += 1
            run_id = row['run_id']
            total_cost = row['best_cost']
            penalty_factor = row['linear_penalization']
            
            full_area = row['x_max'] * row['y_max']   
            best_param = json.loads(row['best_param'])
            ellipses_area = 0
            
            for i in range(0, len(best_param), 5):
                A, B, C = best_param[i+2], best_param[i+3], best_param[i+4]
                det = A * C - B**2
                ellipses_area += np.pi / np.sqrt(det)
            
            area_percent = (full_area - ellipses_area) / full_area
            raw_cost = total_cost - penalty_factor * area_percent

            print(f"{run_id:<8} | {total_cost:<12.6f} | {penalty_factor:<12.2f} | {raw_cost:<12.10f} | {f'{area_percent*100:.4f}%':<18}")

            if plot:
                if penalty_factor not in best_per_penalty or total_cost < best_per_penalty[penalty_factor][0]:
                    best_per_penalty[penalty_factor] = (total_cost, area_percent)
        
        print("-" * 75)
        print(f'Total runs: {n_runs}')

        if plot:
            sorted_penalties = sorted(best_per_penalty.keys())
            areas = [best_per_penalty[p][1] for p in sorted_penalties]

            plt.figure(figsize=(10, 6))
            plt.plot(sorted_penalties, areas, marker='o', linestyle='-')
            plt.xlabel('Penalty Factor')
            plt.ylabel('Area % Remaining')
            plt.title(f'Penalty Factor vs Area % Remaining ({table_name})')
            plt.grid(True, which="both", ls="-", alpha=0.5)
            plt.ylim(0, 1)
            plt.show()
    except sqlite3.OperationalError as e:
        print(f"Error accessing database: {e}")
    finally:
        conn.close()

def calculate_ellipse_areas_from_db(db_path, table_name):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    try:
        query = f"SELECT run_id, best_param, x_max, y_max, linear_penalization, best_cost FROM {table_name} WHERE best_param IS NOT NULL ORDER BY linear_penalization ASC, best_cost ASC"
        cursor.execute(query)
        rows = cursor.fetchall()
        
        print(f"{'Run ID':<10} | {'Penalization':<15} | {'Total Ellipses Area':<20} | {'Unfilled Area %':<15}")
        print("-" * 72)
        
        for run_id, best_param_json, x_max, y_max, linear_penalization, best_cost in rows:
            try:
                best_param = json.loads(best_param_json)
            except json.JSONDecodeError:
                print(f"{run_id:<10} | Invalid JSON format")
                continue
            
            total_area = 0.0
            rect_area = x_max * y_max
            
            for i in range(0, len(best_param), 5):
                if i + 4 < len(best_param):
                    A, B, C = best_param[i+2], best_param[i+3], best_param[i+4]
                    det = (A * C) - (B ** 2)
                    if det > 0:
                        total_area += np.pi / np.sqrt(det)
                    else:
                        print(f"Warning: Ellipse {i//5 + 1} in run {run_id} has a non-positive determinant.")
            
            unfilled_percent = ((rect_area - total_area) / rect_area) * 100 if rect_area > 0 else 0.0
            print(f"{run_id:<10} | {linear_penalization:<15.2f} | {total_area:<20.6f} | {f'{unfilled_percent:.2f}%':<15}")
    except sqlite3.OperationalError as e:
        print(f"Database error: {e}")
    finally:
        conn.close()

# =============================================================================
# 2. PLOTTING DOMAINS & ELLIPSES
# =============================================================================

def plot_domains(db_path, table_name):
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    try:
        cursor.execute(f"SELECT run_id, best_param, x_max, y_max, linear_penalization FROM {table_name} ORDER BY linear_penalization ASC, best_cost ASC")
        rows = cursor.fetchall()

        fig, axes = plt.subplots(4, 4, figsize=(10, 10))
        axes = axes.flatten()
        idx = 0
        for row in rows:
            if idx >= len(axes): break
            if row['x_max'] != 3.0: continue
            
            ax = axes[idx]
            idx += 1

            ax.set_xlim(0.0, row['x_max'])
            ax.set_ylim(0.0, row['y_max'])
            ax.set_aspect('equal')
            ax.grid(True, linestyle='--', alpha=0.6)
            
            colors = ['red', 'blue', 'green', 'purple', 'orange', 'cyan']
            best_param = json.loads(row['best_param'])

            for i in range(0, len(best_param), 5):
                xc, yc = best_param[i], best_param[i+1]
                A, B, C = best_param[i+2], best_param[i+3], best_param[i+4]
                
                M = np.array([[A, B], [B, C]])
                w, v = np.linalg.eigh(M)
                sort_idx = w.argsort()
                w, v = w[sort_idx], v[:, sort_idx]
                
                a, b = 1.0 / np.sqrt(w[0]), 1.0 / np.sqrt(w[1])
                theta = np.degrees(np.arctan2(v[1, 0], v[0, 0]))

                ellipse_patch = MplEllipse(xy=(xc, yc), width=2*a, height=2*b, angle=theta,
                                           edgecolor=colors[(i//5) % len(colors)], facecolor='none', lw=2)
                ax.add_patch(ellipse_patch)
                ax.plot(xc, yc, 'k+', markersize=5)

            ax.set_title(f"Run {row['run_id']} - Pen: {row['linear_penalization']}")
            ax.set_xlabel("x")
            ax.set_ylabel("y")

        plt.tight_layout()
        plt.show()
    except sqlite3.OperationalError as e:
        print(f"Error accessing database: {e}")
    finally:
        conn.close()

def plot_initial_and_best_params(db_path, table_name, run_id):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute(f"SELECT initial_params, best_param, x_max, y_max FROM {table_name} WHERE run_id = ?", (run_id,))
    row = cursor.fetchone()
    conn.close()
    
    if not row:
        print(f"No valid run found for run_id {run_id} in {table_name}.")
        return
        
    initial_params = json.loads(row[0])
    best_params = json.loads(row[1])
    x_max, y_max = float(row[2]), float(row[3])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    def draw_ellipses(ax, params, title):
        ax.set_xlim(0.0, x_max)
        ax.set_ylim(0.0, y_max)
        ax.set_aspect('equal')
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.set_title(title)
        colors = ['red', 'blue', 'green', 'purple', 'orange', 'cyan']
        
        for i in range(0, len(params), 5):
            xc, yc = params[i], params[i+1]
            A, B, C = params[i+2], params[i+3], params[i+4]
            M = np.array([[A, B], [B, C]])
            w, v = np.linalg.eigh(M)
            a, b = 1.0 / np.sqrt(w[0]), 1.0 / np.sqrt(w[1])
            theta = np.degrees(np.arctan2(v[1, 0], v[0, 0]))
            color = colors[(i // 5) % len(colors)]
            
            ellipse = MplEllipse(xy=(xc, yc), width=2*a, height=2*b, angle=theta,
                                 edgecolor=color, facecolor='none', lw=2)
            ax.add_patch(ellipse)
            ax.plot(xc, yc, marker='+', color=color, markersize=5)
            
    draw_ellipses(ax1, initial_params, "Initial Parameters")
    draw_ellipses(ax2, best_params, f"Best Parameters (Run {run_id})")
    
    plt.tight_layout()
    plt.show()

# =============================================================================
# 3. COST HISTORY & IMPROVEMENTS
# =============================================================================

def parse_run_ids(args_list):
    if not args_list:
        return []
    
    joined = " ".join(args_list)
    # Convert 'to' or ' - ' into a clean hyphenated format
    joined = re.sub(r'\s*to\s*', '-', joined, flags=re.IGNORECASE)
    joined = re.sub(r'\s*-\s*', '-', joined)
    
    # Convert 'by' or 'step' into a colon format bound to the range
    joined = re.sub(r'\s*(?:by|step)\s*', ':', joined, flags=re.IGNORECASE)
    
    # Split on commas or spaces
    tokens = re.split(r'[,\s]+', joined)
    result = set()
    
    for token in tokens:
        if not token:
            continue
        if '-' in token:
            step = 1
            if ':' in token:
                parts = token.split(':')
                if len(parts) == 2 and parts[1].isdigit():
                    step = int(parts[1])
                    token = parts[0]
                else:
                    raise ValueError(f"Invalid step format: '{token}'")
            
            parts = token.split('-')
            if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                result.update(range(int(parts[0]), int(parts[1]) + 1, step))
            else:
                raise ValueError(f"Invalid range format: '{token}'")
        elif token.isdigit():
            result.add(int(token))
        else:
            raise ValueError(f"Invalid format: '{token}'. Expected a number or range.")
            
    return sorted(list(result))

def plot_cost_history(db_path, table_name, run_ids, penalty=None, save_path=None):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute(f"PRAGMA table_info({table_name})")
    columns = [row[1] for row in cursor.fetchall()]
    
    if 'track_file_name' not in columns:
        print(f"Error: Column 'track_file_name' does not exist in table '{table_name}'.")
        conn.close()
        return
        
    if penalty is not None:
        if 'linear_penalization' not in columns:
            print(f"Error: Column 'linear_penalization' does not exist in table '{table_name}'.")
            conn.close()
            return

        cursor.execute(f"SELECT DISTINCT run_id FROM {table_name} WHERE track_file_name IS NOT NULL AND linear_penalization = ? ORDER BY run_id ASC", (penalty,))
        penalty_runs = [row[0] for row in cursor.fetchall()]
        
        if not run_ids:
            run_ids = penalty_runs
        else:
            run_ids = sorted(list(set(run_ids).intersection(penalty_runs)))
            
        if not run_ids:
            print(f"No runs found matching the given criteria (penalty = {penalty}).")
            conn.close()
            return
    elif not run_ids:
        cursor.execute(f"SELECT DISTINCT run_id FROM {table_name} WHERE track_file_name IS NOT NULL ORDER BY run_id ASC")
        run_ids = [row[0] for row in cursor.fetchall()]
        if not run_ids:
            print(f"No runs with tracking files found in table '{table_name}'.")
            conn.close()
            return

    plt.figure(figsize=(10, 6))
    for run_id in run_ids:
        cursor.execute(f"SELECT track_file_name FROM {table_name} WHERE run_id = ?", (run_id,))
        row = cursor.fetchone()
        if not row or not row[0] or not os.path.exists(row[0]): continue
            
        df = pd.read_csv(row[0])
        if 'cost' in df.columns: # SA
            current_costs = df['cost'].where(df['accepted']).ffill()
        else: # DSA
            current_costs = df['best_cost']

        best_costs = current_costs.cummin()
        temperatures = df['temp']

        p = plt.plot(temperatures, best_costs, linewidth=2, label=f'Best Cost (Run {run_id})')
        plt.plot(temperatures, current_costs, color=p[0].get_color(), alpha=0.3, linewidth=1)

    conn.close()
    
    title_run_str = ", ".join(map(str, run_ids))
    if len(title_run_str) > 40:
        title_run_str = f"{len(run_ids)} runs"
        
    penalty_str = f" (Penalty: {penalty})" if penalty is not None else ""
    plt.title(f"Optimization Cost History (Run IDs: {title_run_str}){penalty_str}")
    plt.xlabel("Temperature")
    plt.ylabel("Cost")
    plt.xscale('log')
    plt.gca().invert_xaxis()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path)
        plt.close()
    else:
        plt.show()

def print_solution_improvements(db_path, table_name):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute(f"PRAGMA table_info({table_name})")
    cols = [row[1] for row in cursor.fetchall()]
    if 'track_file_name' not in cols or 'best_cost' not in cols:
        print(f"Error: Required columns not found in '{table_name}'.")
        return

    cursor.execute(f"SELECT run_id, track_file_name, best_cost FROM {table_name} ORDER BY best_cost ASC")
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
            initial_cost = float(pd.read_csv(track_file_name, nrows=1)['cost'].iloc[0])
            improvement_pct = ((initial_cost - best_cost) / initial_cost) * 100 if initial_cost != 0 else 0.0
            print(f"{run_id:8} | {initial_cost:15.4f} | {best_cost:15.4f} | {improvement_pct:17.2f}%")
        except Exception:
            print(f"{run_id:8} | {'Error':>15} | {best_str} | {'N/A':>18}")

# =============================================================================
# 3.5. RECOMPUTE COST
# =============================================================================

def _recompute_worker(geometric_info, h, heat_sources, base_temp, num_ellipses, best_param):
    from square_solver import SquareSolver, EllipseBundle, Ellipse
    import time
    
    sqs = SquareSolver(geometric_info, h, heat_sources, base_temp)
    ellipses = EllipseBundle(geometric_info, h, num_ellipses)
    
    for i in range(num_ellipses):
        idx = i * 5
        if idx + 4 < len(best_param):
            ellipses.add(Ellipse(
                best_param[idx], best_param[idx+1], best_param[idx+2], best_param[idx+3], best_param[idx+4]
            ))
            
    big_area = geometric_info['x_max'] * geometric_info['y_max']
    percent_area = (big_area - ellipses.area()) / big_area
    
    start_time = time.perf_counter()
    raw_objective = sqs.solve(ellipses)
    elapsed = time.perf_counter() - start_time
    
    return raw_objective, percent_area, elapsed

def recompute_cost_from_scratch(db_path, table_name, run_id):
    """
    Recomputes the PDE solve from scratch for the best_param of a given run_id.
    Useful for verifying if the recorded best cost matches a reproducible evaluation.
    """
    try:
        from square_solver import SquareSolver, EllipseBundle, Ellipse
    except ImportError:
        print("Error: Could not import square_solver. Ensure the C++ pybind module is compiled and accessible.")
        return

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    
    try:
        cursor.execute(f"SELECT * FROM {table_name} WHERE run_id = ?", (run_id,))
        row = cursor.fetchone()
    except sqlite3.OperationalError as e:
        print(f"Database error: {e}")
        conn.close()
        return
        
    conn.close()

    if not row:
        print(f"Run ID {run_id} not found in table '{table_name}'.")
        return
        
    if not row['best_param']:
        print(f"Run ID {run_id} has no best_param.")
        return
        
    try:
        best_param = json.loads(row['best_param'])
    except json.JSONDecodeError:
        print(f"Invalid JSON format in best_param for Run ID {run_id}.")
        return
        
    keys = row.keys()
    
    # Recover geometric properties (with fallbacks if missing from DB schema)
    geometric_info = {
        'x_max': float(row['x_max']) if 'x_max' in keys else 1.0,
        'y_max': float(row['y_max']) if 'y_max' in keys else 1.0,
        'MW_x': float(row['MW_x']) if 'MW_x' in keys else 0.3,
        'ME_x': float(row['ME_x']) if 'ME_x' in keys else 0.7
    }
    
    h = float(row['h']) if 'h' in keys else 0.02
    heat_sources = float(row['heat_sources']) if 'heat_sources' in keys else 10.0
    base_temp = float(row['base_temp']) if 'base_temp' in keys else 0.0
    num_ellipses = int(row['num_ellipses']) if 'num_ellipses' in keys else len(best_param) // 5
    penalization = float(row['linear_penalization']) if 'linear_penalization' in keys else 0.0
    stored_cost = float(row['best_cost']) if 'best_cost' in keys and row['best_cost'] is not None else None

    print(f"Recomputing cost for Run ID {run_id} (h={h}, penalty={penalization})...")
    
    try:
        with ProcessPool(max_workers=1, context=spawn_ctx) as pool:
            future = pool.schedule(
                _recompute_worker,
                args=(geometric_info, h, heat_sources, base_temp, num_ellipses, best_param), # type: ignore
                timeout=3.0
            )
            raw_objective, percent_area, elapsed = future.result()
        
        total_cost = raw_objective + (percent_area * penalization)
        
        print("-" * 50)
        print(f"Raw PDE Objective:    {raw_objective:.6f}")
        print(f"Unfilled Area %:      {percent_area * 100:.4f}%")
        print(f"Compute Time:         {elapsed:.4f} s")
        print("-" * 50)
        if stored_cost is not None:
            print(f"Stored Best Cost:     {stored_cost:.6f}")
        print(f"Recomputed Cost:      {total_cost:.6f}")
        
        if stored_cost is not None:
            diff = abs(stored_cost - total_cost)
            print(f"Difference:           {diff:.6e}")
            
    except TimeoutError:
        print("Failed during recomputation: Solver timed out.")
    except ProcessExpired:
        print("Failed during recomputation: Solver crashed (ProcessExpired).")
    except Exception as e:
        print(f"Failed during recomputation: {e}")

# =============================================================================
# 4. WATCH OPTIMIZATION
# =============================================================================

def watch(db_path, refresh_rate=2.0, max_age_minutes=10):
    max_age_seconds = max_age_minutes * 60
    while True:
        active_runs = []
        if os.path.exists(db_path):
            try:
                conn = sqlite3.connect(db_path)
                cursor = conn.cursor()
                cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
                tables = [r[0] for r in cursor.fetchall() if not r[0].startswith("sqlite_")]
                for table in tables:
                    cursor.execute(f"PRAGMA table_info({table})")
                    cols = [r[1] for r in cursor.fetchall()]
                    if 'track_file_name' in cols and 'runtime' in cols:
                        cursor.execute(f"SELECT run_id, track_file_name FROM {table} WHERE track_file_name IS NOT NULL AND runtime IS NULL")
                        for run_id, track_file in cursor.fetchall():
                            if track_file and os.path.exists(track_file):
                                mtime = os.path.getmtime(track_file)
                                if time.time() - mtime < max_age_seconds:
                                    active_runs.append({'table': table, 'run_id': run_id, 'file': track_file})
            except sqlite3.Error: pass
            finally: conn.close() #type: ignore
        
        os.system('cls' if os.name == 'nt' else 'clear')
        print(f"=== Optimization Watcher ===\nDB: {db_path} | Time: {datetime.now().strftime('%H:%M:%S')}")
        print("-" * 88)
        print(f"{'Table':<26} | {'Run':<5} | {'Iter':<6} | {'Temp':<10} | {'Cur Cost':<12} | {'Best Cost':<12}")
        print("-" * 88)
        
        if not active_runs: 
            print("No active runs found.")
        else:
            for run in sorted(active_runs, key=lambda x: x['run_id']):
                try:
                    df = pd.read_csv(run['file'])
                    if len(df) == 0: continue
                    last_row = df.iloc[-1]
                    
                    temp = last_row.get('temp', 0.0)
                    if 'best_cost' in df.columns:
                        cur_cost = best_cost = last_row.get('best_cost', 0.0)
                    else:
                        cur_cost = last_row.get('cost', 0.0)
                        best_cost = df[df['is_best'] == True]['cost'].iloc[-1] if 'is_best' in df.columns else df['cost'].min()
                    
                    print(f"{run['table']:<26} | {run['run_id']:<5} | {len(df):<6} | {temp:<10.4f} | {cur_cost:<12.4f} | {best_cost:<12.4f}")
                except Exception:
                    print(f"{run['table']:<26} | {run['run_id']:<5} | Error reading tracking file")
        try: time.sleep(refresh_rate)
        except KeyboardInterrupt: break

# =============================================================================
# 5. INTERACTIVE CLI MODE
# =============================================================================

def interactive_mode():
    print("=== Shape Optimizer DB CLI (Interactive Mode) ===")
    db_path = input("Enter database path [experiments.db]: ").strip() or "experiments.db"
    
    while True:
        if not os.path.exists(db_path):
            print(f"Database '{db_path}' not found.")
            db_path = input("Enter database path [experiments.db]: ").strip() or "experiments.db"
            if not os.path.exists(db_path): return
        
        try:
            conn = sqlite3.connect(db_path)
            cursor = conn.cursor()
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%';")
            tables = [r[0] for r in cursor.fetchall()]
            conn.close()
        except sqlite3.Error as e:
            print(f"Failed to read database: {e}")
            return

        if not tables:
            print(f"No valid tables found in '{db_path}'.")
            return
            
        print(f"\n--- Database: {db_path} ---")
        print("Available tables:")
        for i, t in enumerate(tables, 1):
            print(f"  {i}. {t}")
        
        t_idx = input("\nSelect table number (or 'q' to quit): ").strip()
        if t_idx.lower() == 'q': break
        if not t_idx.isdigit() or not (1 <= int(t_idx) <= len(tables)):
            print("Invalid selection.")
            continue
        
        table = tables[int(t_idx)-1]

        while True:
            print(f"\n--- Table: {table} ---")
            actions = [
                ("Analyze area/penalization results", "analyze"),
                ("Calculate unfilled area percentages", "area"),
                ("Plot a grid of best solution domains", "plot-domains"),
                ("Compare initial vs. best ellipses", "compare-params"),
                ("Plot cost tracking history", "cost-history"),
                ("Print table of cost improvements", "improvements"),
                ("Recompute cost from scratch for a run", "recompute-cost"),
                ("Continuously monitor active optimizations", "watch"),
                ("Change table or database", "back"),
                ("Quit", "exit")
            ]
            for i, (desc, _) in enumerate(actions, 1):
                print(f"  {i}. {desc}")
            
            a_idx = input("\nSelect action number: ").strip()
            if not a_idx.isdigit() or not (1 <= int(a_idx) <= len(actions)):
                print("Invalid action.")
                continue
            
            action = actions[int(a_idx)-1][1]
            if action == "exit": return
            elif action == "back": break
            
            print("\n" + "="*75)
            try:
                match action:
                    case "analyze":
                        plot = input("Plot results? (y/n) [n]: ").strip().lower() == 'y'
                        analyze_results(db_path, table, plot)
                    case "area":
                        calculate_ellipse_areas_from_db(db_path, table)
                    case "plot-domains":
                        plot_domains(db_path, table)
                    case "compare-params":
                        run_id = input("Enter Run ID: ").strip()
                        if run_id.isdigit(): plot_initial_and_best_params(db_path, table, int(run_id))
                        else: print("Invalid Run ID.")
                    case "cost-history":
                        run_ids_str = input("Enter Run IDs (e.g., '1 2 3', '5 to 50 by 5') [all]: ").strip()
                        penalty_str = input("Filter by penalty factor? (Enter value or blank): ").strip()
                        penalty = float(penalty_str) if penalty_str else None
                        parsed_run_ids = parse_run_ids(run_ids_str.split()) if run_ids_str else []
                        plot_cost_history(db_path, table, parsed_run_ids, penalty=penalty)
                    case "improvements":
                        print_solution_improvements(db_path, table)
                    case "recompute-cost":
                        run_id = input("Enter Run ID: ").strip()
                        if run_id.isdigit(): recompute_cost_from_scratch(db_path, table, int(run_id))
                        else: print("Invalid Run ID.")
                    case "watch":
                        refresh = float(input("Refresh rate in seconds [2.0]: ").strip() or 2.0)
                        max_age = int(input("Max age in minutes [10]: ").strip() or 10)
                        watch(db_path, refresh, max_age)
            except Exception as e: print(f"Error executing action: {e}")
            print("="*75)
            if action != "watch": input("\nPress Enter to continue...")

# =============================================================================
# CLI DEFINITION
# =============================================================================

def main():
    if len(sys.argv) == 1:
        try:
            interactive_mode()
        except KeyboardInterrupt:
            print("\nExiting interactive mode.")
        return

    parser = argparse.ArgumentParser(description="Master DB utility for Shape Optimizer")
    parser.add_argument("--db", type=str, default="experiments.db", help="Path to SQLite database (default: experiments.db)")
    
    subparsers = parser.add_subparsers(dest="command", required=True, help="Available subcommands")
    
    # Analyze
    p_analyze = subparsers.add_parser("analyze", help="Analyze area/penalization results")
    p_analyze.add_argument("table", type=str, help="Table name to analyze")
    p_analyze.add_argument("--plot", action="store_true", help="Plot the results")
    
    # Area
    p_area = subparsers.add_parser("area", help="Calculate unfilled area percentages")
    p_area.add_argument("table", type=str, help="Table name")

    # Domains
    p_domains = subparsers.add_parser("plot-domains", help="Plot a grid of best solution domains")
    p_domains.add_argument("table", type=str, help="Table name")
    
    # Params
    p_params = subparsers.add_parser("compare-params", help="Compare initial vs. best ellipses")
    p_params.add_argument("table", type=str, help="Table name")
    p_params.add_argument("run_id", type=int, help="Run ID to visualize")
    
    # History
    p_history = subparsers.add_parser("cost-history", help="Plot cost tracking history")
    p_history.add_argument("table", type=str, help="Table name")
    p_history.add_argument("run_ids", type=str, nargs="*", help="Run IDs to plot (e.g., '1 2 3', '5 to 50 by 5', '1, 5-10'). If omitted, plots all runs.")
    p_history.add_argument("-p", "--penalty", type=float, help="Filter runs by a specific linear_penalization value")
    
    # Improvements
    p_improv = subparsers.add_parser("improvements", help="Print table of cost improvements")
    p_improv.add_argument("table", type=str, help="Table name")
    
    # Recompute
    p_recompute = subparsers.add_parser("recompute-cost", help="Recompute cost from scratch for a run")
    p_recompute.add_argument("table", type=str, help="Table name")
    p_recompute.add_argument("run_id", type=int, help="Run ID to recompute")

    # Watch
    p_watch = subparsers.add_parser("watch", help="Continuously monitor active optimizations")
    p_watch.add_argument("--refresh", type=float, default=2.0, help="Refresh rate in seconds")
    p_watch.add_argument("--max-age", type=int, default=10, help="Max minutes since last file update")

    args = parser.parse_args()

    match args.command:
        case "analyze":
            analyze_results(args.db, args.table, args.plot)
        case "area":
            calculate_ellipse_areas_from_db(args.db, args.table)
        case "plot-domains":
            plot_domains(args.db, args.table)
        case "compare-params":
            plot_initial_and_best_params(args.db, args.table, args.run_id)
        case "cost-history":
            try:
                parsed_run_ids = parse_run_ids(args.run_ids)
                plot_cost_history(args.db, args.table, parsed_run_ids, penalty=args.penalty)
            except ValueError as e:
                print(f"Error: {e}")
        case "improvements":
            print_solution_improvements(args.db, args.table)
        case "recompute-cost":
            recompute_cost_from_scratch(args.db, args.table, args.run_id)
        case "watch":
            watch(args.db, args.refresh, args.max_age)

if __name__ == "__main__":
    main()
