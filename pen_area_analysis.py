import sqlite3
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse as MplEllipse

def analyze_results(db_path="experiments.db", table_name="results", plot = False):
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    try:
        # Select relevant columns from the results table
        cursor.execute(f"SELECT run_id, best_param, best_cost, linear_penalization, x_max, y_max, best_param, initial_params FROM {table_name} ORDER BY linear_penalization ASC")
        rows = cursor.fetchall()


        print(f"{'Run ID':<8} | {'Total Cost':<12} | {'Penalization':<12} | {'Raw Cost':<12} | {'Area % Remaining':<18}")
        print("-" * 75)

        penalties = []
        areas = []
        n_runs = 0
        for row in rows:
            #ensure the algorithm improved the solution
            if row['best_param'] == row['initial_params']:
                continue

            #variable filtering rule
            #if row['x_max'] != 3.0 or row['linear_penalization'] not in [0.0]+[float(2**i) for i in range(13)]:
            #    continue
            n_runs += 1
            run_id = row['run_id']
            total_cost = row['best_cost']
            penalty_factor = row['linear_penalization']
            
            # Calculate total area of the bounding box
            full_area = row['x_max'] * row['y_max']   
            best_param = json.loads(row['best_param'])
            ellipses_area = 0
            # Each ellipse has 5 parameters: [x, y, A, B, C]
            # Area of ellipse defined by Ax^2 + 2Bxy + Cy^2 = 1 is pi / sqrt(AC - B^2)
            for i in range(0, len(best_param), 5):
                A, B, C = best_param[i+2], best_param[i+3], best_param[i+4]
                det = A * C - B**2
                ellipses_area += np.pi / np.sqrt(det)
            
            area_percent = (full_area - ellipses_area) / full_area
            raw_cost = total_cost - penalty_factor * area_percent

            print(f"{run_id:<8} | {total_cost:<12.6f} | {penalty_factor:<12.1f} | {raw_cost:<12.10f} | {f'{area_percent*100:.4f}%':<18}")

            if plot:
                penalties.append(penalty_factor)
                areas.append(area_percent)
        
        print("-" * 75)
        print(f'Total runs: {n_runs}')

        if plot:
            log_scale = False
            plt.figure(figsize=(10, 6))
            plt.plot(penalties, areas, marker='o', linestyle='-')
            if log_scale:
                plt.xscale('log') # Logarithmic scale for penalty factor
                plt.xlabel('Penalty Factor (log scale)')
            else:
                plt.xlabel('Penalty Factor')
            plt.ylabel('Area % Remaining')
            plt.title('Penalty Factor vs Area % Remaining')
            plt.grid(True, which="both", ls="-", alpha=0.5)
            plt.show()
    

    except sqlite3.OperationalError as e:
        print(f"Error accessing database: {e}")
    finally:
        conn.close()


def plot_domains(db_path="experiments.db", table_name = "results"):
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    try:
        cursor.execute(f"SELECT run_id, best_param, x_max, y_max, linear_penalization FROM {table_name} ORDER BY linear_penalization ASC")
        rows = cursor.fetchall()

        fig, axes = plt.subplots(4, 4, figsize=(10, 10))
        axes = axes.flatten()
        idx = 0
        for row in rows:
            penalty = row['linear_penalization']
            run_id = row['run_id']
            x_max, y_max = row['x_max'], row['y_max']
            best_param = json.loads(row['best_param'])
            #filtering rule
            if x_max != 3.0:
                continue
            
            ax = axes[idx]
            idx += 1

            ax.set_xlim(0.0, x_max)
            ax.set_ylim(0.0, y_max)
            ax.set_aspect('equal')
            ax.grid(True, linestyle='--', alpha=0.6)
            
            colors = ['red', 'blue', 'green', 'purple', 'orange', 'cyan']

            # Process each ellipse (5 parameters: x, y, A, B, C)
            for i in range(0, len(best_param), 5):
                xc, yc = best_param[i], best_param[i+1]
                A, B, C = best_param[i+2], best_param[i+3], best_param[i+4]
                
                # Eigen decomposition to find axes and rotation
                M = np.array([[A, B], [B, C]])
                w, v = np.linalg.eigh(M)
                
                # Smallest eigenvalue -> Major axis
                sort_idx = w.argsort()
                w = w[sort_idx]
                v = v[:, sort_idx]
                
                a = 1.0 / np.sqrt(w[0])
                b = 1.0 / np.sqrt(w[1])
                
                # Orientation
                theta = np.degrees(np.arctan2(v[1, 0], v[0, 0]))

                ellipse_patch = MplEllipse(
                    xy=(xc, yc),
                    width=2*a,
                    height=2*b,
                    angle=theta,
                    edgecolor=colors[(i//5) % len(colors)],
                    facecolor='none',
                    lw=2
                )
                ax.add_patch(ellipse_patch)
                ax.plot(xc, yc, 'k+', markersize=5)

            ax.set_title(f"Run {run_id} - Penalty: {penalty}")
            ax.set_xlabel("x")
            ax.set_ylabel("y")

        plt.tight_layout()
        plt.show()

    except sqlite3.OperationalError as e:
        print(f"Error accessing database: {e}")
    finally:
        conn.close()



if __name__ == '__main__':
    analyze_results(db_path="experiments.db", table_name="x1y1n2lambda16_SA", plot=False)
    #plot_domains()