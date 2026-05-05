import sqlite3
import json
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse as MplEllipse
import numpy as np

def get_params_from_db(db_path, table_name, run_id):
    """
    Given a path to a sqlite database, table name and run_id, 
    get the values of 'initial_params' and 'best_param' in list format,
    as well as the domain bounds x_max and y_max.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    query = f"SELECT initial_params, best_param, x_max, y_max FROM {table_name} WHERE run_id = ?"
    cursor.execute(query, (run_id,))
    row = cursor.fetchone()
    
    conn.close()
    
    if row:
        initial_params = json.loads(row[0])
        best_params = json.loads(row[1])
        x_max = float(row[2])
        y_max = float(row[3])
        return initial_params, best_params, x_max, y_max
    
    return None, None, None, None

def plot_initial_and_best_params(initial_params, best_params, x_max, y_max):
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
            
            a = 1.0 / np.sqrt(w[0])
            b = 1.0 / np.sqrt(w[1])
            theta = np.degrees(np.arctan2(v[1, 0], v[0, 0]))
            
            color = colors[(i // 5) % len(colors)]
            ellipse = MplEllipse(
                xy=(xc, yc), width=2*a, height=2*b, angle=theta,
                edgecolor=color, facecolor='none', lw=2
            )
            ax.add_patch(ellipse)
            ax.plot(xc, yc, marker='+', color=color, markersize=5)
            
    draw_ellipses(ax1, initial_params, "Initial Parameters")
    draw_ellipses(ax2, best_params, "Best Parameters")
    
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    db_path = 'experiments.db'
    table_name = 'results'
    run_id = 17
    plot_initial_and_best_params(*get_params_from_db(db_path, table_name, run_id))