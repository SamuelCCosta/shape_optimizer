import sqlite3
import json
import numpy as np

def calculate_ellipse_areas_from_db(db_path, table_name):
    """
    Extracts ellipse parameters from a SQLite database and calculates
    the total area of the ellipses for each row.
    
    Parameters:
    - db_path: Path to the SQLite database file.
    - table_name: Name of the table containing the 'best_param' column.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    try:
        # Fetch run_id, best_param, x_max, y_max, and linear_penalization from the table
        query = f"SELECT run_id, best_param, x_max, y_max, linear_penalization FROM {table_name} WHERE best_param IS NOT NULL ORDER BY linear_penalization ASC, run_id ASC"
        cursor.execute(query)
        rows = cursor.fetchall()
        
        print(f"{'Run ID':<10} | {'Penalization':<15} | {'Total Ellipses Area':<20} | {'Unfilled Area %':<15}")
        print("-" * 72)
        
        for run_id, best_param_json, x_max, y_max, linear_penalization in rows:
            try:
                best_param = json.loads(best_param_json)
            except json.JSONDecodeError:
                print(f"{run_id:<10} | Invalid JSON format")
                continue
            
            total_area = 0.0
            rect_area = x_max * y_max
            
            # Iterate over the parameters in chunks of 5 (x, y, A, B, C)
            for i in range(0, len(best_param), 5):
                # We need at least 5 elements for a complete ellipse
                if i + 4 < len(best_param):
                    A, B, C = best_param[i+2], best_param[i+3], best_param[i+4]
                    
                    # Determinant of the quadratic form matrix [[A, B], [B, C]] is A*C - B^2
                    det = (A * C) - (B ** 2)
                    
                    if det > 0:
                        area = np.pi / np.sqrt(det)
                        total_area += area
                    else:
                        print(f"Warning: Ellipse {i//5 + 1} in run {run_id} has a non-positive determinant.")
            
            unfilled_percent = 0.0
            if rect_area > 0:
                unfilled_percent = ((rect_area - total_area) / rect_area) * 100
            
            print(f"{run_id:<10} | {linear_penalization:<15.2f} | {total_area:<20.6f} | {f'{unfilled_percent:.2f}%':<15}")
            
    except sqlite3.OperationalError as e:
        print(f"Database error: {e}")
    finally:
        conn.close()

if __name__ == '__main__':
    # Example usage:
    db_file = 'jasminum-remote.db'
    #db_file = 'experiments.db'
    table = 'results'  # Adjust the table name to match your database
    #table = 'x1y1n2lambda16_SA'
    calculate_ellipse_areas_from_db(db_file, table)
