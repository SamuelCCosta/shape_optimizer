import sqlite3
import sys

def rename_scales_to_perturbation(db_path="experiments.db"):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    try:
        # Fetch all tables from the database
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = cursor.fetchall()
        
        for (table_name,) in tables:
            if table_name.startswith("sqlite_"):
                continue  # Skip internal sqlite tables
            
            # Check if 'scales' column exists in this table
            cursor.execute(f"PRAGMA table_info({table_name})")
            columns = [col[1] for col in cursor.fetchall()]
            
            if 'scales' in columns:
                print(f"Renaming 'scales' to 'perturbation' in table: '{table_name}'")
                cursor.execute(f"ALTER TABLE {table_name} RENAME COLUMN scales TO perturbation")
        
        conn.commit()
        print("Successfully updated database schema.")
    except sqlite3.Error as e:
        print(f"Database error occurred: {e}")
    finally:
        conn.close()

if __name__ == '__main__':
    db_file = sys.argv[1] if len(sys.argv) > 1 else 'experiments.db'
    rename_scales_to_perturbation(db_file)
