import sqlite3
import sys

def delete_table(db_path, table_name):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    try:
        print(f"Deleting table '{table_name}' from database '{db_path}'...")
        cursor.execute(f"DROP TABLE IF EXISTS {table_name}")
        conn.commit()
        print("Successfully deleted the table (if it existed).")
    except sqlite3.Error as e:
        print(f"Database error occurred: {e}")
    finally:
        conn.close()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python delete_db_table.py <table_name> [db_path]")
        sys.exit(1)
        
    target_table = sys.argv[1]
    db_file = sys.argv[2] if len(sys.argv) > 2 else 'experiments.db'
    
    delete_table(db_file, target_table)