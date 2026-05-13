import sqlite3
import sys

def rename_table(db_path, old_table_name, new_table_name):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    try:
        print(f"Renaming table '{old_table_name}' to '{new_table_name}' in database '{db_path}'...")
        cursor.execute(f"ALTER TABLE {old_table_name} RENAME TO {new_table_name}")
        conn.commit()
        print("Successfully renamed the table.")
    except sqlite3.Error as e:
        print(f"Database error occurred: {e}")
    finally:
        conn.close()

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python rename_db_table.py <old_table_name> <new_table_name> [db_path]")
        sys.exit(1)
        
    old_name = sys.argv[1]
    new_name = sys.argv[2]
    db_file = sys.argv[3] if len(sys.argv) > 3 else 'experiments.db'
    
    rename_table(db_file, old_name, new_name)
