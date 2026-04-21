import sqlite3

conn = sqlite3.connect('experiments.db')
cursor = conn.cursor()

cursor.execute("SELECT * FROM results WHERE initial_temp = 100")
rows = cursor.fetchall()

for row in rows:
    print(row)

cursor.execute("DELETE FROM results WHERE initial_temp = 100")
conn.commit() 

conn.close()