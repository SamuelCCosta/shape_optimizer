import csv
import sys

def analyze_results(csv_filename):
    total_runs, frontal_fails, delaunay_fails = 0, 0, 0
    frontal_times, delaunay_times, cost_diffs = [], [], []
    
    with open(csv_filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            total_runs += 1
            c_front = float(row['cost_frontal'])
            c_del = float(row['cost_delaunay'])
            t_front = float(row['runtime_frontal'])
            t_del = float(row['runtime_delaunay'])
            
            if t_front == 0.0 or c_front == 0.0:
                frontal_fails += 1
            else:
                frontal_times.append(t_front)
                
            if t_del == 0.0 or c_del == 0.0:
                delaunay_fails += 1
            else:
                delaunay_times.append(t_del)
                
            if t_front > 0.0 and t_del > 0.0:
                cost_diffs.append(abs(c_front - c_del))
                
    print(f"Total runs evaluated: {total_runs}\n" + "-" * 40)
    print(f"Frontal Fails:  {frontal_fails} ({frontal_fails/total_runs*100:.2f}%)")
    print(f"Delaunay Fails: {delaunay_fails} ({delaunay_fails/total_runs*100:.2f}%)\n" + "-" * 40)
    if frontal_times: print(f"Avg Runtime (Frontal):  {sum(frontal_times)/len(frontal_times):.2f} ms")
    if delaunay_times: print(f"Avg Runtime (Delaunay): {sum(delaunay_times)/len(delaunay_times):.2f} ms")
    print("-" * 40)
    if cost_diffs: print(f"Avg Absolute Cost Difference: {sum(cost_diffs)/len(cost_diffs):.4f}")

if __name__ == '__main__':
    analyze_results(sys.argv[1] if len(sys.argv) > 1 else 'test_new_domain.csv')
