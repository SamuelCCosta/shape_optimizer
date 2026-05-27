import csv
import json
import time
import multiprocessing as mp
from pebble import ProcessPool
from square_solver import SquareSolver, EllipseBundle, Ellipse

def solve_worker(method, params, geometric_info, h, heat_source, base_temp, num_ellipses):
    """
    Worker function to build the domain and run the requested solver.
    Runs in an isolated process to protect against C++ aborts.
    """
    bundle = EllipseBundle(geometric_info, h, num_ellipses)
    n_param = 5
    for i in range(num_ellipses):
        idx = i * n_param
        bundle.add(Ellipse(params[idx], params[idx+1], params[idx+2], params[idx+3], params[idx+4]))
    
    sqs = SquareSolver(geometric_info, h, heat_source, base_temp)
    
    start_time = time.perf_counter()
    
    if method == 'frontal':
        cost = sqs.solve_frontal(bundle)
    else:
        cost = sqs.solve(bundle)
        
    elapsed_ms = (time.perf_counter() - start_time) * 1000.0
    
    return cost, elapsed_ms

def test_new_domain(N, num_ellipses, csv_filename):
    geometric_info = {'x_max': 1.0, 'y_max': 1.0, 'MW_x': 0.3, 'ME_x': 0.7}
    h = 0.015
    heat_source = 10.0
    base_temp = 0.0
    timeout_seconds = 5.0 # Max time allowed for a single solver run
    
    spawn_ctx = mp.get_context('spawn')
    
    print(f"Starting test for {N} configurations...")
    print(f"Results will be saved to '{csv_filename}'.\n")
    
    with open(csv_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['cost_frontal', 'cost_delaunay', 'runtime_frontal', 'runtime_delaunay', 'params'])
        
        configs_tested = 0
        seed = 1000 # Starting seed for reproducibility
        
        while configs_tested < N:
            # 1. Generate random configuration
            bundle = EllipseBundle(geometric_info, h, num_ellipses)
            bundle.generate_random(seed=seed)
            seed += 1
            
            if len(bundle.bundle) < num_ellipses:
                continue # Skip if couldn't generate a valid, non-overlapping bundle
            
            params = []
            for e in bundle.bundle:
                params.extend([e.center[0], e.center[1], e.quadratic_form[0,0], e.quadratic_form[0,1], e.quadratic_form[1,1]])
            
            # 2. Evaluate frontal solver
            cost_frontal, runtime_frontal = 0.0, 0.0
            try:
                with ProcessPool(max_workers=1, context=spawn_ctx) as pool:
                    future = pool.schedule(solve_worker, args=['frontal', params, geometric_info, h, heat_source, base_temp, num_ellipses], timeout=timeout_seconds)
                    cost_frontal, runtime_frontal = future.result()
            except Exception:
                cost_frontal, runtime_frontal = 0.0, 0.0
            
            # 3. Evaluate delaunay solver
            cost_delaunay, runtime_delaunay = 0.0, 0.0
            try:
                with ProcessPool(max_workers=1, context=spawn_ctx) as pool:
                    future = pool.schedule(solve_worker, args=['delaunay', params, geometric_info, h, heat_source, base_temp, num_ellipses], timeout=timeout_seconds)
                    cost_delaunay, runtime_delaunay = future.result()
            except Exception:
                cost_delaunay, runtime_delaunay = 0.0, 0.0
                
            # 4. Save results
            writer.writerow([cost_frontal, cost_delaunay, runtime_frontal, runtime_delaunay, json.dumps(params)])
            f.flush() # ensure it writes to disk immediately
            
            configs_tested += 1
            print(f"[{configs_tested}/{N}] Frontal: (Cost: {cost_frontal:8.4f}, Time: {runtime_frontal:7.2f}ms) | Delaunay: (Cost: {cost_delaunay:8.4f}, Time: {runtime_delaunay:7.2f}ms)")

if __name__ == '__main__':
    # Customize N and the number of ellipses here
    test_new_domain(N=100, num_ellipses=4, csv_filename='test_new_domain.csv')

'''
h = 0.02
Total runs evaluated: 100
----------------------------------------
Frontal Fails:  63 (63.00%)
Delaunay Fails: 0 (0.00%)
----------------------------------------
Avg Runtime (Frontal):  116.12 ms
Avg Runtime (Delaunay): 30.67 ms
----------------------------------------
Avg Absolute Cost Difference: 0.0592

h = 0.015
Total runs evaluated: 100
----------------------------------------
Frontal Fails:  77 (77.00%)
Delaunay Fails: 0 (0.00%)
----------------------------------------
Avg Runtime (Frontal):  206.68 ms
Avg Runtime (Delaunay): 53.61 ms
----------------------------------------
Avg Absolute Cost Difference: 0.1402
'''
