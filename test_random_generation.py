from square_solver import *
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse as MplEllipse
import numpy as np
import time
import multiprocessing as mp
from pebble import ProcessPool, ProcessExpired
from concurrent.futures import TimeoutError

def raw_solve_worker(geometric_info, h, heat_source, base_temp, num_ellipses, seed, penalization):
    start_time = time.perf_counter()
    bundle = EllipseBundle(geometric_info, h, num_ellipses)
    bundle.generate_random(seed=seed)
    elapsed_time = time.perf_counter() - start_time
    area = bundle.area()
    
    sqs = SquareSolver(geometric_info, h, heat_source, base_temp)
    try:
        raw_cost = sqs.solve(bundle)
        full_area = geometric_info['x_max'] * geometric_info['y_max']
        area_percent = (full_area - area) / full_area
        penalized_cost = raw_cost + penalization * area_percent
        return elapsed_time, area, False, penalized_cost
    except Exception:
        return elapsed_time, area, True, float('inf')

def test_random_generation(geometric_info, h, num_ellipses, seed):
    # 1. Create the bundle
    bundle = EllipseBundle(geometric_info, h, num_ellipses)
    
    # 2. Generate random ellipses with a specific seed for reproducibility
    print(f"Generating {num_ellipses} random ellipses with {seed=}...")
    start_time = time.perf_counter()
    bundle.generate_random(seed=seed)
    elapsed_time = time.perf_counter() - start_time
    
    # 3. Retrieve the generated ellipses from the bundle
    ellipses = bundle.bundle
    print(f"Successfully generated {len(ellipses)} ellipses in {elapsed_time:.6f} seconds.\n")
    
    # Setup plot
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(0.0, geometric_info['x_max'])
    ax.set_ylim(0.0, geometric_info['y_max'])
    ax.set_aspect('equal')
    ax.grid(True, linestyle='--', alpha=0.6)
    
    colors = ['red', 'blue', 'green', 'purple', 'orange', 'cyan']
    
    # 4. Iterate over the ellipses to print their parameters and plot them
    for i, e in enumerate(ellipses):
        xc, yc = e.center[0], e.center[1]
        A, B, C = e.quadratic_form[0,0], e.quadratic_form[0,1], e.quadratic_form[1,1]
        
        print(f"Ellipse {i+1}: x={xc:.4f}, y={yc:.4f}, A={A:.4f}, B={B:.4f}, C={C:.4f}")
        
        # Eigen decomposition to find axes and rotation for plotting
        M = np.array([[A, B], [B, C]])
        w, v = np.linalg.eigh(M)
        
        a = 1.0 / np.sqrt(w[0])
        b = 1.0 / np.sqrt(w[1])
        theta = np.degrees(np.arctan2(v[1, 0], v[0, 0]))
        
        ellipse_patch = MplEllipse(
            xy=(xc, yc), width=2*a, height=2*b, angle=theta,
            edgecolor=colors[i % len(colors)], facecolor='none', lw=2
        )
        ax.add_patch(ellipse_patch)
        ax.plot(xc, yc, 'k+', markersize=5)
    
    plt.title('Randomly Generated Ellipses via C++')
    plt.show()

def test_average_runtime(geometric_info, h, num_ellipses, seed, iteration_count, penalization):
    print(f"Testing average runtime with {num_ellipses} ellipses over {iteration_count} iterations...")
    total_time = 0.0
    bundle_area = 0.0
    error_count = 0
    timeout_count = 0
    crash_count = 0
    finite_costs = []
    
    spawn_ctx = mp.get_context('spawn')
    heat_source = 10.0
    base_temp = 0.0

    for i in range(iteration_count):
        current_seed = seed + i if seed != 0 else 0
        
        try:
            with ProcessPool(max_workers=1, context=spawn_ctx, max_tasks=1) as pool:
                # 5.0 second timeout, adjust if some valid meshes naturally take longer
                future = pool.schedule(raw_solve_worker, args=[geometric_info, h, heat_source, base_temp, num_ellipses, current_seed, penalization], timeout=5.0)
                elapsed_time, area, solver_error, cost = future.result()
                
                total_time += elapsed_time
                bundle_area += area
                if solver_error:
                    error_count += 1
                elif cost != float('inf'):
                    finite_costs.append(cost)
        except TimeoutError:
            timeout_count += 1
        except ProcessExpired:
            crash_count += 1

    # Calculate average only on runs that successfully returned time/area
    tracked_count = iteration_count - timeout_count - crash_count
    avg_time = total_time / tracked_count if tracked_count > 0 else 0.0
    avg_area = bundle_area / tracked_count if tracked_count > 0 else 0.0
    
    print(f"Average generation time: {avg_time:.6f} seconds.")
    print(f"Average area: {avg_area:.6f}")
    print(f"Solver Errors: {error_count} | Crashes: {crash_count} | Timeouts: {timeout_count} (out of {iteration_count})")
    print()
    return finite_costs

if __name__ == '__main__':
    geometric_info = {'x_max': 1.0, 'y_max': 1.0, 'MW_x': 0.3, 'ME_x': 0.7}
    h = 0.02
    num_ellipses = 2
    seed = 0 # 0 = NO SEED
    linear_pen = 128.0
    iteration_count = 500
    #test_random_generation(geometric_info, h, num_ellipses, seed)

    finite_costs = test_average_runtime(geometric_info, h, num_ellipses, seed, iteration_count, linear_pen)

    finite_costs_avg = sum(finite_costs) / len(finite_costs)
    finite_costs_var = sum((x - finite_costs_avg) ** 2 for x in finite_costs) / (len(finite_costs))
    finite_costs_std = finite_costs_var ** 0.5
    finite_costs_min = min(finite_costs)
    finite_costs_max = max(finite_costs)

    print(f'Avg cost: {finite_costs_avg}, std: {finite_costs_std}, min: {finite_costs_min}, max: {finite_costs_max}')

    if finite_costs:
        plt.figure(figsize=(8, 6))
        plt.hist(finite_costs, bins=30, color='skyblue', edgecolor='black', alpha=0.7)
        plt.title(f'Distribution of Randomly Generated Costs (n={len(finite_costs)})')
        plt.xlabel('Penalized Cost')
        plt.ylabel('Frequency')
        plt.grid(True, axis='y', linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.show()