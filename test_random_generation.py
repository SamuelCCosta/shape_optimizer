from square_solver import *
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse as MplEllipse
import numpy as np
import time
import multiprocessing as mp
from pebble import ProcessPool, ProcessExpired
from concurrent.futures import TimeoutError

def raw_solve_worker(geometric_info, h, heat_source, base_temp, num_ellipses, seed):
    start_time = time.perf_counter()
    bundle = EllipseBundle(geometric_info, h, num_ellipses)
    bundle.generate_random(seed=seed)
    elapsed_time = time.perf_counter() - start_time
    area = bundle.area()
    
    sqs = SquareSolver(geometric_info, h, heat_source, base_temp)
    try:
        sqs.solve(bundle)
        return elapsed_time, area, False
    except Exception:
        return elapsed_time, area, True

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

def test_average_runtime(geometric_info, h, num_ellipses, seed, iteration_count):
    print(f"Testing average runtime with {num_ellipses} ellipses over {iteration_count} iterations...")
    total_time = 0.0
    bundle_area = 0.0
    error_count = 0
    timeout_count = 0
    crash_count = 0
    
    spawn_ctx = mp.get_context('spawn')
    heat_source = 10.0
    base_temp = 0.0

    for i in range(iteration_count):
        current_seed = seed + i if seed != 0 else 0
        
        try:
            with ProcessPool(max_workers=1, context=spawn_ctx, max_tasks=1) as pool:
                # 5.0 second timeout, adjust if some valid meshes naturally take longer
                future = pool.schedule(raw_solve_worker, args=[geometric_info, h, heat_source, base_temp, num_ellipses, current_seed], timeout=5.0)
                elapsed_time, area, solver_error = future.result()
                
                total_time += elapsed_time
                bundle_area += area
                if solver_error:
                    error_count += 1
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

if __name__ == '__main__':
    geometric_info = {'x_max': 1.0, 'y_max': 1.0, 'MW_x': 0.3, 'ME_x': 0.7}
    h = 0.02
    num_ellipses = 1
    seed = 0 # 0 = NO SEED
    iteration_count = 100
    '''nums_ellipses = [4,8,12]
    for n_ellipses in nums_ellipses:
        test_average_runtime(geometric_info, h, n_ellipses, seed, iteration_count)'''
    test_average_runtime(geometric_info, h, num_ellipses, seed, iteration_count)