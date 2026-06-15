import time
import multiprocessing as mp
from pebble import ProcessPool, ProcessExpired
from concurrent.futures import TimeoutError
from square_solver import EllipseBundle, SquareSolver
import matplotlib.pyplot as plt
import numpy as np

def raw_solve_worker(geometric_info, h, heat_source, base_temp, num_ellipses, seed):
    start_time = time.perf_counter()
    
    # 1. Generation
    bundle = EllipseBundle(geometric_info, h, num_ellipses)
    bundle.generate_random(seed=seed)
    
    # 2. Early Failure Check (unable to generate non-overlapping ellipses)
    if len(bundle.bundle) < num_ellipses:
        elapsed_time = time.perf_counter() - start_time
        return 'early_fail', elapsed_time, None, None
        
    sqs = SquareSolver(geometric_info, h, heat_source, base_temp)
    
    # 3. Solver Execution
    try:
        cost = sqs.solve(bundle)
        area_filled = (geometric_info['x_max'] * geometric_info['y_max']) - bundle.area()
        elapsed_time = time.perf_counter() - start_time
        return 'suitable', elapsed_time, cost, area_filled
    except Exception:
        elapsed_time = time.perf_counter() - start_time
        return 'solver_fail', elapsed_time, None, None

def test_generation_runtimes(N, geometric_info, h, num_ellipses, seed=1000, timeout=5.0):
    spawn_ctx = mp.get_context('spawn')
    heat_source = 10.0
    base_temp = 0.0

    stats = {
        'early_fail': {'count': 0, 'total_time': 0.0},
        'suitable':   {'count': 0, 'total_time': 0.0},
        'solver_fail':{'count': 0, 'total_time': 0.0},
        'timeout':    {'count': 0, 'total_time': 0.0},
        'crash':      {'count': 0, 'total_time': 0.0}
    }
    computed_costs = []
    computed_areas = []

    print(f"Testing {N} random configurations...")
    
    def print_failure_info(fail_type, current_seed, iteration, extra_msg=""):
        msg = f"\n[!] {fail_type} detected on iteration {iteration} with seed {current_seed}"
        if extra_msg:
            msg += f" ({extra_msg})"
        print(msg)
        fail_bundle = EllipseBundle(geometric_info, h, num_ellipses)
        fail_bundle.generate_random(seed=current_seed)
        params = []
        for e in fail_bundle.bundle:
            params.extend([e.center[0], e.center[1], e.quadratic_form[0,0], e.quadratic_form[0,1], e.quadratic_form[1,1]])
        print(f"Failed Parameters: {params}\n")

    for i in range(N):
        current_seed = seed + i
        start_time_pool = time.perf_counter()
        
        try:
            with ProcessPool(max_workers=1, context=spawn_ctx) as pool:
                future = pool.schedule(
                    raw_solve_worker,
                    args=[geometric_info, h, heat_source, base_temp, num_ellipses, current_seed],
                    timeout=timeout
                )
                result_type, elapsed_time, cost, area_filled = future.result()
            
            elapsed_pool = time.perf_counter() - start_time_pool
            stats[result_type]['count'] += 1
            stats[result_type]['total_time'] += elapsed_pool
            if cost is not None:
                computed_costs.append(cost)
                computed_areas.append(area_filled)
            elif result_type != 'suitable':
                print_failure_info(result_type, current_seed, i + 1)
                
        except TimeoutError:
            elapsed_pool = time.perf_counter() - start_time_pool
            stats['timeout']['count'] += 1
            stats['timeout']['total_time'] += elapsed_pool
            print_failure_info('Timeout', current_seed, i + 1)
        except ProcessExpired:
            elapsed_pool = time.perf_counter() - start_time_pool
            stats['crash']['count'] += 1
            stats['crash']['total_time'] += elapsed_pool
            print_failure_info('Crash', current_seed, i + 1)
        except Exception as e:
            elapsed_pool = time.perf_counter() - start_time_pool
            stats['solver_fail']['count'] += 1
            stats['solver_fail']['total_time'] += elapsed_pool
            print_failure_info('Exception', current_seed, i + 1, str(e))

        if (i + 1) % 10 == 0:
            print(f"Processed {i + 1}/{N} configurations...")

    print("\n--- Runtime Statistics ---")
    for category, data in stats.items():
        count = data['count']
        if count > 0:
            avg_time = data['total_time'] / count
            print(f"{category:<12}: {count:>5} cases, Avg Runtime: {avg_time * 1000:.2f} ms")
        else:
            print(f"{category:<12}: {count:>5} cases, Avg Runtime: N/A")

    if computed_costs:
        avg_cost = sum(computed_costs) / len(computed_costs)
        print(f"\n--- Cost Statistics ---")
        print(f"Count: {len(computed_costs)} | Min: {min(computed_costs):.4f} | Max: {max(computed_costs):.4f} | Avg: {avg_cost:.4f}")

        plt.figure(figsize=(8, 6))
        plt.scatter(computed_areas, computed_costs, alpha=0.7, edgecolors='k', label='Data points')
        
        if len(computed_areas) > 1:
            m, b = np.polyfit(computed_areas, computed_costs, 1)
            x_vals = np.array([min(computed_areas), max(computed_areas)])
            plt.plot(x_vals, m * x_vals + b, color='red', linewidth=2, label=f'Linear fit: y = {m:.4f}x + {b:.4f}')
            
        plt.title('Objective Value vs. Area Filled')
        plt.xlabel('Area Filled (Square Area - Ellipses Area)')
        plt.ylabel('Objective Value (Cost)')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()
        plt.tight_layout()
        plt.show()

    return computed_costs, computed_areas

if __name__ == '__main__':
    geometric_info = {'x_max': 1.0, 'y_max': 1.0, 'MW_x': 0.3, 'ME_x': 0.7}
    h = 0.0125
    num_ellipses = 4
    N = 250
    seed = 123654
    test_generation_runtimes(N, geometric_info, h, num_ellipses, seed=seed)