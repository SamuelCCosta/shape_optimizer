from square_solver import *
import numpy as np
from pebble import concurrent
import time

#h = 0.02
heat_source = 10.0
base_temp = 0.0
num_ellipses = 2
geometric_info = {'x_max' : 1.0, 'y_max' : 1.0, 'MW_x' : 0.3, 'ME_x' : 0.7}
big_area = geometric_info['x_max'] * geometric_info['y_max']

penalization = 8.0
radius = 0.05
AC_param = (1 / radius) ** 2
params = [0.5, 0.5, AC_param, 0.0, AC_param]
print(f'{penalization=}')

#h_values = [0.02, 0.018, 0.01, 0.0075, 0.005]
h_values = [0.02, 0.01, 0.003]
params = [0.17192565444874897, 0.2848518988652481, 147.39626194489895, 25.286903310697127, 166.98236756042007, 0.5857396208311669, 0.5182635065978133, 7.566452987234112, 2.3577718294147756, 5.814000668355892]


@concurrent.process(timeout = 30.0)
def cost(params, h):
    sqs : SquareSolver = SquareSolver(geometric_info, h,
                                      heat_source, base_temp)
    ellipses = EllipseBundle(geometric_info, h, num_ellipses)

    n_param = 5

    try:
        for i in range(num_ellipses):
            idx = i * n_param
            ellipses.add(Ellipse(params[idx], params[idx + 1], params[idx + 2],
                                 params[idx + 3], params[idx + 4]))
        
        percent_area = (big_area - ellipses.area()) / big_area

        return sqs.solve(ellipses), percent_area
    except ValueError: #configuração inválida
        return float('inf'), float('inf')
    except Exception as e:
        print(f'With params \n{params} \nError: {e}')
        return float('inf'), float('inf')

def get_gradient(params, h, delta=1e-5):
    """Computes the numerical gradient of the penalized cost using forward differences."""
    gradient = []
    
    try:
        base_obj, base_area = cost(params, h).result()
        base_cost = base_obj + (base_area * penalization)
    except Exception:
        return [float('nan')] * len(params)
        
    for i in range(len(params)):
        perturbed_params = list(params)
        perturbed_params[i] += delta
        
        try:
            pert_obj, pert_area = cost(perturbed_params, h).result()
            if pert_obj == float('inf'):
                gradient.append(float('nan'))
            else:
                pert_cost = pert_obj + (pert_area * penalization)
                gradient.append((pert_cost - base_cost) / delta)
        except Exception:
            gradient.append(float('nan'))
            
    return gradient

if __name__ == '__main__':
    print(f"{'h':>8} | {'Time (s)':>10} | {'Cost':>12} | {'Penalized':>12}")
    print("-" * 55)

    for h_test in h_values:
        start_time = time.perf_counter()
        
        try:
            future = cost(params, h_test)
            objective, percent_area = future.result() # pyright: ignore
            
            elapsed = time.perf_counter() - start_time
            penalized_cost = objective + (percent_area * penalization)
            
            print(f"{h_test:8.3f} | {elapsed:10.4f} | {objective:12.4f} | {penalized_cost:12.4f}")
            
            grad = get_gradient(params, h_test)
            print(f"         -> Gradient: {[round(g, 4) if not np.isnan(g) else 'NaN' for g in grad]}")
        except Exception as e:
            elapsed = time.perf_counter() - start_time
            # Captura Crash de C++ (ProcessExpired, BrokenProcessPool, etc)
            print(f"{h_test:8.4f} | {elapsed:10.4f} | {'CRASH/ERR':>12} | {'N/A':>12}")

    
    
    
    #objective, percent_area = cost(params).result() # pyright: ignore[reportAttributeAccessIssue]


'''
    print(f'Percent of total area = {percent_area * 100 :.3f}%')
    print(f'Non-penalized cost: {objective:.4f}')
    print(f'Penalized (final) cost: {objective + percent_area * penalization:.4f}')
    '''