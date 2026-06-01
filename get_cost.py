from square_solver import *
import numpy as np
from pebble import concurrent
import time

#h = 0.02
heat_source = 10.0
base_temp = 0.0
num_ellipses = 0
geometric_info = {'x_max' : 1.0, 'y_max' : 1.0, 'MW_x' : 0.3, 'ME_x' : 0.7}
big_area = geometric_info['x_max'] * geometric_info['y_max']

penalization = 0.0
radius = 0.05
AC_param = (1 / radius) ** 2
params = [0.5, 0.5, AC_param, 0.0, AC_param]
print(f'{penalization=}')

#h_values = [0.02, 0.018, 0.01, 0.0075, 0.005]
h_values = [0.015]

params = [
  0.8618206759053694,
  0.638244200440503,
  225.67493248788125,
  -160.57422243346917,
  280.5987404955669,
  0.8379835207017033,
  0.3939957829236339,
  192.45436487876822,
  -21.333606792881376,
  78.24251380072735,
  0.15718368156532972,
  0.11260201636280878,
  398.10059552791836,
  -241.77576084751604,
  451.12919820660574,
  0.3854623436369899,
  0.4835325070889105,
  14.240792213872059,
  7.025599246051849,
  8.336150579822723
]
params=[]


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

def get_gradient(params, h, delta=1e-6):
    """Computes the numerical gradient of the penalized cost using forward differences."""
    gradient = []
    
    try:
        base_obj, base_area = cost(params, h).result() #type: ignore
        base_cost = base_obj + (base_area * penalization)
    except Exception:
        return [float('nan')] * len(params)
        
    for i in range(len(params)):
        perturbed_params = list(params)
        perturbed_params[i] += delta
        
        try:
            pert_obj, pert_area = cost(perturbed_params, h).result() #type: ignore
            if pert_obj == float('inf'):
                gradient.append(float('nan'))
            else:
                pert_cost = pert_obj + (pert_area * penalization)
                gradient.append((pert_cost - base_cost) / delta)
        except Exception:
            gradient.append(float('nan'))
            
    return gradient

def get_gradient_central(params, h, delta=1e-6):
    """Computes the numerical gradient of the penalized cost using central differences."""
    gradient = []
    
    for i in range(len(params)):
        params_plus = list(params)
        params_minus = list(params)
        params_plus[i] += delta
        params_minus[i] -= delta
        
        try:
            plus_obj, plus_area = cost(params_plus, h).result() #type: ignore
            minus_obj, minus_area = cost(params_minus, h).result() #type: ignore
            
            if plus_obj == float('inf') or minus_obj == float('inf'):
                gradient.append(float('nan'))
            else:
                plus_cost = plus_obj + (plus_area * penalization)
                minus_cost = minus_obj + (minus_area * penalization)
                gradient.append((plus_cost - minus_cost) / (2 * delta))
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
            calculate_gradient = False
            if calculate_gradient:
                grad = get_gradient(params, h_test)
                print(f"  -> Forward Gradient: {[round(g, 4) if not np.isnan(g) else 'NaN' for g in grad]}")
                grad_central = get_gradient_central(params, h_test)
                print(f"  -> Central Gradient: {[round(g, 4) if not np.isnan(g) else 'NaN' for g in grad_central]}")
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