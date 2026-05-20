from square_solver import *
import numpy as np
from pebble import concurrent
import time

#h = 0.02
heat_source = 10.0
base_temp = 0.0
num_ellipses = 4
geometric_info = {'x_max' : 1.0, 'y_max' : 1.0, 'MW_x' : 0.3, 'ME_x' : 0.7}
big_area = geometric_info['x_max'] * geometric_info['y_max']

penalization = 8.0
radius = 0.05
AC_param = (1 / radius) ** 2
params = [0.5, 0.5, AC_param, 0.0, AC_param]
print(f'{penalization=}')

#h_values = [0.02, 0.018, 0.01, 0.0075, 0.005]
h_values = [0.02]

params = [
  0.7140744550877801,
  0.8084752847524989,
  323.910908739876,
  -64.06973330334914,
  358.9984573566974,
  0.8770311056021156,
  0.20938680977638013,
  297.5637774573435,
  -22.75745330312562,
  121.67303229256659,
  0.8700388561011926,
  0.6911367346200966,
  515.9118472380815,
  130.6976341989308,
  325.5308045524322,
  0.33735004032166704,
  0.4867808895759328,
  15.266611718741553,
  4.793988330098576,
  6.174042863735902
]


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