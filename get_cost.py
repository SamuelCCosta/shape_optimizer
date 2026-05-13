from square_solver import *
import numpy as np
from pebble import concurrent
import time

#h = 0.02
heat_source = 10.0
base_temp = 0.0
num_ellipses = 1
geometric_info = {'x_max' : 1.0, 'y_max' : 1.0, 'MW_x' : 0.3, 'ME_x' : 0.7}
big_area = geometric_info['x_max'] * geometric_info['y_max']

penalization = 5.0
radius = 0.05
AC_param = (1 / radius) ** 2
params = [0.5, 0.5, AC_param, 0.0, AC_param]
print(f'{penalization=}')

#h_values = [0.02, 0.018, 0.01, 0.0075, 0.005]
h_values = [0.02]


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