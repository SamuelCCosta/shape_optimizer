from square_solver import *
import numpy as np
import math
import random
import multiprocessing as mp
from pebble import concurrent, ProcessExpired
from concurrent.futures import TimeoutError

spawn_ctx = mp.get_context('spawn')

h = 0.02
heat_source = 10.0
base_temp = 0.0
penalization = 50.0
num_ellipses = 1
geometric_info = {'x_max' : 1.0, 'y_max' : 1.0, 'MW_x' : 0.3, 'ME_x' : 0.7}
export_domain = False
export_solution = False

@concurrent.process(timeout = 2.0, mp_context = spawn_ctx)
def cost(params):
    sqs : SquareSolver = SquareSolver(geometric_info, h,
                                      heat_source, base_temp,
                                      penalization,
                                      export_domain, export_solution)
    ellipses = EllipseBundle(geometric_info, h, num_ellipses)
    n_param = 5 #each ellipse needs 5 parameters
    try:
        for i in range(num_ellipses):
            idx = i * n_param
            ellipses.add(Ellipse(params[idx], params[idx + 1], params[idx + 2],
                                 params[idx + 3], params[idx + 4]))
        return sqs.solve(ellipses)
    except ValueError: #configuração inválida
        return float('inf')
    except Exception as e:
        print(f'With params \n{params} \nError: {e}')
        return float('inf')
    
def simmulated_annealing(initial_params : list[float],
                         scales, initial_temp, cooling_rate,
                         max_iter):
    current_params = np.array(initial_params.copy())
    current_cost : float = cost(current_params).result() # pyright: ignore[reportAttributeAccessIssue]

    best_params = current_params.copy()
    best_cost = current_cost

    temp = initial_temp

    i = 0
    for i in range(max_iter):
        if i % 10 == 0:
            print(i)
        neigh_cost = float('inf')
        neighbour_params = np.zeros(num_ellipses * 5, dtype= float)
        #find an admissible solution
        while(neigh_cost == float('inf')):
            noise = np.random.normal(0.0, 1.0, num_ellipses * 5) * scales
            neighbour_params = current_params + noise
            future_neigh = cost(neighbour_params)
            try:
                neigh_cost = future_neigh.result() # pyright: ignore[reportAttributeAccessIssue]
            except TimeoutError:
                print('Took too long')
                neigh_cost = float('inf')
            except ProcessExpired as e:
                print('Worker crash')
                neigh_cost = float('inf')
        
        delta = neigh_cost - current_cost
        if delta < 0:
            prob = 1.0
        else:
            try: #se por acaso der overflow
                prob = math.exp(-delta/temp)
            except OverflowError:
                prob = 0.0
        
        # annealing test
        if random.random() < prob:
            current_params = neighbour_params
            current_cost = neigh_cost

            # check if is a minimum
            if current_cost < best_cost:
                best_params = current_params.copy()
                best_cost = current_cost
        
        temp *= cooling_rate

        if temp < 1e-08:
            print('temp zero reached')
            break
    
    if i == (max_iter - 1) :
        print(f'max_iter {max_iter} reached, temperature = {temp}')

    return best_params

if __name__ == '__main__':
    starting_parameters = [0.4528497255555963, 0.5195018072130883, 255.24891698211133, -38.35196654400833, 262.5530482957416]


    initial_temperature = 100
    cooling_rate = 0.98
    max_iterations = 1250
    scales = np.array([0.1, 0.1, 5.0, 25.0, 5.0] * num_ellipses, dtype=float)

    best_parameters = simmulated_annealing(starting_parameters, scales,
                                           initial_temperature,
                                           cooling_rate,
                                           max_iterations)
    
    print(best_parameters)
    print(f'Best cost: {cost(best_parameters).result()}') # pyright: ignore[reportAttributeAccessIssue]