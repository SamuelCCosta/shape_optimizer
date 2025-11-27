from square_solver import *
import numpy as np
import math
import random
import gc

#multiprocessing.set_start_method('spawn')
#provavelmente comportamentos estranhos com fork (default do linux)

h = 0.02
heat_source = 10.0
base_temp = 0.0
penalization = 10.0
num_ellipses = 4

cfg = DomainConfig(1.0, 1.0, 0.3, 0.7, h, num_ellipses)

def cost(params, sq_solver : SquareSolver):    
    ellipses = EllipseBundle(cfg)
    n_param = 5

    try:
        for i in range(num_ellipses):
            idx = i * n_param
            ellipses.add(Ellipse(params[idx], params[idx + 1], params[idx + 2], params[idx + 3], params[idx + 4]))
        
        return sq_solver.solve(ellipses)
        
    except ValueError:
        return float('inf')

def simmulated_annealing(sq_solver : SquareSolver, initial_params, scales, initial_temp = 100.0, cooling_rate = 0.90, max_iter = 1000):
    # Initialization
    current_params = np.array(initial_params, dtype=float)
    current_cost = cost(initial_params, sq_solver)

    best_params = initial_params.copy()
    best_cost = current_cost

    temp = initial_temp

    print(f"Initial cost: {current_cost}")

    failed = 0 #track how many failed configurations happened
    for i in range(max_iter):
        noise = np.random.normal(0, 1.0, num_ellipses * 5) * scales
        neigh_params = current_params + noise
        print(neigh_params)
        neigh_cost = cost(neigh_params, sq_solver)

        #definir probabilidades para o teste
        if neigh_cost == float('inf'): #configuração inválida
            prob = 0.0
            failed += 1
        else:
            delta = neigh_cost - current_cost
            if delta < 0: #improves
                prob = 1.0
            else:
                try: #se por acaso der overflow
                    prob = math.exp(-delta/temp)
                except OverflowError:
                    prob = 0.0
        
        # teste e verificar se encontrámos valor ótimo
        if random.random() < prob:
            current_params = neigh_params
            current_cost = neigh_cost

            if current_cost < best_cost:
                best_params = current_params.copy()
                best_cost = current_cost
                print(f'Best cost found: {best_cost:.6f} at temp = {temp:.2f}')
        
        #cool down
        temp *= cooling_rate
        if temp < 1e-8:
            i_count = i+1
            print(f'i = {i}, failed = {failed}, fail% = {failed/i_count * 100:.2f}%')
            break
    
    print(f'i = {max_iter}, failed = {failed}, fail% = {failed/max_iter * 100:.2f}')
    return best_params


if __name__ == "__main__":
    initial_ellipses = [ 
        4.29481844e-01, 6.11469557e-01, 1.19159636e+02,  5.89029810e-01, 8.84212799e+01,
        5.27868318e-01, 4.12821412e-01, 2.29649199e+02, -2.40292652e+01, 9.25764037e+01,
        4.23352361e-01, 8.60151240e-01, 8.40236125e+01, -2.27213249e+01, 1.16611590e+02,
        5.20663979e-01, 9.53702110e-02, 1.09547987e+02, -2.83772329e+01, 2.08775829e+02]

    sqs = SquareSolver(cfg, heat_source, base_temp, penalization, False)

    scales = np.array([0.05, 0.05, 5.0, 5.0, 5.0] * num_ellipses, dtype=float)
    best_params = simmulated_annealing(sqs, initial_ellipses, scales, 100, 0.95, 1000)
    
    #Só um SquareSolver pode estar ativo (por processo)
    del sqs
    gc.collect()

    print(best_params)
    sqs_print = SquareSolver(cfg, heat_source, base_temp, penalization, True)
    cost(best_params, sqs_print)
    



    