from square_solver import *
import numpy as np
import math
import random
import signal

#multiprocessing.set_start_method('spawn')
#provavelmente comportamentos estranhos com fork (default do linux)

h = 0.02
heat_source = 10.0
base_temp = 0.0
penalization = 30.0
num_ellipses = 4

cfg = DomainConfig(1.0, 1.0, 0.3, 0.7, h, num_ellipses)

sqs = SquareSolver(cfg, heat_source, base_temp, penalization, False)

class TimeoutError(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutError("Timed out!")

def cost(params, sq_solver : SquareSolver):
    signal.signal(signal.SIGALRM, timeout_handler)
    ellipses = EllipseBundle(cfg)
    n_param = 5
    signal.alarm(2) # 2 seconds
    try:
        for i in range(num_ellipses):
            idx = i * n_param
            ellipses.add(Ellipse(params[idx], params[idx + 1], params[idx + 2], params[idx + 3], params[idx + 4]))
        result = sq_solver.solve(ellipses)
        signal.alarm(0)
        return result
    except ValueError: #configuração inválida
        signal.alarm(0)
        return float('inf')
    except TimeoutError: #timeout devido à malha frontal
        print(f'Timeout with parameters {params}')
        return float('inf')
    except: #unrelated (por exemplo FiniteElement a falhar)
        signal.alarm(0)
        print(f'Something went wrong with parameters {params}')
        return float('inf')


def simmulated_annealing(sq_solver : SquareSolver, initial_params, scales, initial_temp, cooling_rate, max_iter):
    # Initialization
    current_params = np.array(initial_params, dtype=float)
    current_cost = cost(initial_params, sq_solver)

    best_params = initial_params.copy()
    best_cost = current_cost

    temp = initial_temp

    print(f"Initial cost: {current_cost}")

    failed = 0 #track how many failed configurations happened
    for i in range(max_iter):
        if i % 10 == 0:
            print(f'i = {i}')

        noise = np.random.normal(0, 1.0, num_ellipses * 5) * scales
        neigh_params = current_params + noise
        #print(neigh_params.tolist())
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
        
        #cool down if config was valid
        if neigh_cost != float('inf'):
            temp *= cooling_rate

        if temp < 1e-8:
            i_count = i+1
            print(f'i = {i}, failed = {failed}, fail% = {failed/i_count * 100:.2f}%')
            break
    
    #print(f'i = {max_iter}, failed = {failed}, fail% = {failed/max_iter * 100:.2f}')
    return best_params


if __name__ == "__main__":
    initial_ellipses = [   0.53446306,    0.50259374,  243.77396884,  -33.96889496,  261.82802112,
    0.42235164,    0.70059529,  233.99220896,  -43.95489961,  288.09370986,
    0.56719944,    0.34218724,  289.65241533, -106.94255163,  203.10652993,
    0.46513083,    0.8759024,   208.72228596,   71.36566829,  130.41681181]

    scales = np.array([0.05, 0.05, 5.0, 5.0, 5.0] * num_ellipses, dtype=float)
    best_params = simmulated_annealing(sqs, initial_ellipses, scales, 250, 0.98, 5000)

    print('Best parameters:')
    print(best_params.tolist())
    print(f'Objective: {cost(best_params,sqs)}')

    #PROBLEMA: square_solver não quer fazer mais que uma resolução



    