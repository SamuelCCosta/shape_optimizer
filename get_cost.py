from square_solver import *
import numpy as np
import math
import random
from pebble import concurrent, ProcessExpired
from concurrent.futures import TimeoutError

h = 0.02
heat_source = 10.0
base_temp = 0.0
penalization = 50.0
num_ellipses = 1
geometric_info = {'x_max' : 1.0, 'y_max' : 1.0, 'MW_x' : 0.3, 'ME_x' : 0.7}
export_domain = False
export_solution = False


@concurrent.process(timeout = 2)
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

if __name__ == '__main__':
    params = [0.4528497255555963, 0.5195018072130883, 255.24891698211133, -38.35196654400833, 262.5530482957416]

    print(f'Cost: {cost(params).result()}') # pyright: ignore[reportAttributeAccessIssue]