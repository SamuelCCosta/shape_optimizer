from square_solver import *
import numpy as np
import math
import random
from pebble import concurrent, ProcessExpired
from concurrent.futures import TimeoutError

h = 0.02
heat_source = 10.0
base_temp = 0.0
penalization = 0.0
num_ellipses = 4
geometric_info = {'x_max' : 1.0, 'y_max' : 1.0, 'MW_x' : 0.3, 'ME_x' : 0.7}
export_domain = False
export_solution = False


@concurrent.process(timeout = 2.0)
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
            
        print(ellipses.area(), 'pen = 30', 30.0*(1-ellipses.area()))
        return sqs.solve(ellipses)
    except ValueError: #configuração inválida
        return float('inf')
    except Exception as e:
        print(f'With params \n{params} \nError: {e}')
        return float('inf')

if __name__ == '__main__':
    params = [   0.53446306,    0.50259374,  243.77396884,  -33.96889496,  261.82802112,
    0.42235164,    0.70059529,  233.99220896,  -43.95489961,  288.09370986,
    0.56719944,    0.34218724,  289.65241533, -106.94255163,  203.10652993,
    0.46513083,    0.8759024,   208.72228596,   71.36566829,  130.41681181]

    print(f'Cost: {cost(params).result()}') # pyright: ignore[reportAttributeAccessIssue]