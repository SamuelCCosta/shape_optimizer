from square_solver import *
import numpy as np
import math
import random

h = 0.02
heat_source = 10.0
base_temp = 0.0
penalization = 10.0
export_mesh = True
num_ellipses = 4

cfg = DomainConfig(1.0, 1.0, 0.3, 0.7, h, num_ellipses)

sqs = SquareSolver(cfg, heat_source, base_temp, penalization, export_mesh)

def cost(params):    
    ellipses = EllipseBundle(cfg)
    n_param = 5

    try:
        for i in range(num_ellipses):
            id = i * n_param
            ellipses.add(Ellipse(params[id], params[id + 1], params[id + 2], params[id + 3], params[id + 4]))
        
        return sqs.solve(ellipses)
        
    except ValueError:
        return float('inf')


if __name__ == '__main__':
    local_optima = [   0.53446306,    0.50259374,  243.77396884,  -33.96889496,  261.82802112,
    0.42235164,    0.70059529,  233.99220896,  -43.95489961,  288.09370986,
    0.56719944,    0.34218724,  289.65241533, -106.94255163,  203.10652993,
    0.46513083,    0.8759024,   208.72228596,   71.36566829,  130.41681181]
    value = cost(local_optima)
    print(f'Cost: {value}')
