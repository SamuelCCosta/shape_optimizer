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
    local_optima = [ 0.429481844, 0.611469557, 119.159636, 0.58902981, 88.4212799, 
    0.527868318, 0.412821412, 229.649199, -24.0292652, 92.5764037, 
    0.423352361, 0.86015124, 84.0236125, -22.7213249, 116.61159, 
    0.520663979, 0.095370211, 109.547987, -28.3772329, 208.775829]
    
    value = cost(local_optima)
    print(f'Cost: {value}')
