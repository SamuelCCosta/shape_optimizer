from square_solver import *
import numpy as np
import math
import random
from pebble import concurrent, ProcessExpired
from concurrent.futures import TimeoutError

h = 0.02
heat_source = 10.0
base_temp = 0.0
penalization = 100.0
num_ellipses = 4
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
    params = [0.4528497255555963, 0.5195018072130883, 255.24891698211133, -38.35196654400833, 262.5530482957416,
              0.45017099450198184, 0.6923857353204883, 231.7094153654055, -47.15190162486053, 287.10515606531953,
              0.5265268828586713, 0.33421315950352537, 293.72495839816185, -112.08640183253189, 204.51894409046088,
              0.4319381994260433, 0.8739065566217921, 209.82444979164484, 67.87935546812955, 128.8843669760333]

    print(f'Cost: {cost(params).result()}')