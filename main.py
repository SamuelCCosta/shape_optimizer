from square_solver import *
import numpy as np
import math
import random

#multiprocessing.set_start_method('spawn')
#provavelmente comportamentos estranhos com fork (default do linux)

h = 0.02
heat_source = 10.0
base_temp = 0.0
penalization = 100.0
num_ellipses = 4
geometric_info = {'x_max' : 1.0, 'y_max' : 1.0, 'MW_x' : 0.3, 'ME_x' : 0.7}
export_domain = False
export_solution = False

def main():
    sqs = SquareSolver(geometric_info, h, heat_source, base_temp, penalization, export_domain, export_solution)

    bundle = EllipseBundle(geometric_info, h, num_ellipses)
    bundle.add(Ellipse(0.5, 0.5, 13, 0 ,13))

    print(sqs.solve(bundle))

    return 0

if __name__ == '__main__':
    main()