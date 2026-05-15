from square_solver import *
import numpy as np
from pebble import concurrent
import time

@concurrent.process(timeout = 30.0)
def cost_function(sqs_params, ellipse_bundle_params, linear_penalization, state) -> float:
        '''Get full cost, including (linear) penalization.'''
        num_ellipses = ellipse_bundle_params['num_ellipses']
        full_area = sqs_params["geometric_config"]['x_max'] * sqs_params["geometric_config"]['y_max']
        sqs = SquareSolver(**sqs_params)
        ellipses = EllipseBundle(**ellipse_bundle_params)
        n_param = 5 #ellipse parameter number
        try:
            for i in range(num_ellipses):
                idx = i * n_param
                ellipses.add(Ellipse(state[idx], state[idx + 1], state[idx + 2],
                                    state[idx + 3], state[idx + 4]))
            area_percent = (full_area - ellipses.area()) / full_area
            return sqs.solve(ellipses) + linear_penalization * area_percent
        except: #invalid configuration
            return float('inf')


if __name__ == '__main__':
    geometric_params = {
        "geometric_config" : {'x_max': 1.0, 'y_max': 1.0, 'MW_x': 0.3, 'ME_x': 0.7},
        "h" : 0.02,
        "heat_sources" : 10.0,
        "base_temp" : 0.0,
        "num_ellipses" : 2
    }
    
    linear_pen = 8.0

    sqs_params = {k: v for k, v in geometric_params.items() if k not in ("num_ellipses")}
    ellipse_bundle_params = {k: geometric_params[k] for k in ["geometric_config", "h", "num_ellipses"]}

    state = [0.17192565444874897, 0.2848518988652481, 147.39626194489895, 25.286903310697127, 166.98236756042007, 0.5857396208311669, 0.5182635065978133, 7.566452987234112, 2.3577718294147756, 5.814000668355892]


    print(cost_function(sqs_params, ellipse_bundle_params, linear_pen, state).result())