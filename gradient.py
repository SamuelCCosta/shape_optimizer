"""Adaptive finite-difference gradient descent. Edit the settings below and run python gradient.py."""

from math import isfinite, sqrt


def gradient_descent(cost, initial_params, learning_rate=0.01, delta=1e-5,
                     max_iterations=100, tolerance=1e-6, max_backtracks=12,
                     growth_factor=1.2, min_learning_rate=1e-10,
                     max_learning_rate=1.0):
    """Minimize a scalar cost; invalid configurations should return infinity.

    Halve unsuccessful steps; start the next iteration at the accepted step
    times growth_factor, within the learning-rate bounds. A small sufficient
    decrease (Armijo condition) is required to accept a step.
    """
    if not isfinite(delta) or delta <= 0:
        raise ValueError('delta must be finite and positive')
    if (not all(isfinite(v) for v in (min_learning_rate, learning_rate, max_learning_rate))
            or not 0 < min_learning_rate <= learning_rate <= max_learning_rate):
        raise ValueError('Require 0 < min_learning_rate <= learning_rate <= max_learning_rate')
    if not isfinite(growth_factor) or growth_factor < 1:
        raise ValueError('growth_factor must be finite and at least 1')
    if max_backtracks < 0:
        raise ValueError('max_backtracks must be nonnegative')

    params = list(initial_params)
    current_cost = cost(params)
    if not isfinite(current_cost):
        raise ValueError('The initial parameters must have a finite cost')
    print(f'Initial cost: {current_cost:.10g}', flush=True)

    for iteration in range(1, max_iterations + 1):
        gradient = []
        for i in range(len(params)):
            perturbed = params.copy()
            perturbed[i] += delta
            perturbed_cost = cost(perturbed)
            if isfinite(perturbed_cost):
                derivative = (perturbed_cost - current_cost) / delta
            else:
                # Use a backward difference if the forward probe is invalid.
                perturbed[i] = params[i] - delta
                perturbed_cost = cost(perturbed)
                derivative = (current_cost - perturbed_cost) / delta
            if not isfinite(derivative):
                print(f'Stopped: cannot estimate derivative for parameter {i}.')
                return params, current_cost
            gradient.append(derivative)

        gradient_squared = sum(g * g for g in gradient)
        gradient_norm = sqrt(gradient_squared)
        if gradient_norm <= tolerance:
            print('Stopped: gradient norm is below tolerance.')
            break

        step = learning_rate
        accepted = False
        for backtracks in range(max_backtracks + 1):
            candidate = [p - step * g for p, g in zip(params, gradient)]
            candidate_cost = cost(candidate)
            required_decrease = 1e-4 * step * gradient_squared
            if (isfinite(candidate_cost) and candidate_cost < current_cost
                    and current_cost - candidate_cost >= required_decrease):
                accepted = True
                break
            if step <= min_learning_rate:
                break
            step = max(min_learning_rate, step * 0.5)
        if not accepted:
            print('Stopped: no improving step found; try a different delta or learning_rate.')
            break

        improvement = current_cost - candidate_cost
        params, current_cost = candidate, candidate_cost
        learning_rate = min(max_learning_rate, step * growth_factor)
        print(f'{iteration:4d}: cost={current_cost:.10g}, '
              f'gradient_norm={gradient_norm:.6g}, step={step:.6g}, '
              f'backtracks={backtracks}, next_step={learning_rate:.6g}', flush=True)
        if improvement <= tolerance:
            print('Stopped: cost improvement is below tolerance.')
            break
    else:
        print('Stopped: reached max_iterations.')

    return params, current_cost


if __name__ == '__main__':
    from new_get_cost import cost_function

    # Same flattened format as get_cost.py: [x, y, A, B, C] per ellipse.
    params = [0.3392969344746265,0.13954606163739616,442.92542158852655,25.590001936055845,110.86831946689942,
  0.5012864378485594,0.0985756063722968,431.8846232657603,116.66065713084939,349.85687439958383,
  0.4827787298211445,0.5865753255297326,32.81204342585769,-0.35716103899947116,6.4908184623175424,
  0.660961903812361,0.09596265214670122,275.8120498248556,107.60075807191163,228.6757397684626
]
    geometric_info = {'x_max': 1.0, 'y_max': 1.0, 'MW_x': 0.3, 'ME_x': 0.7}
    h = 0.02
    heat_source = 10.0
    base_temp = 0.0
    penalization = 8.0

    learning_rate = 5e-4  # Initial trial step; adapts after each accepted update.
    max_learning_rate = 1.0  # Upper bound on the adaptive step.
    min_learning_rate = 1e-10
    growth_factor = 1.2  # Try a 20% larger step after success; halve on failure.
    delta = 1e-4  # Finite-difference step; very small steps may amplify mesh noise.
    max_iterations = 100
    tolerance = 1e-6

    if not params or len(params) % 5:
        raise ValueError('Provide five parameters (x, y, A, B, C) per ellipse')
    sqs_params = dict(geometric_config=geometric_info, h=h,
                      heat_sources=heat_source, base_temp=base_temp)
    bundle_params = dict(geometric_config=geometric_info, h=h,
                         num_ellipses=len(params) // 5)

    def cost(state):
        try:
            return cost_function(sqs_params, bundle_params, penalization, state).result()
        except Exception:
            # A solver timeout/crash is treated like an invalid configuration.
            return float('inf')

    best_params, best_cost = gradient_descent(
        cost, params, learning_rate=learning_rate, delta=delta,
        max_iterations=max_iterations, tolerance=tolerance,
        growth_factor=growth_factor, min_learning_rate=min_learning_rate,
        max_learning_rate=max_learning_rate,
    )
    print(f'Best cost: {best_cost:.10g}')
    print(f'Best params: {best_params}')
