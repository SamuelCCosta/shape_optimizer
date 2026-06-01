import math

def evaluate_objective(x_max, y_max, MW_x, ME_x, max_terms=1000):
    # Constant term contribution
    source_length = MW_x + (x_max - ME_x)
    j_const = (10.0 * y_max / x_max) * (source_length ** 2)
    
    j_series = 0.0
    for n in range(1, max_terms + 1):
        n_pi = n * math.pi
        
        arg_w = n_pi * MW_x / x_max
        arg_e = n_pi * ME_x / x_max
        arg_y = n_pi * y_max / x_max
        
        tanh_term = math.tanh(arg_y)
        sin_diff = math.sin(arg_w) - math.sin(arg_e)
        
        term = (20.0 * (x_max ** 2) * tanh_term / (n_pi ** 3)) * (sin_diff ** 2)
        j_series += term
        
            
    return j_const + j_series

if __name__ == '__main__':
    x_max = 1.0
    y_max = 1.0
    MW_x = 0.3
    ME_x = 0.7
    max_terms = 10000

    result = evaluate_objective(x_max, y_max, MW_x, ME_x, max_terms)
    
    print(f"Objective with no ellipses: {result}")

"""
THEORETICAL FOUNDATIONS & TERM EXPLANATIONS:

Separation of variables on [0, x_max] x [0, y_{max}] yields the solution:
u(x,y) = A0*y + sum_{n=1}^inf A_n * sinh(n*pi*y/x_max) * cos(n*pi*x/x_max)

Applying the boundary condition u_y(x, y_max) = g(x) via orthogonality, and 
integrating u(x, y_max) over Gamma_sources yields J(u).

CODE TERMS:
- j_const: Represents the n=0 background mode. It integrates the constant 
  temperature profile driven by the average heat flux over the total source length.
  Formula: (10 * y_max / x_max) * (MW_x + x_max - ME_x)^2

- sin_diff: [sin(n*pi*MW_x/x_max) - sin(n*pi*ME_x/x_max)]. Arises from 
  integrating the spatial harmonic cos(n*pi*x/x_max) over the source intervals.

- tanh_term: tanh(n*pi*y_max/x_max). Represents the aspect ratio effect, 
  originating from the ratio of the field potential (sinh) to its normal flux 
  derivative (cosh) evaluated at the upper boundary.

- term & j_series: 'term' computes the isolated quantitative weight of the 
  n-th Fourier mode to the objective functional. 'j_series' accumulates them.
  Contributions decay at a rate proportional to 1/n^3.
"""