from square_solver import *

# Constants
heat_source = 10.0
base_temp = 0.0
export_mesh = False

ellipses = EllipseBundle()

try:
    ellipses.add(Ellipse(0.5, 0.5, 78.1888, -37.7726, 34.5663)) # a = 0.28, b = 0.1, theta = pi/3
    ellipses.add(Ellipse(0.7, 0.31, 204.082, 0.0, 12.755)) # a = 0.28, b = 0.07, theta = pi/2
    ellipses.add(Ellipse(0.19, 0.59, 45.0817, -25.1834, 200.0929)) # a = 0.156, b = 0.07, theta = pi/20
    #ellipses.add(Ellipse(0.19, 0.17, 45.0817, -25.1834, 200.0929)) # a = 0.156, b = 0.07, theta = pi/20

    # Call the function directly
    value = objective(heat_source, base_temp, ellipses, export_mesh)
    #para já a função objetivo é o integral de linha de cima + 10 * volume
    print(f"Objective: {value}")

except ValueError as e:
    print(f"Error: {e}")