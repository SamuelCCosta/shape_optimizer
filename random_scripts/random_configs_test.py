import numpy as np

def generate_configuration(x_max, y_max, h, num_ellipses):
    rng = np.random.default_rng()
    min_axis_length = 0.01
    y_min, y_max = h + min_axis_length, y_max - h - min_axis_length
    x_min, x_max = h + min_axis_length, x_max - h - min_axis_length
    
    output = []
    for _ in range(num_ellipses):
        xc = rng.uniform(x_min, x_max)
        yc = rng.uniform(y_min, y_max)
        dist_to_barrier = np.min([np.abs(xc - x_min), np.abs(xc - x_max), np.abs(yc - y_min), np.abs(yc - y_max)])
        large_axis = rng.uniform(min_axis_length, dist_to_barrier)
        small_axis = rng.uniform(min_axis_length, large_axis)
        theta = rng.uniform(0, 2 * np.pi)
        # Rotation matrix
        c, s = np.cos(theta), np.sin(theta)
        R = np.array([[c, -s], [s, c]])
        
        # Diagonal matrix of eigenvalues (1/a^2 and 1/b^2)
        D = np.diag([1.0 / (large_axis**2), 1.0 / (small_axis**2)])
        
        # Quadratic form matrix M = R * D * R^T
        M = R @ D @ R.T
        
        A = float(M[0, 0])
        B = float(M[0, 1])
        C = float(M[1, 1])
        
        output.extend([xc, yc, A, B, C])
    
    return output


if __name__ == '__main__':
    x_max = 1.0
    y_max = 1.0
    h = 0.02
    num_ellipses = 4
    print(generate_configuration(x_max, y_max, h, num_ellipses))
