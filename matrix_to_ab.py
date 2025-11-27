import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# The input list provided
raw_data = [ 
  1.22818673e-01,  5.83517766e-01,  1.12806305e+02,  1.21366115e+01, 8.96040620e+01,
  5.91629030e-01,  6.40340026e-01,  2.04159081e+02, -4.35976084e+01, 1.12926963e+02,
  5.42716448e-01,  7.59968556e-01,  1.00550021e+02, -1.30988244e+01, 1.17660615e+02,
  3.40343817e-01,  1.44695342e-01,  9.75079343e+01, -1.34730584e+01,  2.13321885e+02]

#as elipses 2 e 3 intersetam *e* passam os testes

def get_ellipse_parameters(data):
    """
    Parses a list of ellipse parameters and returns geometric properties.
    Output format per ellipse: (a, b, theta_radians, center_x, center_y)
    """
    ellipses = []
    
    # Iterate through list in chunks of 5
    for i in range(0, len(data), 5):
        # Unpack elements
        xc, yc = data[i], data[i+1]
        A, B, C = data[i+2], data[i+3], data[i+4]
        
        M = np.array([[A, B], [B, C]])
        
        # 2. Eigen Decomposition
        w, v = np.linalg.eigh(M)
        
        # 3. Sort Eigenvalues
        # The Major axis (a) corresponds to the SMALLEST eigenvalue
        # The Minor axis (b) corresponds to the LARGEST eigenvalue
        idx = w.argsort()
        w = w[idx]
        v = v[:, idx] # Reorder eigenvectors to match
        
        lambda_small = w[0]
        lambda_large = w[1]
        
        a = 1.0 / np.sqrt(lambda_small) # Semi-major
        b = 1.0 / np.sqrt(lambda_large) # Semi-minor
        
        # 5. Calculate Orientation (Theta)
        v_major = v[:, 0] 
        theta = np.arctan2(v_major[1], v_major[0])
        if theta < 0:
            theta += np.pi
        
        ellipses.append({
            'a': a,
            'b': b,
            'theta': theta,
            'center': (xc, yc)
        })
        
    return ellipses


if __name__ == '__main__':
    # Run the function
    results = get_ellipse_parameters(raw_data)

    # Print results
    print(f"{'Major (a)':<12} | {'Minor (b)':<12} | {'Theta (rad)':<12} | {'Center (x, y)'}")
    print("-" * 65)

    for res in results:
        print(f"{res['a']:<12.6f} | {res['b']:<12.6f} | {res['theta']:<12.6f} | {res['center']}")

    fig, ax = plt.subplots(figsize=(6, 6))

    
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_aspect('equal')
    ax.grid(True, linestyle='--', alpha=0.6)

    colors = ['red', 'blue', 'green', 'purple']

    for i, p in enumerate(results):
        angle_deg = np.degrees(p['theta'])
        
        # Matplotlib Ellipse takes: xy, width, height, angle
        # Note: width = 2*a, height = 2*b
        e = Ellipse(
            xy=p['center'], 
            width=p['a'] * 2, 
            height=p['b'] * 2, 
            angle=angle_deg,
            edgecolor=colors[i % len(colors)],
            facecolor='none',
            lw=2,
            label=f'Ellipse {i+1}'
        )
        ax.add_patch(e)
        
        # Plot the center point for reference
        ax.plot(p['center'][0], p['center'][1], 'k+', markersize=5)

    plt.title('Visualized Ellipses')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.show()