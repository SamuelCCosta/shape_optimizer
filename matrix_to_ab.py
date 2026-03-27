import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# The input list provided
raw_data =   [   0.53446306,    0.50259374,  243.77396884,  -33.96889496,  261.82802112,
    0.42235164,    0.70059529,  233.99220896,  -43.95489961,  288.09370986,
    0.56719944,    0.34218724,  289.65241533, -106.94255163,  203.10652993,
    0.46513083,    0.8759024,   208.72228596,   71.36566829,  130.41681181]


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
    print_res = False
    if print_res:
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

    #ax.plot(0.428928, 0.442505, 'k+', markersize=5)
    #ax.plot(0.427547, 0.492422, 'k+', markersize=5)
    plt.title('Visualized Ellipses')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.show()