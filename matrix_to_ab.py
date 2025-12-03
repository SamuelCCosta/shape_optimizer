import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# The input list provided
raw_data = [0.09846759459531274, 0.8779145670536423, 360.2709520038292, 112.38584975360911, 393.63233210725235, 0.5056591695680825, 0.3407869436316506, 162.84307818263258, -62.60305881718613, 274.44694581695825, 0.766743884499224, 0.46129767698468077, 186.3704649771612, -144.589190636172, 227.2284762469888, 0.7819328169643374, 0.6512709822028184, 309.4108476099627, 195.55068464492143, 338.3886609589487]


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