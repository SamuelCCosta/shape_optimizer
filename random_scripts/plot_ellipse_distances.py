import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize # type: ignore

# Define the centers
c1 = np.array([0.3, 0.5])
c2 = np.array([0.7, 0.5])

# Define two different quadratic forms
Q1 = np.array([[27.777, -16.666], [-16.666, 27.777]])
Q2 = np.array([[27.777, -16.666], [-16.666, 27.777]])

# Decompose each Q into parametrize matrix
def get_transform_matrix(Q):
    w, R = np.linalg.eigh(Q)
    return R @ np.diag(1.0 / np.sqrt(w))

M1 = get_transform_matrix(Q1)
M2 = get_transform_matrix(Q2)

t_vals = np.linspace(0, 2 * np.pi, 100)
T1, T2 = np.meshgrid(t_vals, t_vals)
Z = np.zeros_like(T1)

# Calculate distance ||p1 - p2|| for each pair
for i in range(100):
    for j in range(100):
        pt1 = c1 + M1 @ np.array([np.cos(T1[i, j]), np.sin(T1[i, j])])
        pt2 = c2 + M2 @ np.array([np.cos(T2[i, j]), np.sin(T2[i, j])])
        Z[i, j] = np.linalg.norm(pt1 - pt2)

#find minimum
import numpy as np

def adjacent_angles(theta: float) -> list[float]:
    angles = [0.0, np.pi / 2, np.pi, 3 * np.pi / 2]
    
    # Range of arctan2 is [-pi, pi], ensures positive index wrap-around
    index = int(np.floor(theta * 2.0 / np.pi))
    index = (index % 4 + 4) % 4
    
    return [angles[index], angles[(index + 1) % 4], theta]

def starting_parameters(c1: np.ndarray, M1: np.ndarray, 
                        c2: np.ndarray, M2: np.ndarray) -> tuple[float, float]:
    diff = c2 - c1

    dir1 = np.linalg.inv(M1) @ diff
    t1_center = np.arctan2(dir1[1], dir1[0])

    dir2 = np.linalg.inv(M2) @ (-diff)
    t2_center = np.arctan2(dir2[1], dir2[0])

    def point_at(c, M, t):
        return c + M @ np.array([np.cos(t), np.sin(t)])

    theta1, theta2 = t1_center, t2_center
    min_start_dist_sq = np.sum((point_at(c1, M1, theta1) - point_at(c2, M2, theta2))**2)

    for t1 in adjacent_angles(t1_center):
        p1 = point_at(c1, M1, t1)
        for t2 in adjacent_angles(t2_center):
            p2 = point_at(c2, M2, t2)
            dist_sq = np.sum((p1 - p2)**2)
            
            if dist_sq <= min_start_dist_sq:
                min_start_dist_sq = dist_sq
                theta1 = t1
                theta2 = t2

    return theta1 % (2 * np.pi), theta2 % (2 * np.pi)


def dist_func(t):
    pt1 = c1 + M1 @ np.array([np.cos(t[0]), np.sin(t[0])])
    pt2 = c2 + M2 @ np.array([np.cos(t[1]), np.sin(t[1])])
    return np.linalg.norm(pt1 - pt2)

starting_point = starting_parameters(c1, M1, c2, M2)
print(f'{starting_point=}')

res = minimize(dist_func, x0=starting_point) 
t1_spec, t2_spec = res.x
print(f'{t1_spec=}, {t2_spec=}')
z_spec = res.fun

# Plotting
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(T1, T2, Z, cmap='plasma', edgecolor='none', alpha=0.8) # type: ignore
ax.scatter(t1_spec, t2_spec, z_spec, color='cyan', s=100, label='Specific Point (Min)', edgecolors='black', depthshade=False) # type: ignore
ax.scatter(starting_point[0], starting_point[1], 0.0, color='red', s=100, label='Starting Point', edgecolors='black', depthshade=False) # type: ignore

ax.set_xlabel('t1 (Ellipse 1)')
ax.set_ylabel('t2 (Ellipse 2)')
ax.set_zlabel('Distance ||p1 - p2||') # type: ignore


fig.colorbar(surf, shrink=0.5, aspect=5, label='Distance')
plt.show()