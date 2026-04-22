import numpy as np
import matplotlib.pyplot as plt

# Define the centers
c1 = np.array([0.3, 0.5])
c2 = np.array([0.7, 0.5])

# Define two different quadratic forms
Q1 = np.array([[27.777, -16.666], [-16.666, 27.777]])
Q2 = np.array([[27.777, -16.666], [-16.666, 27.777]])

# Decompose each Q into valores próprios and vetores próprios
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

# Plotting
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(T1, T2, Z, cmap='plasma', edgecolor='none') # type: ignore

ax.set_xlabel('t1 (Ellipse 1)')
ax.set_ylabel('t2 (Ellipse 2)')
ax.set_zlabel('Distance ||p1 - p2||') # type: ignore

fig.colorbar(surf, shrink=0.5, aspect=5, label='Distance')
plt.show()