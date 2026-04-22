import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle
import numpy as np

fig, ax = plt.subplots(figsize=(8, 8))

e1 = Ellipse((0.3, 0.5), width=0.6, height=0.3, angle=45, edgecolor='b', fc='none', lw=2, label='E1')
e2 = Ellipse((0.7, 0.5), width=0.6, height=0.3, angle=45, edgecolor='r', fc='none', lw=2, label='E2')

ax.add_patch(e1)
ax.add_patch(e2)

domain = Rectangle((0.02, 0.02), 0.96, 0.96, edgecolor='k', fc='none', ls='--', label='Domain')
ax.add_patch(domain)

def plot_tips(cx, cy, a, b, ang_deg, col):
    rad = np.radians(ang_deg)
    c, s = np.cos(rad), np.sin(rad)
    ax.plot([cx + a*c, cx - a*c], [cy + a*s, cy - a*s], 'o', color=col)
    ax.plot([cx - b*s, cx + b*s], [cy + b*c, cy - b*c], 'o', color=col)

plot_tips(0.3, 0.5, 0.3, 0.15, 45, 'b')
plot_tips(0.7, 0.5, 0.3, 0.15, 45, 'r')

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect('equal')
ax.grid(True, ls=':')
ax.legend()
plt.show()