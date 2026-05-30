import matplotlib.pyplot as plt
import numpy as np

def plot_segment_statistics():
    # Data extracted from the solver output
    steps = [0, 1, 3, 5, 10, 20]
    
    min_vals = [0.00969412, 0.0112946, 0.0106927, 0.0105606, 0.0104564, 0.0105588]
    max_vals = [0.0392085, 0.038831, 0.0375823, 0.0373872, 0.0373203, 0.0373041]
    avg_vals = np.array([0.0233855, 0.0233166, 0.0232983, 0.0232933, 0.0232888, 0.0232865])
    med_vals = [0.0253147, 0.0251546, 0.0248922, 0.0248073, 0.0247472, 0.024766]
    std_vals = np.array([0.00318602, 0.00295168, 0.00293004, 0.00292813, 0.00292922, 0.00293115])

    plt.figure(figsize=(10, 6))

    # Plot the lines
    plt.plot(steps, max_vals, marker='o', label='Max', color='red', linewidth=2)
    plt.plot(steps, med_vals, marker='s', label='Median', color='orange', linewidth=2)
    plt.plot(steps, avg_vals, marker='^', label='Average', color='green', linewidth=2)
    plt.plot(steps, min_vals, marker='v', label='Min', color='blue', linewidth=2)

    # Add a shaded region for the standard deviation around the average
    plt.fill_between(steps, avg_vals - std_vals, avg_vals + std_vals, 
                     color='green', alpha=0.2, label='Avg ± Std Dev')

    # Formatting the plot
    plt.title('Mesh Segment Length Statistics vs. Smooth Steps', fontsize=14)
    plt.xlabel('Smooth Steps', fontsize=12)
    plt.ylabel('Segment Length', fontsize=12)
    plt.xticks(steps)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()

    plt.savefig('segment_length_stats.png', dpi=300)
    print("Plot saved successfully to 'segment_length_stats.png'")

if __name__ == "__main__":
    plot_segment_statistics()
