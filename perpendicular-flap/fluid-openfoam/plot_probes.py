
import matplotlib.pyplot as plt
import numpy as np
import re

def plot_probed_data(file_path):
    """
    Reads and plots coordinate data for a vertex from a file.
    """
    times = []
    vertex_coords = []

    # Regular expression to parse the coordinate tuples
    coord_regex = re.compile(r'\(([^)]+)\)')

    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue

            try:
                times.append(float(parts[0]))
                
                # Find all coordinate tuples in the line
                coords_found = coord_regex.findall(line)
                
                if len(coords_found) >= 2:
                    # Vertex coordinates are the second tuple
                    vertex_str = coords_found[1].split()
                    vertex_coords.append([float(v) for v in vertex_str])
                else:
                    print(f"Warning: Skipping malformed line: {line.strip()}")

            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse line '{line.strip()}'. Error: {e}")
                continue

    if not times:
        print("Error: No valid data found in the file.")
        return

    # Convert lists to numpy arrays for easier slicing
    vertex_coords = np.array(vertex_coords)

    # --- Plot for Vertex ---
    n_points = len(times)
    # Use the 'magma' colormap for a clear visual transition
    colors = plt.cm.magma(np.linspace(0, 1, n_points))

    plt.figure(figsize=(14, 7))
    # Use a neutral color for the line to emphasize the points
    plt.plot(vertex_coords[:, 0], times, color='grey', alpha=0.5, zorder=1)
    # Color scatter points by their order
    scatter = plt.scatter(vertex_coords[:, 0], times, c=np.arange(n_points), cmap='magma', zorder=2)
    
    # Add a color bar to the side
    cbar = plt.colorbar(scatter, label='Step number')
    cbar.set_label('Step number', size=12)

    plt.axhline(y=0.1, color='b', linestyle='--', label='Time = 0.1s')
    plt.xlabel('Horizontal Position')
    plt.ylabel('Time (s)')
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.savefig('vertex_coordinate.png')
    print("Saved plot to vertex_coordinate.png")
    plt.show()


if __name__ == "__main__":
    plot_probed_data('probedLocations.dat')
