import matplotlib.pyplot as plt
import numpy as np
import re
import os

def read_probed_data(file_path):
    """
    Reads coordinate data for a vertex from a file, taking the last entry for each timestep.
    """
    data_dict = {}

    # Regular expression to parse the coordinate tuples
    coord_regex = re.compile(r'\(([^)]+)\)')

    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue

            try:
                time = float(parts[0])
                
                # Find all coordinate tuples in the line
                coords_found = coord_regex.findall(line)
                
                if len(coords_found) >= 2:
                    # Vertex coordinates are the second tuple
                    vertex_str = coords_found[1].split()
                    vertex_coord = [float(v) for v in vertex_str]
                    data_dict[time] = vertex_coord
                else:
                    print(f"Warning: Skipping malformed line: {line.strip()}")

            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse line '{line.strip()}'. Error: {e}")
                continue

    if not data_dict:
        print(f"Error: No valid data found in the file: {file_path}")
        return None, None

    # Sort data by time
    sorted_times = sorted(data_dict.keys())
    times = np.array(sorted_times)
    vertex_coords = np.array([data_dict[t] for t in sorted_times])
    
    return times, vertex_coords

def get_dt_from_name(case_name):
    """Extracts the timestep value from the case folder name."""
    match = re.search(r'(\d+\.\d+)', case_name)
    if match:
        return float(match.group(1))
    return float('inf') # Should not happen with current naming

def plot_comparison(case_data):
    """
    Plots the horizontal displacement from multiple cases on the same graph,
    ordered by decreasing timestep.
    """
    plt.figure(figsize=(14, 7))
    
    # Sort cases by timestep in descending order
    sorted_cases = sorted(case_data.keys(), key=get_dt_from_name, reverse=True)

    for case_name in sorted_cases:
        data = case_data[case_name]
        times, vertex_coords = data
        if times is not None and vertex_coords is not None:
            plt.plot(times, vertex_coords[:, 0], label=f'{case_name} - Horizontal Displacement')

    plt.xlabel('Time (s)')
    plt.ylabel('Horizontal Position')
    plt.title('Comparison of Horizontal Displacement Over Time')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    
    output_filename = 'displacement_comparison.png'
    plt.savefig(output_filename)
    print(f"Saved comparison plot to {output_filename}")
    plt.show()

def calculate_and_print_l2_norm_diff(base_case_data, fine_case_data, base_case_name, fine_case_name):
    """
    Interpolates, calculates, and prints the normalized L2 norm of the difference.
    """
    times_base, coords_base = base_case_data
    times_fine, coords_fine = fine_case_data

    if times_base is None or times_fine is None:
        print(f"Cannot calculate L2 norm for {fine_case_name} vs {base_case_name} due to missing data.")
        return

    # Interpolate the base case data onto the fine case's time steps
    coords_base_interp = np.zeros_like(coords_fine)
    for i in range(coords_base.shape[1]): # for each dimension (x, y, z)
        coords_base_interp[:, i] = np.interp(times_fine, times_base, coords_base[:, i])

    # Calculate the L2 norm of the difference vector
    l2_norm_diff = np.linalg.norm(coords_fine - coords_base_interp)

    # Normalize by the number of time steps in the fine case
    normalized_l2_norm = l2_norm_diff / len(times_fine)

    print(f"\n--- Convergence Analysis ---")
    print(f"Comparing '{fine_case_name}' against '{base_case_name}':")
    print(f"Normalized L2 Norm of displacement difference: {normalized_l2_norm}")

if __name__ == "__main__":
    case_folders = ["reference-dt-0.01", "dt-0.005", "dt-0.001", "dt-0.02"]
    case_data = {}

    for folder in case_folders:
        file_path = os.path.join(folder, 'probedLocations.dat')
        if os.path.exists(file_path):
            print(f"Reading data from: {file_path}")
            times, coords = read_probed_data(file_path)
            case_data[folder] = (times, coords)
        else:
            print(f"Error: Could not find file {file_path}")

    if len(case_data) > 1:
        plot_comparison(case_data)
        
        # Define pairs for comparison (base, fine)
        comparison_pairs = [
            ("dt-0.02", "reference-dt-0.01"),
            ("reference-dt-0.01", "dt-0.005"),
            ("dt-0.005", "dt-0.001"),
            ("reference-dt-0.01", "dt-0.001")
        ]

        for base, fine in comparison_pairs:
            if base in case_data and fine in case_data:
                calculate_and_print_l2_norm_diff(
                    case_data[base], 
                    case_data[fine],
                    base,
                    fine
                )
    else:
        print("Could not read data for enough cases to make a comparison.")