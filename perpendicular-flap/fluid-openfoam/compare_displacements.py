import matplotlib.pyplot as plt
import numpy as np
import re
import os

def read_probed_data(file_path):
    """
    Reads coordinate data for a vertex from a file, taking the last entry for each timestep.
    """
    data_dict = {}
    coord_regex = re.compile(r'\(([^)]+)\)')
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts: continue
            try:
                time = float(parts[0])
                coords_found = coord_regex.findall(line)
                if len(coords_found) >= 2:
                    vertex_str = coords_found[1].split()
                    vertex_coord = [float(v) for v in vertex_str]
                    data_dict[time] = vertex_coord
            except (ValueError, IndexError):
                continue
    if not data_dict: return None, None
    sorted_times = sorted(data_dict.keys())
    times = np.array(sorted_times)
    vertex_coords = np.array([data_dict[t] for t in sorted_times])
    return times, vertex_coords

def get_sort_keys(case_name):
    """
    Generates a tuple for sorting cases for the plot.
    """
    numbers = [float(n) for n in re.findall(r'(\d+\.\d+)', case_name)]
    is_subcycling = 'subcycling' in case_name
    group = 1 if is_subcycling else 0
    val1 = -numbers[0] if numbers else 0
    val2 = -numbers[1] if len(numbers) > 1 and is_subcycling else 0
    return (group, val1, val2)

def plot_comparison(case_data):
    """
    Plots the horizontal displacement from multiple cases on the same graph.
    """
    plt.figure(figsize=(14, 7))
    sorted_cases = sorted(case_data.keys(), key=get_sort_keys)
    for case_name in sorted_cases:
        data = case_data[case_name]
        times, vertex_coords = data
        if times is not None and vertex_coords is not None:
            plt.plot(times, vertex_coords[:, 0], label=f'{case_name}')
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
    if times_base is None or times_fine is None: return

    coords_base_interp = np.zeros_like(coords_fine)
    for i in range(coords_base.shape[1]):
        coords_base_interp[:, i] = np.interp(times_fine, times_base, coords_base[:, i])
    
    l2_norm_diff = np.linalg.norm(coords_fine - coords_base_interp)
    normalized_l2_norm = l2_norm_diff / len(times_fine)

    print(f"\n--- Comparing '{fine_case_name}' against '{base_case_name}' ---")
    print(f"Normalized L2 Norm: {normalized_l2_norm:.2e}")

if __name__ == "__main__":
    case_folders = [
        "reference-dt-0.01", #"dt-0.005", 
        #"dt-0.001", "dt-0.02",
        # "dt-0.02-subcycling-dt-0.01", "dt-0.01-subcycling-dt-0.005",
        # "dt-0.01-subcycling-dt-0.00334", "dt-0.005-subcycling-dt-0.0025",
        "dt-0.01-subcycling-dt-0.001-waveform-0",
        # "dt-0.005-subcycling-dt-0.0025-waveform-0",
        # "dt-0.01-subcycling-dt-0.005-waveform-0",
        "dt-0.01-subcycling-dt-0.001-waveform-1",
        "dt-0.01-subcycling-dt-0.001-waveform-2",
    ]
    case_data = {}
    for folder in case_folders:
        file_path = os.path.join(folder, 'probedLocations.dat')
        if os.path.exists(file_path):
            print(f"Reading data from: {file_path}")
            times, coords = read_probed_data(file_path)
            case_data[folder] = (times, coords)

    if len(case_data) > 1:
        plot_comparison(case_data)
        
        # Define all pairs for comparison (base, fine)
        comparison_pairs = [
            # Timestep convergence vs finer case
            ("dt-0.02", "reference-dt-0.01"),
            ("reference-dt-0.01", "dt-0.005"),
            ("dt-0.005", "dt-0.001"),
            
            # Timestep cases vs reference
            ("reference-dt-0.01", "dt-0.02"),
            ("reference-dt-0.01", "dt-0.001"),

            # Subcycling cases vs reference
            ("reference-dt-0.01", "dt-0.02-subcycling-dt-0.01"),
            ("reference-dt-0.01", "dt-0.01-subcycling-dt-0.005"),
            ("reference-dt-0.01", "dt-0.01-subcycling-dt-0.00334"),
            ("reference-dt-0.01", "dt-0.005-subcycling-dt-0.0025"),
            ("reference-dt-0.01", "dt-0.01-subcycling-dt-0.001-waveform-0"),
            ("reference-dt-0.01", "dt-0.01-subcycling-dt-0.001-waveform-1"),

            # Subcycling cases vs each other (ordered by plot sort key)
            ("dt-0.02-subcycling-dt-0.01", "dt-0.01-subcycling-dt-0.005"),
            ("dt-0.01-subcycling-dt-0.005", "dt-0.01-subcycling-dt-0.00334"),
            ("dt-0.01-subcycling-dt-0.00334", "dt-0.005-subcycling-dt-0.0025"),
        ]

        print("\n\n--- CONVERGENCE ANALYSIS ---")
        for base, fine in comparison_pairs:
            if base in case_data and fine in case_data:
                calculate_and_print_l2_norm_diff(
                    case_data[base], case_data[fine], base, fine
                )
    else:
        print("Could not read data for enough cases to make a comparison.")
