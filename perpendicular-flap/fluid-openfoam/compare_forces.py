import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import re
import os

def read_force_data(file_path):
    """
    Reads max force magnitude data from a file, taking the last entry for each timestep.
    """
    if not os.path.exists(file_path):
        print(f"Warning: Could not find file {file_path}")
        return None, None
        
    try:
        forces_raw = np.loadtxt(file_path)
        if forces_raw.ndim == 1:
            forces_raw = forces_raw.reshape(1, -1)
    except ValueError:
        print(f"Warning: Could not parse file {file_path}.")
        return None, None

    force_timesteps_raw = forces_raw[:, 0]
    
    force_dict = {}
    for i, t in enumerate(force_timesteps_raw):
        force_dict[t] = forces_raw[i, 2]

    timesteps = np.array(sorted(force_dict.keys()))
    force_max_magnitudes = np.array([force_dict[t] for t in timesteps])
    
    return timesteps, force_max_magnitudes

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

def plot_force_comparison(case_data):
    """
    Plots the smoothed max force magnitude from multiple cases on the same graph.
    """
    plt.figure(figsize=(14, 7))
    
    start_time = 1.0

    # Find the minimum number of timesteps to use as a baseline for smoothing
    min_len = float('inf')
    for case_name in case_data:
        times, forces = case_data[case_name]
        if times is not None and len(times) > 0:
            min_len = min(min_len, len(times))

    if min_len == float('inf'):
        print("Error: No valid data found to plot.")
        return

    sorted_cases = sorted(case_data.keys(), key=get_sort_keys)

    for case_name in sorted_cases:
        times, forces = case_data[case_name]
        if times is not None and forces is not None and len(times) > 0:
            
            # --- Apply Smoothing ---
            window_size = int(round(len(times) / min_len))
            if window_size > 1:
                smoothed_forces = np.convolve(forces, np.ones(window_size)/window_size, mode='valid')
                # Adjust times array to match the shortened smoothed_forces array
                smoothed_times = times[window_size-1:]
            else:
                smoothed_forces = forces
                smoothed_times = times

            # --- Apply Start Time ---
            start_indices = np.where(smoothed_times >= start_time)[0]
            if start_indices.size > 0:
                start_index = start_indices[0]
                z_order = 2 if 'subcycling' not in case_name else 1
                plt.plot(smoothed_times[start_index:], smoothed_forces[start_index:], label=f'{case_name}', zorder=z_order)
            else:
                print(f"Warning: No data points for case '{case_name}' after time {start_time}s.")

    plt.xlabel('Time (s)')
    plt.ylabel('Smoothed Max Force Magnitude')
    plt.title('Comparison of Smoothed Max Force Magnitude Over Time')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.ylim(12.5, 16.5)

    output_filename = 'force_comparison.png'
    plt.savefig(output_filename)
    print(f"Saved comparison plot to {output_filename}")
    plt.show()

if __name__ == "__main__":
    case_folders = [
        "reference-dt-0.01", "dt-0.005", "dt-0.02",
        "dt-0.02-subcycling-dt-0.01", "dt-0.01-subcycling-dt-0.005",
        "dt-0.01-subcycling-dt-0.00334", "dt-0.005-subcycling-dt-0.0025",
        "dt-0.01-subcycling-dt-0.001-waveform-0",
        "dt-0.005-subcycling-dt-0.0025-waveform-0",
        "dt-0.01-subcycling-dt-0.005-waveform-0"
    ]
    case_data = {}

    for folder in case_folders:
        file_path = os.path.join(folder, 'preciceForceWrite.dat')
        print(f"Reading data from: {file_path}")
        times, forces = read_force_data(file_path)
        case_data[folder] = (times, forces)

    plot_force_comparison(case_data)
