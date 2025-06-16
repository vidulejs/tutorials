import matplotlib.pyplot as plt
import numpy as np

def parse_force_file(filepath):
    """
    Reads an OpenFOAM-style force file and extracts the numerical values.
    
    Args:
        filepath (str): The path to the force file.
        
    Returns:
        list: A list of floats containing the force values, or None if the file is not found.
    """
    forces = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                # Check if the line contains the expected text 'average force='
                if 'average force=' in line:
                    # Split the line at '=' and take the second part
                    value_str = line.split('=')[1]
                    # Convert the string to a float and add it to our list
                    forces.append(float(value_str.strip()))
    except FileNotFoundError:
        print(f"Error: The file '{filepath}' was not found.")
        return None
    except Exception as e:
        print(f"An error occurred while reading {filepath}: {e}")
        return None
        
    return forces

# --- Main Script ---

# 1. Define the paths to your force files
file_subcycling = 'euler-subcycling3/force'
file_reference = 'euler-reference/force'

# 2. Parse the data from the files
forces_subcycling = parse_force_file(file_subcycling)
forces_reference = parse_force_file(file_reference)

# 3. Check if data was loaded successfully before plotting
if forces_subcycling is not None and forces_reference is not None:
    
    # Create an index for the x-axis (e.g., 0, 1, 2, 3...)
    # This assumes each data point is a sequential step.
    x_subcycling = np.arange(len(forces_subcycling))
    x_reference = np.arange(len(forces_reference))
    
    # 4. Create the plot
    plt.figure(figsize=(12, 7)) # Create a figure with a nice size
    
    # Plot both series
    plt.plot(x_subcycling, forces_subcycling, label='Subcycling (3 steps)', marker='o', linestyle='-', markersize=4)
    plt.plot(x_reference, forces_reference, label='Reference', marker='x', linestyle='--', markersize=4)
    
    # 5. Add labels, title, legend, and a grid for readability
    plt.xlabel('Iteration / Time Step')
    plt.ylabel('Average Force')
    plt.legend()
    plt.ylim(0, 500)
    plt.grid(True)
    
    # 6. Show the plot
    plt.show()
