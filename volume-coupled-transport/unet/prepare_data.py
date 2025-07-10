import numpy as np
import os

# --- Configuration ---
# These values should match the grid dimensions from your simulation.
# Based on your scripts, they are:
nx = 60
ny = 40
output_filename = 'training_data_gridded.npz'
input_filename = '../fluid-deepxde/training_data.npz'

print(f"Loading pointwise data from: {input_filename}")
try:
    data_npz = np.load(input_filename)
    X_data = data_npz['X_train'].astype(np.float32)
    Y_data = data_npz['Y_train'].astype(np.float32)
except FileNotFoundError:
    print(f"Error: The file '{input_filename}' was not found.")
    print("Please ensure the path is correct and you have run the necessary data generation scripts.")
    exit()

print("Original X data shape:", X_data.shape)
print("Original Y data shape:", Y_data.shape)

time_steps = np.unique(X_data[:, 2])
num_time_steps = len(time_steps)
points_per_timestep = X_data.shape[0] // num_time_steps

print(f"Found {num_time_steps} unique time steps.")
print(f"Points per time step: {points_per_timestep}")

if points_per_timestep != nx * ny:
    print(f"Error: The number of points per time step ({points_per_timestep}) "
          f"does not match the expected grid size of nx*ny ({nx*ny}).")
    exit()

num_variables = Y_data.shape[1]
gridded_Y = Y_data.reshape(num_time_steps, ny, nx, num_variables)
gridded_Y = np.transpose(gridded_Y, (1, 2, 0, 3))
gridded_Y = np.expand_dims(gridded_Y, axis=0)

print(f"Successfully reshaped Y data into a gridded format.")
print("New gridded Y data shape:", gridded_Y.shape)

print("Normalizing data...")

means = np.mean(gridded_Y, axis=(0, 1, 2, 3), keepdims=True)
stds = np.std(gridded_Y, axis=(0, 1, 2, 3), keepdims=True)

stds = np.where(stds == 0, 1, stds)

normalized_gridded_Y = (gridded_Y - means) / stds

print("Data normalization complete.")
print("Channel means:", means.squeeze())
print("Channel stds:", stds.squeeze())
print("Final normalized data shape:", normalized_gridded_Y.shape)

np.savez(output_filename, 
         data=normalized_gridded_Y,
         means=means,
         stds=stds)
print(f"Gridded and normalized data saved to: {output_filename}")

