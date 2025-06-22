#!/usr/bin/env python
# coding: utf-8

# In[1]:


import precice
import deepxde as dde
import numpy as np

# pde setup
from pde_geom import *
from net import *

precice_config = '../precice-config.xml'
solver_process_index = 0
solver_process_size = 1

participant_name = "Feedback"
mesh_name = "Feedback-Mesh"
participant = precice.Participant(participant_name, precice_config, solver_process_index, solver_process_size)


# In[2]:


# should match openfoam mesh
x_min, x_max = 0.0, domain_length
y_min, y_max = 0.0, domain_width
nx, ny = 60, 40

x_centers = np.linspace(x_min + (x_max - x_min)/(2*nx), x_max - (x_max - x_min)/(2*nx), nx)
y_centers = np.linspace(y_min + (y_max - y_min)/(2*ny), y_max - (y_max - y_min)/(2*ny), ny)

X_grid, Y_grid = np.meshgrid(x_centers, y_centers)
coordinates_precice = np.vstack([X_grid.flatten(), Y_grid.flatten()]).T

vertex_ids_precice = participant.set_mesh_vertices(mesh_name, coordinates_precice)

participant.initialize()
time = 0.0


# In[3]:


dt_precice = participant.get_max_time_step_size()


# In[4]:


V_dataset = []
P_dataset = []
T_dataset = []


# In[5]:


velocity_buffer = np.zeros((len(vertex_ids_precice)))


# In[6]:


participant.read_data


# In[7]:


while participant.is_coupling_ongoing():
    
    time += dt_precice
        
    Velocity_ground_truth = participant.read_data(mesh_name, "Velocity", vertex_ids_precice, 0)
    Pressure_ground_truth = participant.read_data(mesh_name, "Pressure", vertex_ids_precice, 0)

    print(f"Time: {time:.3f}, dt_precice: {dt_precice:.3f}")

    V_dataset.append(np.copy(Velocity_ground_truth))
    P_dataset.append(np.copy(Pressure_ground_truth))
    T_dataset.append(time)

    participant.advance(dt_precice)

    dt_precice = participant.get_max_time_step_size()
    
participant.finalize()
print("preCICE coupling finished.")


# In[8]:


dt_precice = 0.005


# In[9]:


import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots(figsize=(10, 6))

def magnitude(vector_field):
    reshaped = vector_field.reshape((ny, nx, 2))
    x = reshaped[:, :, 0]
    y = reshaped[:, :, 1]
    return np.sqrt(x**2 + y**2)

# Function to update the plot for each frame
def update(frame):
    ax.clear()
    velocity_mag = magnitude(V_dataset[frame])
    im = ax.imshow(velocity_mag, cmap='coolwarm', origin='lower', 
                   extent=(x_min, x_max, y_min, y_max),
                   vmin=0, vmax=2)  # Consistent color range
    
    # Add obstacle representation
    ax.add_patch(plt.Rectangle((obstacle_x0, obstacle_y0), 
                              obstacle_x1-obstacle_x0, 
                              obstacle_y1-obstacle_y0,
                              color='white', alpha=0.7))
    
    ax.set_title(f'Velocity Magnitude at t = {T_dataset[frame]:.2f}s')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    return im,

update(len(V_dataset) - 1)

# Add colorbar
cbar = plt.colorbar(update(len(V_dataset) - 1)[0])
cbar.set_label('Velocity Magnitude')

plt.tight_layout()
plt.show()


# In[10]:


fps = 20 
total_frames = len(V_dataset)
frames_to_skip = total_frames // (fps * 10)

ani = FuncAnimation(fig, update, frames=range(0, total_frames, frames_to_skip), 
                    interval=1000 // fps, blit=True)

import os
if not os.path.exists('velocity_magnitude_animation.mp4'):
    print("Saving animation...")
    ani.save('velocity_magnitude_animation.mp4', writer='ffmpeg', fps=fps)


# # Reshape the dataset for training

# In[11]:


# X_train (x, y, t) coordinates input
(x, y) = (coordinates_precice[:, 0:1], coordinates_precice[:, 1:2])
# tile the coordinates over time
(x, y) = np.tile(x, (len(T_dataset), 1)), np.tile(y, (len(T_dataset), 1))
t = np.tile(np.array(T_dataset)[:, None], (len(coordinates_precice), 1))

X_train = np.hstack((x, y, t))

print(f"X_train shape: {X_train.shape}")
# y_train (Vx, Vy, P) outputs

all_V = np.vstack(V_dataset)
all_P = np.hstack(P_dataset)
all_P = all_P[:, None] # Reshape P to be a column vector

# Combine them to get the final (u, v, p) output array
Y_train = np.hstack((all_V, all_P))
print(f"Y_train shape: {Y_train.shape}")

# sample part of the data
sample_size = 2**16
indices = np.random.choice(len(X_train), sample_size, replace=False)
X_train = X_train[indices]
Y_train = Y_train[indices]


# # Define model

# In[12]:


unsupervised_loss_terms


# In[13]:


import torch
import gc
torch.cuda.empty_cache()
gc.collect()


# In[14]:


supervised_loss = dde.icbc.PointSetBC(X_train, Y_train)
loss_terms = unsupervised_loss_terms + [supervised_loss]
loss_weights = unsupervised_loss_weights + [10.0]

data = dde.data.TimePDE(
    geom,
    navier_stokes_pde,
    loss_terms,
    num_domain=num_domain_points,
    num_boundary=num_boundary_points,
    num_initial=num_initial_points,
    num_test=num_test_points
)

model = dde.Model(data, net)

lr = 1e-4
model.compile(
    "adam",
    lr=lr,
    loss_weights=loss_weights
)

# # Load  trained model
# model_path = "navier_stokes_channel_obstacle_model-600.pt"

# if not os.path.exists(model_path):
#     print(f"Error: Model file {model_path} not found!")
#     sys.exit(1)
# model.restore(save_path=model_path, verbose=1)


# In[15]:


def model_predict(x_mesh, y_mesh, time=0.0):
    n_points_mesh = len(x_mesh.flatten())
    t_values_mesh = np.full(n_points_mesh, time)
    
    # Input to model: (x, y, t)
    points_to_predict = np.vstack((x_mesh.flatten(), y_mesh.flatten(), t_values_mesh)).T
    
    # Model predicts (u, v, p)
    predictions = model.predict(points_to_predict)
    print(f"Predictions shape: {predictions.shape}")
    
    u_pred = predictions[:, 0].reshape(x_mesh.shape)
    v_pred = predictions[:, 1].reshape(x_mesh.shape)
    p_pred = predictions[:, 2].reshape(x_mesh.shape)

    Velocity = np.column_stack([u_pred.flatten(), v_pred.flatten()])
    Pressure = p_pred.flatten()
    
    return Velocity, Pressure

Velocity, Pressure = model_predict(X_grid, Y_grid, time=0.0)


# In[16]:


# checkpoint_callback = dde.callbacks.ModelCheckpoint(
#     os.path.join(training_output_dir, f"{model_base_name}"),
#     save_better_only=True,
#     period=50
# )

iterations = 25000

losshistory, train_state = model.train(
    iterations=iterations,
    display_every=100
    # callbacks=[checkpoint_callback]
)


# In[28]:


model_path = os.path.join(training_output_dir, f"{model_base_name}_final")
model.save(model_path)
print(f"Model saved to: {model_path}")


# In[21]:


dde.saveplot(
    losshistory, 
    train_state, 
    issave=True, 
    isplot=True,
    output_dir=training_output_dir
)


# In[27]:


Velocity, Pressure = model_predict(X_grid, Y_grid, time=1.0)

# Plot the predicted velocity magnitude
velocity_magnitude = magnitude(Velocity)
plt.figure(figsize=(10, 6))
plt.contourf(X_grid, Y_grid, velocity_magnitude, levels=50, cmap='jet')
plt.colorbar(label='Velocity Magnitude')
plt.xlabel('X')
plt.ylabel('Y')
plt.savefig(os.path.join(training_output_dir, 'velocity_magnitude.png'))
plt.show()


# In[29]:


get_ipython().system("jupyter nbconvert --to script 'training-participant.ipynb'")


# In[ ]:




