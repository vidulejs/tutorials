#!/usr/bin/env python
# coding: utf-8

# In[1]:


import precice
import deepxde as dde
import numpy as np

# pde setup
from pde_geom import *
from net import *

precice_config = '../precice-config-training.xml'
solver_process_index = 0
solver_process_size = 1

participant_name = "Neural"
mesh_name = "Neural-Mesh"
participant = precice.Participant(participant_name, precice_config, solver_process_index, solver_process_size)


# In[2]:


data = define_pde(geom, navier_stokes_pde, loss_terms)
model = dde.Model(data, net)

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


# In[3]:


def model_predict(x_mesh, y_mesh, time=0.0):
    n_points_mesh = len(x_mesh.flatten())
    t_values_mesh = np.full(n_points_mesh, time)
    
    # Input to model: (x, y, t)
    points_to_predict = np.vstack((x_mesh.flatten(), y_mesh.flatten(), t_values_mesh)).T
    
    # Model predicts (u, v, p)
    predictions = model.predict(points_to_predict)
    
    u_pred = predictions[:, 0].reshape(x_mesh.shape)
    v_pred = predictions[:, 1].reshape(x_mesh.shape)
    p_pred = predictions[:, 2].reshape(x_mesh.shape)

    Velocity = np.column_stack([u_pred.flatten(), v_pred.flatten()])
    Pressure = p_pred.flatten()
    
    return Velocity, Pressure


# In[ ]:


# should match openfoam mesh
x_min, x_max = 0.0, domain_length
y_min, y_max = 0.0, domain_width
nx, ny = 60, 20

x_centers = np.linspace(x_min + (x_max - x_min)/(2*nx), x_max - (x_max - x_min)/(2*nx), nx)
y_centers = np.linspace(y_min + (y_max - y_min)/(2*ny), y_max - (y_max - y_min)/(2*ny), ny)

X_grid, Y_grid = np.meshgrid(x_centers, y_centers)
coordinates_precice = np.vstack([X_grid.flatten(), Y_grid.flatten()]).T

vertex_ids_precice = participant.set_mesh_vertices(mesh_name, coordinates_precice)

Velocity, Pressure = model_predict(X_grid, Y_grid, time=0.0)

dt_precice = participant.initialize()
time = 0.0


# In[ ]:


# while participant.is_coupling_ongoing():
    
#     # The PINN model can predict at any t, so we advance by dt_precice.
#     print(f"Advancing PINN from t={time:.5f} (dt={dt_precice:.5f})")
#     time += dt_precice

#     Velocity, Pressure = model_predict(X_grid, Y_grid, time=time)
        
#     Velocity_ground_truth = participant.read_data(mesh_name, "Velocity", vertex_ids_precice, Velocity)
#     Pressure_ground_truth = participant.read_data(mesh_name, "Pressure", vertex_ids_precice, Pressure)

#     participant.advance(dt_precice)

#     dt_precice = participant.get_max_time_step_size()
    
# participant.finalize()
# print("preCICE coupling finished.")


# In[ ]:


V_dataset = []
P_dataset = []
T_dataset = []


# In[ ]:


while participant.is_coupling_ongoing():
    
    time += dt_precice
        
    Velocity_ground_truth = participant.read_data(mesh_name, "Velocity", vertex_ids_precice, Velocity)
    Pressure_ground_truth = participant.read_data(mesh_name, "Pressure", vertex_ids_precice, Pressure)

    V_dataset.append(Velocity_ground_truth)
    P_dataset.append(Pressure_ground_truth)
    T_dataset.append(time)

    participant.advance(dt_precice)

    dt_precice = participant.get_max_time_step_size()
    
participant.finalize()
print("preCICE coupling finished.")


# In[2]:


get_ipython().system("jupyter nbconvert --to script 'training-participant.ipynb'")

