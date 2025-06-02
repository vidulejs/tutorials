import numpy as np
import precice
import time
import sys
import os
import deepxde as dde

def main():
    precice_config = '../precice-config.xml'
    solver_process_index = 0
    solver_process_size = 1
    
    # pde setup
    from pde import geom, timedomain, net, data, domain_length, domain_width
    
    # Load  trained model
    model = dde.Model(data, net)
    model.compile("adam", lr=1e-3)
    
    model_path = "navier_stokes_channel_obstacle_model-600.pt"
    if not os.path.exists(model_path):
        print(f"Error: Model file {model_path} not found!")
        sys.exit(1)
    model.restore(save_path=model_path, verbose=1)

    def create_velocity_field(x_mesh, y_mesh, current_time=0.0):
        n_points_mesh = len(x_mesh.flatten())
        t_values_mesh = np.full(n_points_mesh, current_time)
        
        # Input to model: (x, y, t)
        points_to_predict = np.vstack((x_mesh.flatten(), y_mesh.flatten(), t_values_mesh)).T
        
        # Model predicts (u, v, p)
        predictions = model.predict(points_to_predict)
        
        # We only need u and v
        u_pred = predictions[:, 0].reshape(x_mesh.shape)
        v_pred = predictions[:, 1].reshape(x_mesh.shape) # Use x_mesh.shape for both
        
        return u_pred, v_pred

    # Mesh for preCICE coupling (can be different from DeepXDE's internal points)
    # Using the domain_length from the NS model
    x_min_precice, x_max_precice = 0.0, domain_length
    y_min_precice, y_max_precice = 0.0, domain_width
    nx_precice, ny_precice = 60, 20 # Or match your preCICE setup

    x_centers_precice = np.linspace(x_min_precice + (x_max_precice - x_min_precice)/(2*nx_precice), x_max_precice - (x_max_precice - x_min_precice)/(2*nx_precice), nx_precice)
    y_centers_precice = np.linspace(y_min_precice + (y_max_precice - y_min_precice)/(2*ny_precice), y_max_precice - (y_max_precice - y_min_precice)/(2*ny_precice), ny_precice)
    
    X_precice, Y_precice = np.meshgrid(x_centers_precice, y_centers_precice)
    coordinates_precice = np.vstack([X_precice.flatten(), Y_precice.flatten()]).T
    
    # Initialize velocity field at t=0 for preCICE
    U_init, V_init = create_velocity_field(X_precice, Y_precice, current_time=0.0)
    velocities_to_write = np.column_stack([U_init.flatten(), V_init.flatten()])
    
    participant_name = "Fluid" # Should match your precice-config.xml
    mesh_name = "Fluid-Mesh"   # Should match your precice-config.xml
    participant = precice.Participant(participant_name, precice_config, solver_process_index, solver_process_size)
    
    vertex_ids_precice = participant.set_mesh_vertices(mesh_name, coordinates_precice)
    
    dt_precice = participant.initialize() # dt_precice is the first coupling time step
    dt_precice = 0.005
    current_sim_time = 0.0
    
    # Checkpoint variables
    velocities_checkpoint = None
    time_checkpoint = None

    while participant.is_coupling_ongoing():
        if participant.requires_writing_checkpoint():
            print(f"Writing checkpoint at t={current_sim_time:.5f}")
            velocities_checkpoint = velocities_to_write.copy()
            time_checkpoint = current_sim_time
        
        # dt_precice is the coupling window size suggested by preCICE
        # The PINN model can predict at any t, so we advance by dt_precice.
        target_time = current_sim_time + dt_precice
        print(f"Advancing PINN from t={current_sim_time:.5f} to t={target_time:.5f} (dt={dt_precice:.5f})")
        
        U_new, V_new = create_velocity_field(X_precice, Y_precice, current_time=target_time)
        velocities_to_write = np.column_stack([U_new.flatten(), V_new.flatten()])
        
        # Assuming "Velocity" is the data name in precice-config.xml for vector data
        participant.write_data(mesh_name, "Velocity", vertex_ids_precice, velocities_to_write)
        
        participant.advance(dt_precice) # Advance preCICE by the coupling window
        current_sim_time = target_time # Update simulation time
        
        if participant.requires_reading_checkpoint():
            print(f"Reading checkpoint, restoring to t={time_checkpoint:.5f}")
            velocities_to_write = velocities_checkpoint.copy()
            current_sim_time = time_checkpoint
            # Need to re-write the checkpointed data if the advance was reverted
            participant.write_data(mesh_name, "Velocity", vertex_ids_precice, velocities_to_write)


        # Get the next suggested coupling time step size for the *next* iteration
        if participant.is_coupling_ongoing(): # Check again before getting max time step
             dt_precice = participant.get_max_time_step_size()
        
    participant.finalize()
    print("preCICE coupling finished.")

if __name__ == "__main__":
    main()