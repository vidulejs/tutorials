import numpy as np
import precice
import time
import sys

def main():
    precice_config = '../precice-config.xml'
    solver_process_index = 0
    solver_process_size = 1

    def create_velocity_field(x, y):
        u = 4.0 * y * (2.0 - y)
        v = np.zeros_like(u)
        return u, v

    x_min, x_max = 0.0, 6.0
    y_min, y_max = 0.0, 2.0
    nx, ny = 60, 20
    
    x_centers = np.linspace(x_min + (x_max - x_min)/(2*nx), x_max - (x_max - x_min)/(2*nx), nx)
    y_centers = np.linspace(y_min + (y_max - y_min)/(2*ny), y_max - (y_max - y_min)/(2*ny), ny)
    
    X, Y = np.meshgrid(x_centers, y_centers)
    
    coordinates = np.vstack([X.flatten(), Y.flatten()]).T
    
    U, V = create_velocity_field(X, Y)
    velocities = np.column_stack([U.flatten(), V.flatten()])
    
    n_vertices = coordinates.shape[0]
    
    participant_name = "Fluid"
    mesh_name = "Fluid-Mesh"
    participant = precice.Participant(participant_name, precice_config, solver_process_index, solver_process_size)
    
    vertex_ids = participant.set_mesh_vertices(mesh_name, coordinates)
    
    dt = participant.initialize()
    time = 0
    
    while participant.is_coupling_ongoing():
        dt = participant.get_max_time_step_size()
        
        print(f"Advancing time {time + dt:.5f}...")
        
        if participant.requires_writing_checkpoint():
            velocities_checkpoint = velocities.copy()
            time_checkpoint = time
        
        participant.write_data(mesh_name, "Velocity", vertex_ids, velocities)
        
        participant.advance(dt)
        
        time += dt
        
        if participant.requires_reading_checkpoint():
            velocities = velocities_checkpoint.copy()
            time = time_checkpoint
        
    participant.finalize()
    print("Finished.")

if __name__ == "__main__":
    main()