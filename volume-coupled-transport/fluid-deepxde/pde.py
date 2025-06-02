import deepxde as dde
import numpy as np
import matplotlib.pyplot as plt
import os

# --- 1. Parameters ---
nu = 0.0002
inlet_kinematic_pressure = 1.0
outlet_kinematic_pressure = 0.0

domain_length = 6.0
domain_width = 2.0
time_final = 10.0

obstacle_x0 = 2.0
obstacle_x1 = 3.0
obstacle_y0 = 0.0
obstacle_y1 = 1.0

# --- 2. Define the Navier-Stokes Equations ---
def navier_stokes_pde(x_coords, y_solution):
    """
    2D Incompressible Navier-Stokes equations (kinematic pressure formulation):
    Continuity: du/dx + dv/dy = 0
    Momentum_x: du/dt + u*du/dx + v*du/dy + dp/dx - nu*(d^2u/dx^2 + d^2u/dy^2) = 0
    Momentum_y: dv/dt + u*dv/dx + v*dv/dy + dp/dy - nu*(d^2v/dx^2 + d^2v/dy^2) = 0

    Inputs:
        x_coords: (x, y, t) tensor for coordinates
        y_solution: (u, v, p) tensor for velocity components and kinematic pressure
    """
    u = y_solution[:, 0:1]
    v = y_solution[:, 1:2]
    p = y_solution[:, 2:3]

    # First derivatives
    du_x = dde.grad.jacobian(y_solution, x_coords, i=0, j=0)
    du_y = dde.grad.jacobian(y_solution, x_coords, i=0, j=1)
    du_t = dde.grad.jacobian(y_solution, x_coords, i=0, j=2)

    dv_x = dde.grad.jacobian(y_solution, x_coords, i=1, j=0)
    dv_y = dde.grad.jacobian(y_solution, x_coords, i=1, j=1)
    dv_t = dde.grad.jacobian(y_solution, x_coords, i=1, j=2)

    dp_x = dde.grad.jacobian(y_solution, x_coords, i=2, j=0)
    dp_y = dde.grad.jacobian(y_solution, x_coords, i=2, j=1)

    # Second derivatives (Laplacians)
    du_xx = dde.grad.hessian(y_solution, x_coords, component=0, i=0, j=0)
    du_yy = dde.grad.hessian(y_solution, x_coords, component=0, i=1, j=1)

    dv_xx = dde.grad.hessian(y_solution, x_coords, component=1, i=0, j=0)
    dv_yy = dde.grad.hessian(y_solution, x_coords, component=1, i=1, j=1)

    # Continuity equation
    continuity_eq = du_x + dv_y
    # X-momentum equation
    x_momentum_eq = du_t + (u * du_x + v * du_y) + dp_x - nu * (du_xx + du_yy)
    # Y-momentum equation
    y_momentum_eq = dv_t + (u * dv_x + v * dv_y) + dp_y - nu * (dv_xx + dv_yy)

    return [continuity_eq, x_momentum_eq, y_momentum_eq]

# --- 3. Define the Domain ---

geom_channel = dde.geometry.Rectangle([0, 0], [domain_length, domain_width]) # domain_length is now 6.0
geom_obstacle = dde.geometry.Rectangle([obstacle_x0, obstacle_y0], [obstacle_x1, obstacle_y1])
geom_spatial = dde.geometry.CSGDifference(geom_channel, geom_obstacle)
timedomain = dde.geometry.TimeDomain(0, time_final)
# Spatio-temporal domain
geom = dde.geometry.GeometryXTime(geom_spatial, timedomain)

# --- 4. Define Initial and Boundary Conditions ---

def initial_condition_zeros(x_coords_t):
    return np.zeros((x_coords_t.shape[0], 1))
ic_u = dde.icbc.IC(geom, initial_condition_zeros, lambda _, on_initial: on_initial, component=0)
ic_v = dde.icbc.IC(geom, initial_condition_zeros, lambda _, on_initial: on_initial, component=1)
ic_p = dde.icbc.IC(geom, initial_condition_zeros, lambda _, on_initial: on_initial, component=2)

def boundary_inlet(x_coords_t, on_boundary):
    return on_boundary and np.isclose(x_coords_t[0], 0)
def boundary_outlet(x_coords_t, on_boundary):
    return on_boundary and np.isclose(x_coords_t[0], domain_length) # domain_length is now 6.0
def boundary_solid_walls(x_coords_t, on_boundary):
    is_on_inlet = np.isclose(x_coords_t[0], 0)
    is_on_outlet = np.isclose(x_coords_t[0], domain_length) # domain_length is now 6.0
    return on_boundary and (not is_on_inlet) and (not is_on_outlet)

bc_inlet_p_dirichlet = dde.icbc.DirichletBC(
    geom, lambda x: inlet_kinematic_pressure * np.ones((x.shape[0], 1)), boundary_inlet, component=2
)
bc_inlet_u_neumann = dde.icbc.NeumannBC(
    geom, lambda x: np.zeros((x.shape[0], 1)), boundary_inlet, component=0
)
bc_inlet_v_neumann = dde.icbc.NeumannBC(
    geom, lambda x: np.zeros((x.shape[0], 1)), boundary_inlet, component=1
)
bc_outlet_p_dirichlet_new = dde.icbc.DirichletBC(
    geom, lambda x: outlet_kinematic_pressure * np.ones((x.shape[0], 1)), boundary_outlet, component=2
)
bc_outlet_u_neumann = dde.icbc.NeumannBC(
    geom, lambda x: np.zeros((x.shape[0], 1)), boundary_outlet, component=0
)
bc_outlet_v_neumann = dde.icbc.NeumannBC(
    geom, lambda x: np.zeros((x.shape[0], 1)), boundary_outlet, component=1
)
bc_walls_u_noslip = dde.icbc.DirichletBC(
    geom, lambda x: np.zeros((x.shape[0], 1)), boundary_solid_walls, component=0
)
bc_walls_v_noslip = dde.icbc.DirichletBC(
    geom, lambda x: np.zeros((x.shape[0], 1)), boundary_solid_walls, component=1
)
bc_walls_p_neumann = dde.icbc.NeumannBC(
    geom, lambda x: np.zeros((x.shape[0], 1)), boundary_solid_walls, component=2
)
bcs_and_ics_pressure_driven = [
    ic_u, ic_v, ic_p,
    bc_inlet_p_dirichlet, bc_inlet_u_neumann, bc_inlet_v_neumann,
    bc_outlet_p_dirichlet_new, bc_outlet_u_neumann, bc_outlet_v_neumann,
    bc_walls_u_noslip, bc_walls_v_noslip, bc_walls_p_neumann
]

# --- 5. Create the PDE Problem with NEW BCs and UPDATED GEOMETRY ---
data = dde.data.TimePDE(
    geom, # This now uses the updated domain_length
    navier_stokes_pde,
    bcs_and_ics_pressure_driven,
    num_domain=20000, # You might want to adjust these based on the new domain sizeR
    num_boundary=5000,
    num_initial=4000,
)

# --- 6. Neural Network Architecture ---
# Input: (x, y, t) -> 3 features
# Output: (u, v, p) -> 3 features
layer_size = [3] + [64] * 6 + [3]  # Example: 3 input, 6 hidden layers of 64 neurons, 3 output
activation = "tanh"
initializer = "Glorot uniform" # or "Glorot normal"
net = dde.nn.FNN(layer_size, activation, initializer)

# --- 7. Model Training ---
if __name__ == "__main__":
    # Create the model
    model = dde.Model(data, net)

    # Loss weights:
    # PDE: [cont, x-mom, y-mom]
    # ICs: [ic_u, ic_v, ic_p]
    # Inlet: [bc_iu, bc_iv, bc_ip_N]
    # Outlet: [bc_ou_N, bc_ov_N, bc_op_D]
    # Walls: [bc_wu_D, bc_wv_D, bc_wp_N]
    # Total = 3 (PDEs) + 3 (ICs) + 3 (Inlet) + 3 (Outlet) + 3 (Walls) = 15 loss terms
    loss_weights = [
        10, 5, 5,              # PDE residuals
        10, 10, 100,           # ICs u, v, p
        10, 10, 100,          # Inlet u(D), v(D), p(N)
        10, 10, 100,          # Outlet u(N), v(N), p(D)
        10, 10, 100           # Walls u(D), v(D), p(N)
    ]

    previous_model_checkpoint_path = "navier_stokes_channel_obstacle_model-10000.pt"

    if not os.path.exists(previous_model_checkpoint_path):
        print(f"Error: Previous model checkpoint {previous_model_checkpoint_path} not found!")
        print("Please ensure the base model was trained on the correct (old) geometry or correct the path.")
        print("If the old model was for a different domain size, fine-tuning might be challenging.")
        exit()
        print("Proceeding without loading pre-trained weights (training from scratch).")
        model_loaded = False
    else:
        print(f"Loading weights from: {previous_model_checkpoint_path}")
        try:
            model.restore(previous_model_checkpoint_path, verbose=1)
            print("Model weights restored.")
            model_loaded = True
        except Exception as e:
            print(f"Could not restore model: {e}")
            print("Proceeding without loading pre-trained weights (training from scratch).")
            model_loaded = False

    model.compile("adam", lr=1e-5, loss_weights=loss_weights)
    losshistory_2, train_state_2 = model.train(iterations=600, display_every=200)

    # Save the model
    model.save("navier_stokes_channel_obstacle_model")

    # Plot loss history
    dde.saveplot(losshistory_2, train_state_2, issave=True, isplot=True, output_dir="NS_forward_results")

    print("Training finished. Model and plots saved in 'NS_forward_results' directory.")

    # --- 8. Post-processing and Visualization (Example) ---
    print("\nGenerating sample plot at t=time_final/2...")
    t_plot = time_final / 2.0
    nx, ny = 100, 50 # Resolution for plotting grid
    x_plot_coords = np.linspace(0, domain_length, nx)
    y_plot_coords = np.linspace(0, domain_width, ny)
    X_plot, Y_plot = np.meshgrid(x_plot_coords, y_plot_coords)
    
    # Create a mask for the obstacle region
    obstacle_mask = (X_plot >= obstacle_x0) & (X_plot <= obstacle_x1) & \
                    (Y_plot >= obstacle_y0) & (Y_plot <= obstacle_y1)

    x_flat = X_plot.flatten()
    y_flat = Y_plot.flatten()
    t_flat = np.full_like(x_flat, t_plot)
    
    xyt_topred = np.vstack((x_flat, y_flat, t_flat)).T
    uvp_pred = model.predict(xyt_topred)
    u_pred = uvp_pred[:, 0].reshape(X_plot.shape)
    v_pred = uvp_pred[:, 1].reshape(X_plot.shape)
    p_pred = uvp_pred[:, 2].reshape(X_plot.shape)

    # Apply mask to predictions (set values inside obstacle to NaN for plotting)
    u_pred[obstacle_mask] = np.nan
    v_pred[obstacle_mask] = np.nan
    p_pred[obstacle_mask] = np.nan

    fig, axes = plt.subplots(3, 1, figsize=(10, 12))
    
    # Plot u
    cp0 = axes[0].contourf(X_plot, Y_plot, u_pred, levels=50, cmap="viridis")
    fig.colorbar(cp0, ax=axes[0])
    axes[0].add_patch(plt.Rectangle((obstacle_x0, obstacle_y0), obstacle_x1-obstacle_x0, obstacle_y1-obstacle_y0, facecolor='gray'))
    axes[0].set_title(f"Velocity u (m/s) at t={t_plot:.2f}s")
    axes[0].set_xlabel("x (m)")
    axes[0].set_ylabel("y (m)")
    axes[0].set_aspect('equal')

    # Plot v
    cp1 = axes[1].contourf(X_plot, Y_plot, v_pred, levels=50, cmap="viridis")
    fig.colorbar(cp1, ax=axes[1])
    axes[1].add_patch(plt.Rectangle((obstacle_x0, obstacle_y0), obstacle_x1-obstacle_x0, obstacle_y1-obstacle_y0, facecolor='gray'))
    axes[1].set_title(f"Velocity v (m/s) at t={t_plot:.2f}s")
    axes[1].set_xlabel("x (m)")
    axes[1].set_ylabel("y (m)")
    axes[1].set_aspect('equal')

    # Plot p
    cp2 = axes[2].contourf(X_plot, Y_plot, p_pred, levels=50, cmap="viridis")
    fig.colorbar(cp2, ax=axes[2])
    axes[2].add_patch(plt.Rectangle((obstacle_x0, obstacle_y0), obstacle_x1-obstacle_x0, obstacle_y1-obstacle_y0, facecolor='gray'))
    axes[2].set_title(f"Kinematic Pressure p (m^2/s^2) at t={t_plot:.2f}s")
    axes[2].set_xlabel("x (m)")
    axes[2].set_ylabel("y (m)")
    axes[2].set_aspect('equal')

    plt.tight_layout()
    plt.savefig("NS_forward_results/flow_field_snapshot.png")
    plt.show()