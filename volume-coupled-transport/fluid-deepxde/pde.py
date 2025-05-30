import deepxde as dde
import numpy as np
import matplotlib.pyplot as plt

# --- 1. Parameters ---
# Fluid properties
nu = 0.0002  # Kinematic viscosity (m^2/s)
# rho = 1.0 # Density (kg/m^3)

# Inlet velocity
inlet_u_velocity = 2.0  # m/s

# Domain dimensions
domain_length = 10.0  # m
domain_width = 2.0   # m
time_final = 5.0     # s

# Rectangular obstacle dimensions and position
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
# Spatio-temporal domain
geom_channel = dde.geometry.Rectangle([0, 0], [domain_length, domain_width])
geom_obstacle = dde.geometry.Rectangle([obstacle_x0, obstacle_y0], [obstacle_x1, obstacle_y1])
geom_spatial = dde.geometry.CSGDifference(geom_channel, geom_obstacle)

timedomain = dde.geometry.TimeDomain(0, time_final)
geomtime = dde.geometry.GeometryXTime(geom_spatial, timedomain)

# --- 4. Define Initial and Boundary Conditions ---

# Initial conditions: u=0, v=0, p=0 (kinematic pressure) at t=0
def initial_condition_zeros(x_coords_t):
    return np.zeros((x_coords_t.shape[0], 1))

ic_u = dde.icbc.IC(geomtime, initial_condition_zeros, lambda _, on_initial: on_initial, component=0)
ic_v = dde.icbc.IC(geomtime, initial_condition_zeros, lambda _, on_initial: on_initial, component=1)
ic_p = dde.icbc.IC(geomtime, initial_condition_zeros, lambda _, on_initial: on_initial, component=2) # Kinematic pressure initial ref 0

# Boundary condition predicates
def boundary_inlet(x_coords_t, on_boundary):
    return on_boundary and np.isclose(x_coords_t[0], 0)

def boundary_outlet(x_coords_t, on_boundary):
    return on_boundary and np.isclose(x_coords_t[0], domain_length)

def boundary_solid_walls(x_coords_t, on_boundary): # Covers channel top/bottom and obstacle
    is_on_inlet = np.isclose(x_coords_t[0], 0)
    is_on_outlet = np.isclose(x_coords_t[0], domain_length)
    return on_boundary and (not is_on_inlet) and (not is_on_outlet)

# Inlet BCs
bc_inlet_u = dde.icbc.DirichletBC(
    geomtime,
    lambda x: inlet_u_velocity * np.ones((x.shape[0], 1)),
    boundary_inlet,
    component=0
)
bc_inlet_v = dde.icbc.DirichletBC(
    geomtime,
    lambda x: np.zeros((x.shape[0], 1)),
    boundary_inlet,
    component=1
)
bc_inlet_p_neumann = dde.icbc.NeumannBC( # p: zeroGradient
    geomtime,
    lambda x: np.zeros((x.shape[0], 1)),
    boundary_inlet,
    component=2
)

# Outlet BCs
bc_outlet_p_dirichlet = dde.icbc.DirichletBC( # p: fixedValue (kinematic reference 0)
    geomtime,
    lambda x: np.zeros((x.shape[0], 1)),
    boundary_outlet,
    component=2
)
bc_outlet_u_neumann = dde.icbc.NeumannBC( # u: zeroGradient
    geomtime,
    lambda x: np.zeros((x.shape[0], 1)),
    boundary_outlet,
    component=0
)
bc_outlet_v_neumann = dde.icbc.NeumannBC( # v: zeroGradient
    geomtime,
    lambda x: np.zeros((x.shape[0], 1)),
    boundary_outlet,
    component=1
)

# Solid Walls (Channel top/bottom, Obstacle) BCs
bc_walls_u_noslip = dde.icbc.DirichletBC( # u: noSlip
    geomtime,
    lambda x: np.zeros((x.shape[0], 1)),
    boundary_solid_walls,
    component=0
)
bc_walls_v_noslip = dde.icbc.DirichletBC( # v: noSlip
    geomtime,
    lambda x: np.zeros((x.shape[0], 1)),
    boundary_solid_walls,
    component=1
)
bc_walls_p_neumann = dde.icbc.NeumannBC( # p: zeroGradient
    geomtime,
    lambda x: np.zeros((x.shape[0], 1)),
    boundary_solid_walls,
    component=2
)

# Consolidate BCs and ICs
# Order matters for loss_weights later
bcs_and_ics = [
    ic_u, ic_v, ic_p,                       # Initial conditions
    bc_inlet_u, bc_inlet_v, bc_inlet_p_neumann, # Inlet
    bc_outlet_u_neumann, bc_outlet_v_neumann, bc_outlet_p_dirichlet, # Outlet
    bc_walls_u_noslip, bc_walls_v_noslip, bc_walls_p_neumann # Solid walls
]

# --- 5. Create the PDE Problem ---
data = dde.data.TimePDE(
    geomtime,
    navier_stokes_pde,
    bcs_and_ics,
    num_domain=20000,    # Collocation points for PDE residual in the domain
    num_boundary=4000,   # Points for boundary conditions
    num_initial=4000,    # Points for initial conditions
    # solution=None,     # No analytical solution provided for training
    # num_test=5000      # Optional: for evaluating PDE residuals on a test set during/after training
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

    # Loss weights: [PDE_cont, PDE_xm, PDE_ym, BC_iu, BC_iv, BC_op, BC_nsu, BC_nsv, IC_u, IC_v, IC_p]
    # Order corresponds to definition in navier_stokes_pde and bcs_and_ics list
    loss_weights = [
        2, 1, 1,          # PDE residuals (continuity, x-mom, y-mom)
        100, 100,          # Inlet u, v
        100,               # Outlet p
        100, 100,          # No-slip u, v
        100, 100, 100     # IC u, v, p
    ]
    # Compile the model
    model.compile("adam", lr=1e-3, loss_weights=loss_weights)

    # # Train the model
    # # For a real run, iterations should be much higher (e.g., 50000-200000+)
    # # and consider a learning rate scheduler.
    # losshistory, train_state = model.train(iterations=20000, display_every=1000) # Initial coarse training

    # Refine with a smaller learning rate
    model.compile("adam", lr=1e-4, loss_weights=loss_weights)
    losshistory_2, train_state_2 = model.train(iterations=10000, display_every=200)
    
    # Combine loss histories for plotting
    # This part needs careful handling if you want to combine them properly
    # For simplicity, we'll just plot the last one or save both separately.
    # To properly combine:
    # losshistory.loss_train = np.vstack((losshistory.loss_train, losshistory_2.loss_train))
    # losshistory.loss_test = np.vstack((losshistory.loss_test, losshistory_2.loss_test)) # if num_test is used
    # losshistory.metrics_test = np.vstack((losshistory.metrics_test, losshistory_2.metrics_test)) # if num_test is used
    # train_state.best_step = train_state_2.best_step # or sum iterations
    # train_state.best_loss_train = train_state_2.best_loss_train
    # etc.

    # Save the model (optional)
    model.save("navier_stokes_channel_obstacle_model")

    # Plot loss history
    # dde.saveplot(losshistory, train_state, issave=True, isplot=True, output_dir="NS_forward_results_1")
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