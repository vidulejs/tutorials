import deepxde as dde
import numpy as np
import matplotlib.pyplot as plt

# Parameters
nu = 0.0002  # Viscosity
inlet_velocity = 2.0
domain_length = 6.0
domain_width = 2.0
time_final = 3.0

obstacle_x0 = 2.0
obstacle_x1 = 3.0
obstacle_y0 = 0.0
obstacle_y1 = 1.0

# Define the 2D Burgers equation
def pde(x, y):
    """
    2D Burgers equation:
    du/dt + u*du/dx + v*du/dy = nu*(d²u/dx² + d²u/dy²)
    dv/dt + u*dv/dx + v*dv/dy = nu*(d²v/dx² + d²v/dy²)
    
    where u and v are the velocity components in x and y directions.
    
    Input:
        x: (x, y, t)
        y: (u, v)
    """

    u, v = y[:, 0:1], y[:, 1:2]
    
    # Get the derivatives
    du_x = dde.grad.jacobian(y, x, i=0, j=0)
    du_y = dde.grad.jacobian(y, x, i=0, j=1)
    du_t = dde.grad.jacobian(y, x, i=0, j=2)
    
    dv_x = dde.grad.jacobian(y, x, i=1, j=0)
    dv_y = dde.grad.jacobian(y, x, i=1, j=1)
    dv_t = dde.grad.jacobian(y, x, i=1, j=2)
    
    du_xx = dde.grad.hessian(y, x, component=0, i=0, j=0)
    du_yy = dde.grad.hessian(y, x, component=0, i=1, j=1)
    
    dv_xx = dde.grad.hessian(y, x, component=1, i=0, j=0)
    dv_yy = dde.grad.hessian(y, x, component=1, i=1, j=1)
    
    # Burgers equations
    eq_u = du_t + u * du_x + v * du_y - nu * (du_xx + du_yy)
    eq_v = dv_t + u * dv_x + v * dv_y - nu * (dv_xx + dv_yy)
    
    return [eq_u, eq_v]

# Define the domain
geom_large = dde.geometry.Rectangle([0, 0], [domain_length, domain_width])
geom_obstacle = dde.geometry.Rectangle([obstacle_x0, obstacle_y0], [obstacle_x1, obstacle_y1])
geom = dde.geometry.CSGDifference(geom_large, geom_obstacle)

timedomain = dde.geometry.TimeDomain(0, time_final)
geomtime = dde.geometry.GeometryXTime(geom, timedomain)

# zero initial condition
def initial_condition(x): 
    return np.zeros_like(x[:, 0:1])  # u and v are both zero at t=0

# inlet boundary (x=0)
def inlet_boundary(x, on_boundary):
    return on_boundary and np.isclose(x[0], 0)

# all other boundaries
def wall_boundary(x, on_boundary):
    return on_boundary and not np.isclose(x[0], 0)


# Create boundary and initial conditions
# Inlet velocity u_x = 1
bc_inlet_x = dde.icbc.DirichletBC(geomtime, lambda y: inlet_velocity * np.ones((len(y), 1)), inlet_boundary, component=0)
bc_inlet_y = dde.icbc.DirichletBC(geomtime, lambda y: np.zeros((len(y), 1)), inlet_boundary, component=1)

bc_walls_x = dde.icbc.DirichletBC(geomtime, lambda y: np.zeros((len(y), 1)), wall_boundary, component=0)
bc_walls_y = dde.icbc.DirichletBC(geomtime, lambda y: np.zeros((len(y), 1)), wall_boundary, component=1)

ic = dde.icbc.IC(geomtime, initial_condition, lambda _, on_initial: on_initial)

# Create the PDE problem
data = dde.data.TimePDE(
    geomtime, 
    pde, 
    [bc_inlet_x, bc_inlet_y, bc_walls_x, bc_walls_y, ic], 
    num_domain=5000,  # No. points in the domain
    num_boundary=1000,  # No. on the boundary
    num_initial=1000,   # No. points at the initial time
    solution=None,  # Analytical solution
    num_test=5000   # No. test points
)

# Neural network architecture
layer_size = [3] + [50] * 4 + [2]  # Input: (x, y, t), Output: (u, v)
activation = "tanh"
initialization = "Glorot normal"
net = dde.nn.FNN(layer_size, activation, initialization)

if __name__ == "__main__":
    # Create the model
    model = dde.Model(data, net)

    # Weights for boundary conditions
    # loss_weights = [1, 1, 2, 2, 0.5]  # Weights for [PDE_u, PDE_v, bc_inlet, bc_walls, ic]

    # Compile the model
    model.compile("adam", lr=1e-3)

    # Train the model
    losshistory, train_state = model.train(iterations=5000)

    # Save the model
    model.save("burgers2d_model")

    dde.saveplot(losshistory, train_state, issave=True, isplot=True)

