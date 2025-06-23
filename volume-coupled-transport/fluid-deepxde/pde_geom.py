import deepxde as dde
import numpy as np

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

geom_channel = dde.geometry.Rectangle([0, 0], [domain_length, domain_width])
geom_obstacle = dde.geometry.Rectangle([obstacle_x0, obstacle_y0], [obstacle_x1, obstacle_y1])
geom_spatial = dde.geometry.CSGDifference(geom_channel, geom_obstacle)
timedomain = dde.geometry.TimeDomain(0, time_final)
geom = dde.geometry.GeometryXTime(geom_spatial, timedomain)

def initial_condition_zeros(x_coords_t):
    return np.zeros((x_coords_t.shape[0], 1))

ic_u = dde.icbc.IC(geom, initial_condition_zeros, lambda _, on_initial: on_initial, component=0)
ic_v = dde.icbc.IC(geom, initial_condition_zeros, lambda _, on_initial: on_initial, component=1)
ic_p = dde.icbc.IC(geom, initial_condition_zeros, lambda _, on_initial: on_initial, component=2)

def boundary_inlet(x_coords_t, on_boundary):
    return on_boundary and np.isclose(x_coords_t[0], 0)

def boundary_outlet(x_coords_t, on_boundary):
    return on_boundary and np.isclose(x_coords_t[0], domain_length)

def boundary_solid_walls(x_coords_t, on_boundary):
    is_on_inlet = np.isclose(x_coords_t[0], 0)
    is_on_outlet = np.isclose(x_coords_t[0], domain_length)
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
bc_outlet_p_dirichlet = dde.icbc.DirichletBC(
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

unsupervised_loss_terms = [
    ic_p, ic_u, ic_v,
    bc_inlet_p_dirichlet, bc_inlet_u_neumann, bc_inlet_v_neumann,
    bc_outlet_p_dirichlet, bc_outlet_u_neumann, bc_outlet_v_neumann,
    bc_walls_p_neumann, bc_walls_u_noslip, bc_walls_v_noslip
]

unsupervised_loss_weights = [
    100, 50, 50,      # PDE residuals: continuity, x-momentum, y-momentum
    1, 1, 1,    # ICs: p, u, v
    10, 10, 10,   # Inlet BCs: p(D), u(N), v(N)
    10, 10, 10,   # Outlet BCs: p(D), u(N), v(N)
    10, 10, 10    # Wall BCs: p(N), u(D), v(D)
]