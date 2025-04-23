import deepxde as dde
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import os

# import model setup
from pde import *

time_points = np.linspace(0, 10, 21)

model = dde.Model(data, net)
model.compile("adam", lr=1e-3)

# Load the trained model
model_path = "burgers2d_model-20000.pt"
if not os.path.exists(model_path):
    print("Wrong model path")
    exit()

model.restore(save_path=model_path, verbose=1)

# Generate random points in the domain for visualization
x_geom = geom.random_points(1000)
plt.figure(figsize=(8, 6))
plt.scatter(x_geom[:, 0], x_geom[:, 1], s=1)
plt.title("Random sampling points in the domain")
plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")
plt.savefig("domain_points.png")
plt.close()

def plot_at_time(t_sample, n_points=100, save_path=None):
    print(f"Visualizing at t = {t_sample}...")
    
    # Create a regular grid for visualization
    x = np.linspace(0, domain_length, n_points)
    y = np.linspace(0, domain_width, n_points)
    X, Y = np.meshgrid(x, y)
    X_flat = X.flatten()
    Y_flat = Y.flatten()
    T_flat = np.ones_like(X_flat) * t_sample
    X_star = np.vstack((X_flat, Y_flat, T_flat)).T
    
    # Predict velocities
    u_pred = model.predict(X_star)
    
    # Reshape for plotting
    u_star = u_pred[:, 0].reshape(n_points, n_points)
    v_star = u_pred[:, 1].reshape(n_points, n_points)
    
    # Calculate velocity magnitude for coloring
    magnitude = np.sqrt(u_star**2 + v_star**2)
    
    # Mask the obstacle region
    mask = np.zeros_like(u_star, dtype=bool)
    for i in range(n_points):
        for j in range(n_points):
            if (obstacle_x0 <= X[i, j] <= obstacle_x1 and 
                obstacle_y0 <= Y[i, j] <= obstacle_y1):
                mask[i, j] = True
    
    # Apply mask
    u_star_masked = np.ma.array(u_star, mask=mask)
    v_star_masked = np.ma.array(v_star, mask=mask)
    magnitude_masked = np.ma.array(magnitude, mask=mask)
    
    # Plot the solution
    fig, axs = plt.subplots(3, 1, figsize=(12, 18))
    
    # U velocity component
    im1 = axs[0].pcolor(X, Y, u_star_masked, cmap="RdBu_r", vmin=-1, vmax=5)
    plt.colorbar(im1, ax=axs[0])
    axs[0].set_title(f"u(x,y,t={t_sample})")
    axs[0].set_xlabel("x")
    axs[0].set_ylabel("y")
    axs[0].set_aspect('equal')
    
    # V velocity component
    im2 = axs[1].pcolor(X, Y, v_star_masked, cmap="RdBu_r", vmin=-1, vmax=1)
    plt.colorbar(im2, ax=axs[1])
    axs[1].set_title(f"v(x,y,t={t_sample})")
    axs[1].set_xlabel("x")
    axs[1].set_ylabel("y")
    axs[1].set_aspect('equal')
    
    # Velocity magnitude
    im3 = axs[2].pcolor(X, Y, magnitude_masked, cmap="RdBu_r", vmin=0, vmax=6)
    plt.colorbar(im3, ax=axs[2])
    axs[2].set_title(f"Velocity magnitude at t={t_sample}")
    axs[2].set_xlabel("x")
    axs[2].set_ylabel("y")
    axs[2].set_aspect('equal')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path)
        print(f"Figure saved to {save_path}")
    
    return fig

os.makedirs("viz", exist_ok=True)

for t in time_points:
    fig = plot_at_time(t, n_points=100, save_path=f"viz/burgers2d_t{t:.1f}.png")
    plt.close(fig)

print("Visualization complete!")
