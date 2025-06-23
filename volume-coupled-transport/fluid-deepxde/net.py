import deepxde as dde
import torch
import torch.nn as nn


# Define the dimensions of the neural network
# [input_dim] + [hidden_dim] * num_hidden_layers + [output_dim]
input_dim = 3   # (x, y, t)
output_dim = 3  # (u, v, p)
hidden_dim = 32
num_layers = 3

activation = "tanh"
initializer = "Glorot uniform"

num_domain_points = 20000
num_boundary_points = 5000
num_initial_points = 4000
num_test_points = None

training_output_dir = "train_dir/"
model_base_name = "ns_obstacle"

layer_size = [input_dim] + [hidden_dim] * num_layers + [hidden_dim*2] * 2 + [hidden_dim] * num_layers + [output_dim]

net = dde.nn.FNN(
	layer_size,
	activation=activation,
	kernel_initializer=initializer)