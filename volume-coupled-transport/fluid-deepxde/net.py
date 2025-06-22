import deepxde as dde

layer_size = [3] + [64] * 6 + [3]
activation = "tanh"
initializer = "Glorot uniform"

num_domain_points = 20000
num_boundary_points = 5000
num_initial_points = 4000
num_test_points = 5000

training_output_dir = "train_dir"
model_base_name = "ns_obstacle"

net = dde.nn.FNN(layer_size, activation, initializer)