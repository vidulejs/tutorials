import deepxde as dde

layer_size = [3] + [64] * 6 + [3]
activation = "tanh"
initializer = "Glorot uniform"

num_domain_points = 20000
num_boundary_points = 5000
num_initial_points = 4000
num_test_points = 5000

iterations = 600

lr = 1e-4

training_output_dir = "train_dir"
model_base_name = "ns_obstacle"

def define_pde(geom, pde_residuals, loss_terms):
    data = dde.data.TimePDE(
        geom,
        pde_residuals,
        loss_terms,
        num_domain=num_domain_points,
        num_boundary=num_boundary_points,
        num_initial=num_initial_points,
        num_test=num_test_points
    )
    return data

net = dde.nn.FNN(layer_size, activation, initializer)

# model = dde.Model(data, net)

# model.compile(
#     "adam",
#     lr=lr,
#     loss_weights=loss_weights
# )

# # --- 7. Training Callbacks ---
# # Model checkpoint callback
# checkpoint_callback = dde.callbacks.ModelCheckpoint(
#     os.path.join(training_output_dir, f"{model_base_name}"),
#     save_better_only=True,
#     period=50
# )

# losshistory, train_state = model.train(
#     iterations=iterations,
#     callbacks=[checkpoint_callback]
# )

# # --- 9. Save the Trained Model ---
# model_path = os.path.join(training_output_dir, f"{model_base_name}_final")
# model.save(model_path)
# print(f"Model saved to: {model_path}")

# # --- 10. Plot Training History ---
# dde.saveplot(
#     losshistory, 
#     train_state, 
#     issave=True, 
#     isplot=True,
#     output_dir=training_output_dir
# )


            # x = geom.uniform_points(num, boundary=True)
            # y_true = ...
            # y_pred = model.predict(x)
            # error= dde.metrics.l2_relative_error(y_true, y_pred)