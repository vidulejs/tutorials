# --- Path Configuration ---
PRETRAINED_MODEL_PATH = 'models/2D_CFD_Turb_M0.1_Eta1e-08_Zeta1e-08_periodic_512_Train_Unet-PF-20.pt'
FINETUNED_MODEL_DIR = 'finetuned_model'
FINETUNED_MODEL_NAME = 'finetuned_model_final_20250705_224038.pt'
GRIDDED_DATA_PATH = 'training_data_gridded.npz'
OUTPUTS_DIR = 'outputs'

# --- Training Hyperparameters ---
EPOCHS = 100
LEARNING_RATE = 1e-3
BATCH_SIZE = 8
INITIAL_STEP = 10
UNROLL_STEP = 5

# --- Model & Simulation Parameters ---
MODEL_RESOLUTION = (256, 256)
TOTAL_SIM_TIME = 10.0
TOTAL_SIM_STEPS = 2000
DT = TOTAL_SIM_TIME / TOTAL_SIM_STEPS

# --- Inference & Visualization Parameters ---
INFERENCE_STEPS = 200
ANIMATION_FPS = 30

# --- Physical Domain (for plotting) ---
DOMAIN_LENGTH = 6.0
DOMAIN_WIDTH = 2.0
OBSTACLE_X0 = 2.0
OBSTACLE_X1 = 3.0
OBSTACLE_Y0 = 0.0
OBSTACLE_Y1 = 1.0
