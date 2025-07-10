import numpy as np
from torch.utils.data import Dataset

class GriddedDataset(Dataset):
    """PyTorch Dataset for gridded simulation data."""
    def __init__(self, data_np, initial_step, unroll_step):
        self.data = data_np
        self.initial_step = initial_step
        self.unroll_step = unroll_step
        self.num_time_steps = data_np.shape[2]

    def __len__(self):
        return self.num_time_steps - self.initial_step - self.unroll_step

    def __getitem__(self, idx):
        t_start = idx
        t_end_input = t_start + self.initial_step
        input_slice = self.data[:, :, t_start:t_end_input, :]
        
        t_end_target = t_end_input + self.unroll_step
        target_sequence = self.data[:, :, t_end_input:t_end_target, :]

        # Reshape using numpy and return numpy arrays
        inp = np.transpose(input_slice, (2, 3, 0, 1)).reshape(-1, self.data.shape[0], self.data.shape[1])
        tar = np.transpose(target_sequence, (2, 3, 0, 1))
        
        return inp, tar
