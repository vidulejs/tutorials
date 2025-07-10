import torch
print(f"Initial Check - CUDA available: {torch.cuda.is_available()}")

import numpy as np
import os
from datetime import datetime
from torch import nn
from torch.utils.data import DataLoader
import torch.nn.functional as F

from src.config import *
from src.data import GriddedDataset
from src.model import UNet2d

def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    os.makedirs(FINETUNED_MODEL_DIR, exist_ok=True)
    os.makedirs(OUTPUTS_DIR, exist_ok=True)
    log_path = os.path.join(OUTPUTS_DIR, f'training_log_{timestamp}.txt')
    
    with open(log_path, 'w') as f:
        f.write("Epoch,Average Loss\n")

    print("Initializing UNet model from scratch...")
    
    in_channels = 3 * INITIAL_STEP
    out_channels = 3
    
    print(f"Model parameters: in_channels={in_channels}, out_channels={out_channels}")

    model = UNet2d(in_channels=in_channels, out_channels=out_channels).to(device)
    model.train()
    
    print("Model initialized successfully from scratch.")

    print(f"Loading and preparing gridded data from '{GRIDDED_DATA_PATH}'...")
    full_data_np = np.load(GRIDDED_DATA_PATH)['data'][0]
    
    train_dataset = GriddedDataset(full_data_np, INITIAL_STEP, UNROLL_STEP)
    train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True, num_workers=0)

    print(f"Starting training for {EPOCHS} epoch(s) with unroll_step={UNROLL_STEP}...")
    optimizer = torch.optim.Adam(model.parameters(), lr=LEARNING_RATE)
    loss_fn = nn.MSELoss(reduction="mean")

    best_loss = float('inf')
    checkpoint_interval = 25

    for epoch in range(EPOCHS):
        total_loss = 0.0
        for i, (xx, yy_truth) in enumerate(train_loader):
            xx, yy_truth = xx.to(device), yy_truth.to(device)
            
            loss = 0
            current_input = F.interpolate(xx, size=MODEL_RESOLUTION, mode='bilinear', align_corners=False)
            
            for t in range(UNROLL_STEP):
                prediction = model(current_input)
                true_next_step = F.interpolate(yy_truth[:, t, :, :, :], size=MODEL_RESOLUTION, mode='bilinear', align_corners=False)
                loss += loss_fn(prediction, true_next_step)
                
                current_input = torch.cat([current_input[:, out_channels:, :, :], prediction], dim=1)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            total_loss += loss.item() / UNROLL_STEP

            if (i + 1) % 10 == 0:
                print(f"  Epoch {epoch+1}/{EPOCHS}, Batch {i+1}/{len(train_loader)}, Avg. Loss: {loss.item() / UNROLL_STEP:.6f}")

        avg_loss = total_loss / len(train_loader)
        log_line = f"{epoch+1},{avg_loss:.6f}\n"
        print(f"Epoch {epoch+1}/{EPOCHS} finished. Average Loss: {avg_loss:.6f}")
        
        with open(log_path, 'a') as f:
            f.write(log_line)

        is_best = avg_loss < best_loss
        if is_best:
            best_loss = avg_loss
            
        if (epoch + 1) % checkpoint_interval == 0 or is_best:
            checkpoint_type = "best" if is_best else f"epoch_{epoch+1}"
            checkpoint_path = os.path.join(FINETUNED_MODEL_DIR, f'model_{checkpoint_type}_{timestamp}.pt')
            print(f"Saving {'best ' if is_best else ''}checkpoint to '{checkpoint_path}'...")
            
            checkpoint_to_save = {
                'model_state_dict': model.state_dict(),
                'epoch': epoch + 1,
                'loss': avg_loss,
                'optimizer_state_dict': optimizer.state_dict(),
                'best_loss': best_loss,
                'timestamp': timestamp
            }
            torch.save(checkpoint_to_save, checkpoint_path)
            print(f"Checkpoint saved. {'New best model!' if is_best else ''}")

    final_model_path = os.path.join(FINETUNED_MODEL_DIR, f'finetuned_model_final_{timestamp}.pt')
    print(f"Saving final model to '{final_model_path}'...")
    final_checkpoint = {
        'model_state_dict': model.state_dict(),
        'epoch': EPOCHS,
        'loss': avg_loss,
        'optimizer_state_dict': optimizer.state_dict(),
        'best_loss': best_loss,
        'timestamp': timestamp
    }
    torch.save(final_checkpoint, final_model_path)
    print("Final model saved successfully.")

if __name__ == "__main__":
    main()
