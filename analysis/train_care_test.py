import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader, random_split
import argparse

CACHE_PATH = ("/uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/cache/"
              "care_cube_pid13_E1e+03-1e+05_R5_thr20_sz0.5_combined.npz")
device = torch.device("cpu")
torch.manual_seed(0); np.random.seed(0)

# ── Command-line arguments ────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="CARE regressor with temporal downsampling")
parser.add_argument("--time_step", type=int, default=2, 
                    help="Sample every N time steps (e.g., 2 or 4)")
parser.add_argument("--agg_method", type=str, default="max", 
                    choices=["max", "mean", "median"],
                    help="Aggregation method for each time-step window")
args = parser.parse_args()

TIME_STEP = args.time_step
AGG_METHOD = args.agg_method

print(f"Time-step sampling: every {TIME_STEP} | Aggregation: {AGG_METHOD}")

# ── Load cached data ──────────────────────────────────────────────────────────
d = np.load(CACHE_PATH, allow_pickle=True)
X = d["X"].astype(np.float32)              # (N, 200, 16, 16), already /150
y = d["y"].astype(np.float32)              # (N, 5): log10E, zen, az, tel_x, tel_z
target_names = list(d["target_names"])
N, T, H, W = X.shape
print(f"Original X shape: {X.shape}, y shape: {y.shape}, target_names: {target_names}")

# ── Downsample time dimension ─────────────────────────────────────────────────
def downsample_time(X_orig, time_step, agg_method="max"):
    """
    Downsample X along time axis by selecting representative frames.
    
    For each window of size `time_step`, compute total PE per spatial location,
    then select/aggregate according to agg_method.
    
    Args:
        X_orig: (N, T, 16, 16) array
        time_step: window size (e.g., 2 or 4)
        agg_method: "max" | "mean" | "median"
    
    Returns:
        X_ds: (N, T_new, 16, 16) downsampled array
        T_new: new time dimension
    """
    N, T, H, W = X_orig.shape
    if time_step == 1:
        return X_orig.copy(), T
    # Pad T to be divisible by time_step
    T_pad = ((T + time_step - 1) // time_step) * time_step
    X_padded = np.pad(X_orig, ((0, 0), (0, T_pad - T), (0, 0), (0, 0)), 
                      mode='constant', constant_values=0)
    
    # Reshape into windows: (N, T_new, time_step, H, W)
    T_new = T_pad // time_step
    X_windowed = X_padded.reshape(N, T_new, time_step, H, W)
    
    if agg_method == "max":
        # For each window, find the frame with max total PE, select that frame
        total_pe = X_windowed.sum(axis=(3, 4))  # (N, T_new, time_step)
        max_idx = np.argmax(total_pe, axis=2)   # (N, T_new) indices
        X_ds = X_windowed[np.arange(N)[:, None], 
                          np.arange(T_new)[None, :], 
                          max_idx]  # (N, T_new, H, W)
    
    elif agg_method == "mean":
        # Average all frames in each window
        X_ds = X_windowed.mean(axis=2)  # (N, T_new, H, W)
    
    elif agg_method == "median":
        # Median of all frames in each window
        X_ds = np.median(X_windowed, axis=2)  # (N, T_new, H, W)
    
    else:
        raise ValueError(f"Unknown agg_method: {agg_method}")
    
    return X_ds.astype(np.float32), T_new

X_ds, T_new = downsample_time(X, TIME_STEP, AGG_METHOD)
print(f"Downsampled X shape: {X_ds.shape} (time: {T} -> {T_new})")

# Treat the downsampled time frames as input CHANNELS to a 2D CNN: (N, T_new, 16, 16)
X_t = torch.from_numpy(X_ds)                  # (N, T_new, H, W)
y_t = torch.from_numpy(y)

# ── Split ─────────────────────────────────────────────────────────────────────
n_val = max(1, int(0.2 * N))
n_train = N - n_val
train_ds, val_ds = random_split(TensorDataset(X_t, y_t), [n_train, n_val],
                                generator=torch.Generator().manual_seed(0))

# ── Standardize targets using TRAIN split only ────────────────────────────────
train_idx = train_ds.indices
y_mean = y_t[train_idx].mean(0)
y_std  = y_t[train_idx].std(0).clamp_min(1e-6)
print("target mean", y_mean.tolist())
print("target std ", y_std.tolist())

def standardize(yb):  return (yb - y_mean) / y_std
def destd(yb):        return yb * y_std + y_mean

train_dl = DataLoader(train_ds, batch_size=16, shuffle=True)
val_dl   = DataLoader(val_ds, batch_size=32, shuffle=False)

# ── Model: small CNN over (T_new channels, 16x16) ──────────────────────────────
class CareRegressor(nn.Module):
    def __init__(self, in_ch, n_out=5, p=0.3):
        super().__init__()
        self.features = nn.Sequential(
            nn.Conv2d(in_ch, 64, 3, padding=1), nn.BatchNorm2d(64), nn.ReLU(),
            nn.Conv2d(64, 64, 3, padding=1),    nn.BatchNorm2d(64), nn.ReLU(),
            nn.MaxPool2d(2),                                    # 16->8
            nn.Conv2d(64, 128, 3, padding=1),   nn.BatchNorm2d(128), nn.ReLU(),
            nn.MaxPool2d(2),                                    # 8->4
            nn.AdaptiveAvgPool2d(1),                            # ->1x1
        )
        self.head = nn.Sequential(
            nn.Flatten(),
            nn.Linear(128, 128), nn.ReLU(), nn.Dropout(p),
            nn.Linear(128, n_out),
        )
    def forward(self, x):
        return self.head(self.features(x))

model = CareRegressor(in_ch=T_new).to(device)
y_mean, y_std = y_mean.to(device), y_std.to(device)

opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-3)
sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, "min", patience=15, factor=0.5)
loss_fn = nn.SmoothL1Loss()   # robust to the small-sample outliers

# ── Train ─────────────────────────────────────────────────────────────────────
EPOCHS = 100
best_val = np.inf; best_state = None
for ep in range(1, EPOCHS + 1):
    model.train(); tr_loss = 0.0
    for xb, yb in train_dl:
        xb, yb = xb.to(device), yb.to(device)
        opt.zero_grad()
        pred = model(xb)
        loss = loss_fn(pred, standardize(yb))
        loss.backward(); opt.step()
        tr_loss += loss.item() * len(xb)
    tr_loss /= n_train

    model.eval(); vl_loss = 0.0
    with torch.no_grad():
        for xb, yb in val_dl:
            xb, yb = xb.to(device), yb.to(device)
            vl_loss += loss_fn(model(xb), standardize(yb)).item() * len(xb)
    vl_loss /= n_val
    sched.step(vl_loss)

    if vl_loss < best_val:
        best_val = vl_loss
        best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
    if ep % 20 == 0 or ep == 1:
        print(f"epoch {ep:3d} | train {tr_loss:.4f} | val {vl_loss:.4f}")

model.load_state_dict(best_state)
print(f"Best val loss {best_val:.4f}")

# ── Evaluate in physical units (per-target MAE) ───────────────────────────────
model.eval()
with torch.no_grad():
    xb, yb = next(iter(DataLoader(val_ds, batch_size=n_val)))
    pred_phys = destd(model(xb.to(device))).cpu()
mae = (pred_phys - yb).abs().mean(0)
print("\nPer-target validation MAE (physical units):")
for name, m in zip(target_names, mae.tolist()):
    print(f"  {name:10s}: {m:.4f}")

# ── Save model + normalization ────────────────────────────────────────────────
model_name = f"model_T{TIME_STEP}_{AGG_METHOD}.pt"
torch.save({
    "state_dict": best_state,
    "y_mean": y_mean.cpu(), "y_std": y_std.cpu(),
    "target_names": target_names, 
    "T_original": T,
    "T_downsampled": T_new,
    "time_step": TIME_STEP,
    "agg_method": AGG_METHOD,
}, CACHE_PATH.replace(".npz", f"_{model_name}"))
print(f"Saved model: {model_name}")