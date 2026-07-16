import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader, random_split

CACHE_PATH = ("/uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/cache/"
              "care_cube_pid13_E1e+03-1e+04_R5_thr20.npz")
device = torch.device("cpu")
torch.manual_seed(0); np.random.seed(0)

# ── Load cached data ──────────────────────────────────────────────────────────
d = np.load(CACHE_PATH, allow_pickle=True)
X = d["X"].astype(np.float32)              # (N, 200, 16, 16), already /150
y = d["y"].astype(np.float32)              # (N, 5): log10E, zen, az, tel_x, tel_z
orig_names = list(d["target_names"])
N, T, H, W = X.shape
print("X", X.shape, "y", y.shape, orig_names)

# ── Build angle-encoded target vector ─────────────────────────────────────────
# New target layout: [log10E, sin_zen, cos_zen, sin_az, cos_az, tel_x, tel_z]
log10E = y[:, 0:1]
zen_rad = np.deg2rad(y[:, 1])
az_rad  = np.deg2rad(y[:, 2])
tel     = y[:, 3:5]

y_enc = np.concatenate([
    log10E,
    np.sin(zen_rad)[:, None], np.cos(zen_rad)[:, None],
    np.sin(az_rad)[:, None],  np.cos(az_rad)[:, None],
    tel,
], axis=1).astype(np.float32)             # (N, 7)

target_names = ["log10E", "sin_zen", "cos_zen", "sin_az", "cos_az",
                "tel_x_m", "tel_z_m"]
N_OUT = y_enc.shape[1]

X_t = torch.from_numpy(X)
y_t = torch.from_numpy(y_enc)
# keep original (deg) targets around for physical-unit MAE evaluation
y_deg = torch.from_numpy(y)

# ── Split ─────────────────────────────────────────────────────────────────────
n_val = max(1, int(0.2 * N))
n_train = N - n_val
full_ds = TensorDataset(X_t, y_t, y_deg)
train_ds, val_ds = random_split(full_ds, [n_train, n_val],
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

# ── Model: small CNN over (T=200 channels, 16x16) ─────────────────────────────
class CareRegressor(nn.Module):
    def __init__(self, in_ch=T, n_out=N_OUT, p=0.3):
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

model = CareRegressor().to(device)
y_mean, y_std = y_mean.to(device), y_std.to(device)

opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-3)
sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, "min", patience=15, factor=0.5)
loss_fn = nn.SmoothL1Loss()   # robust to the small-sample outliers

# ── Helper: decode standardized net output -> physical [log10E, zen_deg, az_deg, tel_x, tel_z]
def decode_to_phys(pred_std):
    p = destd(pred_std)                          # (B, 7) in encoded space
    log10E = p[:, 0]
    zen = torch.atan2(p[:, 1], p[:, 2])          # sin, cos
    az  = torch.atan2(p[:, 3], p[:, 4])
    zen_deg = torch.rad2deg(zen) % 360.0
    az_deg  = torch.rad2deg(az)  % 360.0
    tel_x = p[:, 5]; tel_z = p[:, 6]
    return torch.stack([log10E, zen_deg, az_deg, tel_x, tel_z], dim=1)

# ── Train ─────────────────────────────────────────────────────────────────────
EPOCHS = 200
best_val = np.inf; best_state = None
for ep in range(1, EPOCHS + 1):
    model.train(); tr_loss = 0.0
    for xb, yb, _ in train_dl:
        xb, yb = xb.to(device), yb.to(device)
        opt.zero_grad()
        pred = model(xb)
        loss = loss_fn(pred, standardize(yb))
        loss.backward(); opt.step()
        tr_loss += loss.item() * len(xb)
    tr_loss /= n_train

    model.eval(); vl_loss = 0.0
    with torch.no_grad():
        for xb, yb, _ in val_dl:
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
eval_names = ["log10E", "zen_deg", "az_deg", "tel_x_m", "tel_z_m"]
model.eval()
with torch.no_grad():
    xb, yb_enc, yb_deg = next(iter(DataLoader(val_ds, batch_size=n_val)))
    pred_phys = decode_to_phys(model(xb.to(device))).cpu()

# angular MAE with wraparound for azimuth/zenith
diff = pred_phys - yb_deg
for j in (1, 2):  # zen_deg, az_deg columns -> wrap to [-180,180]
    diff[:, j] = (diff[:, j] + 180.0) % 360.0 - 180.0
mae = diff.abs().mean(0)

print("\nPer-target validation MAE (physical units):")
for name, m in zip(eval_names, mae.tolist()):
    print(f"  {name:10s}: {m:.4f}")

# ── Save model + normalization ────────────────────────────────────────────────
torch.save({
    "state_dict": best_state,
    "y_mean": y_mean.cpu(), "y_std": y_std.cpu(),
    "target_names": target_names, "T": T,
    "encoding": "angles_sincos",
}, CACHE_PATH.replace(".npz", "_model_sincos.pt"))
print("Saved model.")