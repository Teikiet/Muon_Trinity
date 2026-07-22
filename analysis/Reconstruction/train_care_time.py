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
target_names = list(d["target_names"])
N, T, H, W = X.shape
print("X", X.shape, "y", y.shape, target_names)

# Keep the time axis EXPLICIT: (N, T, 1, 16, 16) — T is a sequence, not channels
X_t = torch.from_numpy(X).unsqueeze(2)     # (N, T, 1, H, W)
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

# ── Model: per-frame CNN encoder + GRU over time ──────────────────────────────
class FrameEncoder(nn.Module):
    """Small 2D CNN applied to a single 16x16 frame -> feature vector."""
    def __init__(self, feat_dim=64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv2d(1, 32, 3, padding=1), nn.BatchNorm2d(32), nn.ReLU(),
            nn.MaxPool2d(2),                                    # 16->8
            nn.Conv2d(32, 64, 3, padding=1), nn.BatchNorm2d(64), nn.ReLU(),
            nn.MaxPool2d(2),                                    # 8->4
            nn.Conv2d(64, feat_dim, 3, padding=1), nn.BatchNorm2d(feat_dim), nn.ReLU(),
            nn.AdaptiveAvgPool2d(1),                            # ->1x1
            nn.Flatten(),                                       # (B, feat_dim)
        )
    def forward(self, x):
        return self.net(x)

class CareSeqRegressor(nn.Module):
    """CNN encodes each frame; GRU reads the frame sequence in time order."""
    def __init__(self, feat_dim=64, hidden=128, n_out=5, p=0.3):
        super().__init__()
        self.encoder = FrameEncoder(feat_dim)
        self.gru = nn.GRU(input_size=feat_dim, hidden_size=hidden,
                          num_layers=1, batch_first=True)
        self.head = nn.Sequential(
            nn.Linear(hidden, 128), nn.ReLU(), nn.Dropout(p),
            nn.Linear(128, n_out),
        )
    def forward(self, x):                  # x: (B, T, 1, H, W)
        B, T, C, H, W = x.shape
        x = x.view(B * T, C, H, W)         # fold time into batch for the CNN
        f = self.encoder(x)                # (B*T, feat_dim)
        f = f.view(B, T, -1)               # unfold -> ordered sequence
        out, h = self.gru(f)               # h: (1, B, hidden) last hidden state
        return self.head(h[-1])            # use final time step

model = CareSeqRegressor().to(device)
y_mean, y_std = y_mean.to(device), y_std.to(device)

opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-3)
sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, "min", patience=15, factor=0.5)
loss_fn = nn.SmoothL1Loss()

# ── Train ─────────────────────────────────────────────────────────────────────
EPOCHS = 300
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
torch.save({
    "state_dict": best_state,
    "y_mean": y_mean.cpu(), "y_std": y_std.cpu(),
    "target_names": target_names, "T": T,
}, CACHE_PATH.replace(".npz", "_seqmodel.pt"))
print("Saved model.")