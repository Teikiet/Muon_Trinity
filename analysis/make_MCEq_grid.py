#Generate MCEq Grid run once
from MCEq.core import MCEqRun
import crflux.models as crf
import MCEq.config as config
from MCEq.geometry.density_profiles import MSIS00Atmosphere
from MCEq.geometry.geometry import EarthGeometry
config.kernel_config= 'MKL'
config.e_min = 0.160 # GeV
config.integrator= 'euler'
import numpy as np
import crflux.models as pm
import matplotlib.pyplot as plt
# generate_mceq_grid.py
import pickle
# Initalize MCEq by creating the user interface object MCEqRun
# Initialize custom atmosphere model with coordinates
# Note: Using custom location with latitude/longitude parameters
mag = 0
# Initialize MCEq with custom atmosphere and Frisco Peak location
mceq = MCEqRun(
    # interaction interaction model
    interaction_model='SIBYLL23C',
    # Primary cosmic ray model
    primary_model=(pm.GlobalSplineFitBeta, None),
    # Set to 0° for horizontal muons
    theta_deg=0.,
    # Use custom atmosphere and geometry
    density_model=("CORSIKA", ('USStd', None)),
    #density_model=("MSIS00_IC", ('FriscoPeak', 'January')),
)
#mceq.set_density_model(("CORSIKA", ('USStd', None)))
earth_geom = EarthGeometry()
def get_spectrum(theta, Nh, mag=0, h_obs_m=2944, h_max_m=112800.0, E_min=1e3, E_max=1e7):
    mceq.set_theta_deg(theta)
    mceq.integration_path = None 
    angle_rad = np.radians(theta)
    #L = earth_geom.path_len(angle_rad)
    # l_grid goes L -> 0 (Bottom to Top)
    #L_grid = np.linspace(L, 0, Nh) 
    # h_grid goes h(0) -> h(L) which is TOP to BOTTOM (e.g., 112.8km -> 0km)
    #H_grid = np.array([earth_geom.h(l, angle_rad) for l in (L) - L_grid])
    # count H_grid to observation height:
    #H_grid = H_grid[H_grid>=h_obs_m*1e2] # 2944m in cm
    #H_grid = H_grid[H_grid<=h_max_m*1e2] # 112800m in cm
    #L_grid = L_grid[-len(H_grid):] # Corresponding L_grid for the valid H_grid
    # X_grid goes X(Top) -> X(Bottom) which is 0 -> 1030 (Correct for solver)
    #H_grid = np.linspace((h_max_m*1e2), (h_obs_m*1e2), Nh) # in cm
    #X_grid = mceq.density_model.h2X(H_grid)
    
    min_X = mceq.density_model.h2X(h_max_m*1e2)
    max_X = mceq.density_model.h2X(h_obs_m*1e2)
    X_grid = np.logspace(np.log10(min_X), np.log10(max_X), Nh) #np.linspace(min_X, max_X,Nh)##
    H_grid = mceq.density_model.X2h(X_grid)
    L_grid = np.array([earth_geom.delta_l(h, angle_rad) for h in H_grid])
    mceq.solve(int_grid=X_grid)
    longitudinal_spectrum = [
        mceq.get_solution('total_mu-', grid_idx=idx, mag=mag) + mceq.get_solution('total_mu+', grid_idx=idx, mag=mag)
        for idx in range(len(X_grid))]
    CDF = np.array(longitudinal_spectrum)
    CDF_zero = np.concatenate([np.zeros((1, CDF.shape[1])), CDF], axis=0) #(Nh+1, NE)
    d_CDF = np.diff(CDF_zero, axis=0) #(Nh, NE)
    E_mask = (mceq.e_grid >= E_min) & (mceq.e_grid <= E_max)
    sum_mceq = np.trapezoid(CDF[-1, :][E_mask], mceq.e_grid[E_mask])
    return d_CDF, L_grid, X_grid, H_grid, sum_mceq


# --- Grid settings ---
Nh    = 1000
Nzen  = 100
min_zen, max_zen = 80, 90
min_E,  max_E    = 1, 1e7

cos_edges   = np.linspace(np.cos(np.radians(max_zen)),
                          np.cos(np.radians(min_zen)), Nzen + 1)
zen_centers = np.degrees(np.arccos(0.5*(cos_edges[:-1] + cos_edges[1:])))

# --- Run MCEq at each zenith and collect results ---
grid_data = {
    'zen_centers': zen_centers,
    'cos_edges':   cos_edges,
    'E_GeV':       mceq.e_grid.copy(),
    'min_E':       min_E,
    'max_E':       max_E,
    'Nh':          Nh,
    'per_zenith':  [],   # list of dicts, one per zenith
    'meta': {
        'interaction_model': 'SIBYLL23C',
        'primary_model':     'GlobalSplineFitBeta',
        'density_model':     ('CORSIKA', ('USStd', None)),
    },
}

for theta in zen_centers:
    print(f"Running MCEq at zenith {theta:.3f} deg ...")
    d_CDF, L_grid, X_grid, H_grid, sum_mceq = get_spectrum(
        theta, Nh, E_min=min_E, E_max=max_E
    )
    grid_data['per_zenith'].append({
        'theta':   theta,
        'd_CDF':   d_CDF,        # (Nh-1, nE) muon production dN/dh/dE
        'L_grid':  L_grid,
        'X_grid':  X_grid,
        'H_grid':  H_grid,       # cm, top -> bottom
        'sum_mceq': sum_mceq,    # raw ground flux (no efficiency)
    })

# --- Pickle it ---
out_path = f'/uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/data/mceq_grid_{min_zen}_{max_zen}deg.pkl'
with open(out_path, 'wb') as f:
    pickle.dump(grid_data, f, protocol=pickle.HIGHEST_PROTOCOL)

print(f"Saved MCEq grid to {out_path}")
