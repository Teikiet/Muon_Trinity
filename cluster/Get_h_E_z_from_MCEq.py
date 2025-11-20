from MCEq.core import MCEqRun
import crflux.models as crf
import MCEq.config as config
from MCEq.geometry.geometry import EarthGeometry
config.kernel_config= 'MKL'
config.integrator= 'euler'
import numpy as np
import crflux.models as pm
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

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
def get_spectrum(theta):
    mceq.set_theta_deg(theta)
    angle_rad = np.radians(theta) # Convert zenith angle to radians
    L = earth_geom.path_len(angle_rad) # Calculate the total Path Length for the given zenith angle
    l_grid = np.linspace(L, 0,1000) # Define the grid for the Path Length in cm
    h_grid = np.array([earth_geom.h(l, angle_rad) for l in (L)-l_grid])
    X_grid = mceq.density_model.h2X(h_grid)
    mceq.solve(int_grid=X_grid)
    longitudinal_spectrum = []
    for idx in range(len(X_grid)):
        #print('Reading solution at X = {0:5.2f} g/cm2'.format(X_grid[idx]))
        longitudinal_spectrum.append(mceq.get_solution('total_mu-', grid_idx=idx, mag=mag))# + mceq.get_solution('total_mu+', grid_idx=idx, mag=mag))
    return np.array(longitudinal_spectrum), l_grid, X_grid, h_grid

# 1. Prepare your full 3D flux distribution
flux_3d = []
zenith_angles = np.linspace(60, 88, 10)  # Your zenith range
Area = np.pi*(150)**2 # Define the area of the detector in cm^2
d_theta = np.radians(28)
d_phi = np.radians(28)
Solid_Angle = np.pi*np.tan(d_theta)**2
E = mceq.e_grid  # Energy grid from MCEq
for zenith in zenith_angles:
    m, l, X, h = get_spectrum(zenith)
    muon_flux = m * Solid_Angle * Area
    # Modify energy spectrum from E^-3 to E^-1
    muon_flux_modified = muon_flux * E[np.newaxis, :]**2  # Broadcast over height dimension
    flux_3d.append(muon_flux_modified)

flux_3d = np.array(flux_3d)  # Shape: (n_zenith, n_heights, n_energies)
max_flux = np.max(flux_3d)

# 2. Rejection sampling function
def sample_correlated_muon():
    while True:
        # Sample uniformly from parameter space
        theta_idx = np.random.randint(0, len(zenith_angles))
        h_idx = np.random.randint(0, len(h))
        e_idx = np.random.randint(0, len(E))
        
        # Get flux at this point
        flux_value = flux_3d[theta_idx, h_idx, e_idx]
        
        # Accept/reject
        if np.random.random() < flux_value / max_flux:
            return E[e_idx], h[h_idx], zenith_angles[theta_idx]

# 3. Generate samples
n_samples = 1
samples = [sample_correlated_muon() for _ in range(n_samples)]
energies, heights, zeniths = zip(*samples)
print("Sampled Energies:", energies)
print("Sampled Heights:", heights)
print("Sampled Zeniths:", zeniths)