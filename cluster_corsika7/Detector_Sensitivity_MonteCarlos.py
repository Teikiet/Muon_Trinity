import numpy as np
import random
import matplotlib.pyplot as plt
import math
import sys
import os
from typing import Dict, List, Tuple
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple
from scipy.stats import linregress
from matplotlib.colors import LogNorm
from eventio import IACTFile
from corsikaio import CorsikaParticleFile
import numpy as np
from scipy.spatial.transform import Rotation as R
import argparse
import pickle

"""The x_base is pointing at the sky and y_base is pointing to the North-East direction."""

def find_equivalent_z_rotation(theta_y_deg, theta_z_deg):
    """
    Find equivalent Z rotation after Y and Z rotations.
    Initial vector: [0, 0, 1] (pointing in z direction)
    """
    # Initial vector pointing in z direction
    initial_vector = np.array([0, 0, 1])
    
    # Apply X rotation first, then Y rotation
    R_y = R.from_euler('y', theta_y_deg, degrees=True)
    R_z = R.from_euler('z', theta_z_deg, degrees=True)

    # Combined rotation: R_total = R_z * R_y
    R_total = R_z * R_y

    # Apply to initial vector
    final_vector = R_total.apply(initial_vector)
    #angle between final vector and z-axis
    theta_z_equivalent = np.degrees(np.arccos(final_vector[2]))
    
    return theta_z_equivalent, final_vector

def get_lambda_fit():
    # Data from user
    x = np.array([
        0.7644305772230889,
    1.5288611544461779,
    2.293291731669267,
    3.2761310452418098,
    4.149765990639626,
    5.6786271450858035,
    6.1154446177847115,
    7.316692667706708,
    6.9890795631825275,
    8.190327613104525,
    9.391575663026522,
    10.265210608424338,
    11.24804992199688,
    12.230889235569423,
    13.213728549141965,
    14.087363494539781

    ])

    y = np.array([
        6.674311926605505,
        5.8348623853211015,
        5.277522935779817,
        4.727064220183487,
        4.217889908256881,
        3.729357798165138,
        3.295871559633028,
        2.883027523,
        2.793577981651376,
        2.3738532110091746,
        1.926605504587156,
        1.4862385321100917,
        1.1972477064220184,
        0.805045872,
        0,
        -0.275229358

    ])

    # Keep only positive y
    x_fit = x
    y_fit = np.log(10**y)

    # Linear regression
    slope, intercept, r_value, p_value, std_err = linregress(x_fit, y_fit)

    # Compute lambda
    lambda_fit = -1 / slope
    return lambda_fit

def noise_function(H, lambda_fit = 0.8946330713407988):
    size = H[0].shape[1]
    noise_matrix = []
    for i in range(len(H)):
        noise = np.random.exponential(scale=1/lambda_fit, size=(size,size))
        noise_matrix.append(noise)
    return np.array(noise_matrix)

def camera_response_function(t, t_photon, rise=5, fall=15):
    signal = np.zeros_like(t)
    for i, ti in enumerate(t):
        if ti < t_photon:
            signal[i] = np.exp(-(t_photon - ti)**2 / (2 * rise**2))
        else:
            signal[i] = np.exp(-(ti - t_photon)**2 / (2 * fall**2))
    if np.max(signal) == 0:
        print("Warning: PE signal is zero, check the time range and t0 values.")
        print("t0:", t_photon)
        return signal
    else:
        signal /= np.max(signal)
        #signal *= 24.1 #1PE = 24.1 DC
        #signal *= 0.1 #1DC = 0.1 mV
        return signal #in mV

class VAtmosAbsorption:
    """Atmospheric absorption of Cherenkov photons"""
    
    def __init__(self, model: str, seed: int = 0, source_file: str = ""):
        """
        Initialize atmospheric absorption model
        
        Args:
            model: atmospheric absorption model ("CORSIKA", "kascade", "us76_new", etc.)
            seed: random seed
            source_file: path to source file (optional)
        """
        self.model = model.lower()
        self.observation_level = (294400+1500)/100  # default observation level in [m]
        self.random_state = np.random.RandomState(seed)
        self.min_wave = 300.0  # default CORSIKA values
        self.max_wave = 600.0  # default CORSIKA values
        
        print(f"Atmospheric extinction model: {model}")
        
        # Initialize data structures
        self.coeff: Dict[int, List[float]] = {}
        self.coeff_obs: Dict[int, float] = {}
        self.extint: List[List[float]] = []
        
        # Load appropriate model data
        if self.model == "modtran5":
            self.source_file = source_file
            self._read_extint_m5()
            
        elif self.model in ["corsika"]:
            self.source_file = source_file if source_file else "data/atmabs.dat"
            self._read_corsika_atmabs()
            
        elif self.model == "kascade":
            self.source_file = source_file if source_file else "data/kextint.dat"
            self._read_extint(180)
            
        elif self.model in ["us76_new", "us76.50km", "us76.50", "modtran4"]:
            self.model = "modtran4"
            self.source_file = source_file if source_file else "data/us76.50km.ext"
            self._read_extint(200)
            
        elif self.model in ["us76.23km", "us76.23", "modtran4_2"]:
            self.model = "modtran4_2"
            self.source_file = source_file if source_file else "data/us76.23km.ext"
            self._read_extint_f2(200)
            
        elif self.model == "artemis":
            self.source_file = source_file if source_file else "data/extinction_uv.dat"
            self._read_extint(180)
            
        elif self.model == "noextinction":
            self.source_file = "noExtinction"
            print("VAtmosAbsorption: no atmospheric extinction applied")
            
        else:
            print(f"VAtmosAbsorption::__init__: error, unknown model: {model}")
            sys.exit(-1)

    def set_observation_level(self, obs_level: float) -> None:
        """
        Set observation level
        
        Args:
            obs_level: observation level in [m]
        """
        self.observation_level = obs_level
        
        if self.model == "corsika":
            xobs = int(obs_level / 1000)
            for wavelength, coeffs in self.coeff.items():
                if xobs < len(coeffs) - 1:
                    self.coeff_obs[wavelength] = self._get_linear_interpolate(
                        obs_level / 1000.0, xobs, xobs + 1, 
                        coeffs[xobs], coeffs[xobs + 1]
                    )

    def prob_atm_absorbed(self, wavelength: float, z_emis: float, w_emis: float) -> Tuple[float, float]:
        """
        Calculate survival probability for photons
        
        Args:
            wavelength: wavelength in nm
            z_emis: emission height in m
            w_emis: cos of emission angle (cos theta)
            
        Returns:
            Tuple of (survival probability, optical depth)
        """
        # Fixed altitude steps [m]
        altitude_step = 1000.0
        opt_depth = 0.0
        
        # CORSIKA model
        if self.model == "corsika":
            wl_min = 180.0
            wl_step = 5.0
            
            # Calculate reference wavelength and index for interpolations
            riwl = 1 + int((wavelength - wl_min) / wl_step)
            wli0 = riwl * int(wl_step) + int(wl_min - wl_step)
            wli1 = riwl * int(wl_step) + int(wl_min)
            if wli0 in self.coeff and wli1 in self.coeff:
                # Consider atmospheric extinction
                htkm = z_emis / 1000  # m -> km
                hti0 = int(htkm)
                hti1 = int(htkm) + 1
                
                if hti0 < 0:
                    phi0 = self.coeff[wli0][0]
                    phi1 = self.coeff[wli1][0]
                elif hti1 > 50:
                    if len(self.coeff[wli0]) > 50 and len(self.coeff[wli1]) > 50:
                        phi0 = self.coeff[wli0][50]
                        phi1 = self.coeff[wli1][50]
                    else:
                        phi0 = phi1 = 0.0
                else:
                    # Interpolation in height
                    fx0 = self.coeff[wli0][hti0]
                    fx1 = self.coeff[wli0][hti1]
                    phi0 = self._get_linear_interpolate(htkm, hti0, hti1, fx0, fx1)
                    phi0 -= self.coeff_obs.get(wli0, 0.0)
                    
                    fx0 = self.coeff[wli1][hti0]
                    fx1 = self.coeff[wli1][hti1]
                    phi1 = self._get_linear_interpolate(htkm, hti0, hti1, fx0, fx1)
                    phi1 -= self.coeff_obs.get(wli1, 0.0)
                
                coatex = self._get_linear_interpolate(wavelength, wli0, wli1, phi0, phi1)
                opt_depth = coatex / w_emis if w_emis != 0 else 0.0
                probs = math.exp(-opt_depth)
                return probs, opt_depth
            else:
                print("VAtmosAbsorption::prob_atm_absorbed: coeff. matrix not valid")
                sys.exit(-1)
                
        elif self.model == "noextinction":
            return 1.0, 0.0
            
        # Other models (kascade, modtran4, etc.)
        else:
            hobs = int(self.observation_level)
            z = int(z_emis)
            wave = int(wavelength)
            
            min_wave = 180
            if self.model in ["modtran4", "modtran4_2"]:
                min_wave = 200
            elif self.model == "modtran5":
                min_wave = 205
                
            if wave < min_wave:
                return 0.0, 0.0
                
            # Step size fixed to 5 nm
            iwave1 = (wave - min_wave) // 5
            iwave2 = (wave - min_wave) // 5 + 1
            
            ihobs = hobs // int(altitude_step)
            
            if ihobs + 1 < len(self.extint) and iwave2 < len(self.extint[0]):
                p1 = (self.extint[ihobs][iwave1] + 
                      (self.extint[ihobs + 1][iwave1] - self.extint[ihobs][iwave1]) * 
                      (hobs / altitude_step - ihobs))
                p2 = (self.extint[ihobs][iwave2] + 
                      (self.extint[ihobs + 1][iwave2] - self.extint[ihobs][iwave2]) * 
                      (hobs / altitude_step - ihobs))
                tlow = self._get_linear_interpolate(
                    wave, iwave1 * 5 + min_wave, iwave2 * 5 + min_wave, p1, p2
                )
                
                ihgt = z // int(altitude_step)
                iext_index = ihgt
                
                # Extinction calculated up to 50 km only
                if iext_index > len(self.extint) - 2:
                    iext_index = len(self.extint) - 2
                    
                p1 = (self.extint[iext_index][iwave1] + 
                      (self.extint[iext_index + 1][iwave1] - self.extint[iext_index][iwave1]) * 
                      (z / altitude_step - ihgt))
                p2 = (self.extint[iext_index][iwave2] + 
                      (self.extint[iext_index + 1][iwave2] - self.extint[iext_index][iwave2]) * 
                      (z / altitude_step - ihgt))
                thigh = self._get_linear_interpolate(
                    wave, iwave1 * 5 + min_wave, iwave2 * 5 + min_wave, p1, p2
                )
                
                if w_emis != 0.0:
                    opt_depth = -(tlow - thigh) / w_emis
                    atm_prob = math.exp(-opt_depth)
                else:
                    atm_prob = 0.0
                    
                # Check validity of results
                if not math.isfinite(atm_prob):
                    print(f"VAtmosAbsorption::prob_atm_absorbed not normal {wavelength}\t{atm_prob}")
                    return 0.0, 0.0
                    
                return atm_prob, opt_depth
            else:
                return 0.0, 0.0
                
        return 0.0, 0.0

    def set_wavelength_interval(self, min_wave: float, max_wave: float) -> None:
        """Set wavelength interval"""
        self.min_wave = min_wave
        self.max_wave = max_wave

    def get_wavelength(self, emission_height: float, emission_angle: float) -> float:
        """
        Get random wavelength from Jelly formula and apply atmospheric absorption
        
        Args:
            emission_height: emission height of photon [m]
            emission_angle: emission angle of photon
            
        Returns:
            photon wavelength if photon survives atmosphere, otherwise -1
        """
        lambda_val = 1.0 / (1.0 / self.min_wave - 
                            self.random_state.uniform() * 
                            (1.0 / self.min_wave - 1.0 / self.max_wave))
        prob, _ = self.prob_atm_absorbed(lambda_val, emission_height, emission_angle)
        
        if self.random_state.uniform() > prob:
            return -1.0
            
        return lambda_val

    def _get_linear_interpolate(self, x: float, x0: float, x1: float, y0: float, y1: float) -> float:
        """Linear interpolation"""
        if x0 == x1:
            return 0.0
        return y0 + (y1 - y0) * (x - x0) / (x1 - x0)

    def _read_corsika_atmabs(self) -> None:
        """Read CORSIKA atmospheric absorption file"""
        if not os.path.exists(self.source_file):
            print(f"Atmospheric extinction file not found: {self.source_file}")
            sys.exit(-1)
            
        print(f"VAtmosAbsorption: reading atmospheric extinction file (from CORSIKA): {self.source_file}")
        
        with open(self.source_file, 'r') as f:
            lines = f.readlines()
            
        wl = 0
        i_coeff = []
        n_step = 0
        
        for line in lines:
            line = line.strip()
            if len(line) > 0:
                if len(line) < 5:  # wavelength line
                    if n_step > 0:
                        self.coeff[wl] = i_coeff.copy()
                    wl = int(line)
                    i_coeff.clear()
                    n_step += 1
                elif "EXTINCTION" not in line:
                    values = line.split()
                    i_coeff.extend([float(val) for val in values])
                    
        if n_step > 0:
            self.coeff[wl] = i_coeff

    def _read_extint(self, lambda_min: int) -> None:
        """Read extinction file (kascade format)"""
        lambda_max = 900
        step_size_needed = 5
        lsteps = 1 + (lambda_max - lambda_min) // step_size_needed
        alt_steps = 51
        
        # Initialize extint array
        self.extint = [[0.0 for _ in range(lsteps)] for _ in range(alt_steps)]
        
        if not os.path.exists(self.source_file):
            print(f"Unable to open input file: {self.source_file}")
            sys.exit(-1)
            
        with open(self.source_file, 'r') as f:
            lambda_val = lambda_min
            step = step_size_needed
            
            while lambda_val <= lambda_max:
                # Read wavelength
                if self.model != "artemis":
                    line = f.readline()
                    lambda_in = int(line.split()[0])
                else:
                    line = f.readline()
                    lambda_in = int(line.strip())
                    
                if lambda_val != lambda_in and self.model != "artemis":
                    print(f"VAtmosAbsorption::_read_extint: Error in reading {self.source_file}")
                    print(f"Lambda expected: {lambda_val}, Lambda read: {lambda_in}")
                    
                bin_idx = (lambda_val - lambda_min) // step_size_needed
                
                # Read extinction coefficients
                for i in range(alt_steps):
                    if self.model != "artemis":
                        val = float(f.read(12))  # Assuming fixed-width format
                    else:
                        val = float(f.readline().strip())
                        
                    self.extint[i][bin_idx] = val
                    if not math.isfinite(val) and abs(val) > 1e-5:
                        print(f"Warning: value not normal in {self.source_file} for lambda = {lambda_val} and altitude {i}")
                        self.extint[i][bin_idx] = 200.0
                        
                # Handle step size changes for specific wavelengths
                if self.model != "artemis":
                    if lambda_val == 270:
                        step = 10
                    elif lambda_val == 280:
                        step = 20
                    elif lambda_val == 400:
                        step = 50
                    elif lambda_val == 700:
                        step = 100
                        
                lambda_val += step

    def _read_extint_f2(self, lambda_min: int) -> None:
        """Read extinction file (F2 format)"""
        lambda_max = 900
        step_size_needed = 5
        lsteps = 1 + (lambda_max - lambda_min) // step_size_needed
        alt_steps = 51
        
        # Initialize extint array
        self.extint = [[0.0 for _ in range(lsteps)] for _ in range(alt_steps)]
        
        if not os.path.exists(self.source_file):
            print(f"VAtmosAbsorption::_read_extint_f2, source file not found: {self.source_file}")
            sys.exit(-1)
            
        print(f"\t reading atmospheric extinction file (from kascade): {self.source_file}")
        
        with open(self.source_file, 'r') as f:
            for line in f:
                line = line.strip()
                if len(line) > 0:
                    values = line.split()
                    lambda_val = int(values[0])
                    bin_idx = (lambda_val - lambda_min) // step_size_needed
                    
                    for i in range(alt_steps):
                        if i + 1 < len(values):
                            ext_val = float(values[i + 1])
                            if not math.isfinite(ext_val):
                                ext_val = 1e10
                            self.extint[i][bin_idx] = ext_val

    def _read_extint_m5(self) -> None:
        """Read MODTRAN5 extinction file"""
        lambda_max = 900
        lambda_min = 205
        step_size_needed = 5
        lsteps = 1 + (lambda_max - lambda_min) // step_size_needed
        
        if not os.path.exists(self.source_file):
            print(f"VAtmosAbsorption::_read_extint_m5, source file not found: {self.source_file}")
            sys.exit(-1)
        else:
            print(f"VAtmosAbsorption::_read_extint_m5, reading source file: {self.source_file}")
            
        # Initialize with zeros for 0-1 km and 0-0 km
        self.extint = [[0.0 for _ in range(lsteps)] for _ in range(2)]
        
        with open(self.source_file, 'r') as f:
            for line in f:
                line = line.strip()
                if len(line) > 0:
                    values = line.split()
                    lambda_val = int(values[0])
                    bin_idx = (lambda_val - lambda_min) // step_size_needed
                    
                    i = 2
                    for val_str in values[1:]:
                        ext = float(val_str)
                        minus_ln = -math.log(ext)
                        
                        # Extend extint if necessary
                        while len(self.extint) <= i:
                            self.extint.append([0.0 for _ in range(lsteps)])
                            
                        if bin_idx < len(self.extint[i]):
                            self.extint[i][bin_idx] = minus_ln + self.extint[i-1][bin_idx]
                        else:
                            self.extint[i].append(minus_ln + self.extint[i-1][bin_idx])
                            
                        if not math.isfinite(self.extint[i][bin_idx]):
                            self.extint[i][bin_idx] = 1e10
                            
                        i += 1

def pde_vs_wavelength(wl):
    wavelength = np.array([300, 320, 340, 360, 380, 400, 420, 440, 480, 500, 520, 560, 600, 640, 660, 700])
    pde = np.array([0.28, 0.38, 0.38, 0.38, 0.44, 0.46, 0.48, 0.48, 0.5, 0.5, 0.48, 0.4, 0.34, 0.28, 0.26, 0.22])
    wl1 = int(np.round(wl/20,0)*20)
    index = np.where(wavelength == wl1)[0]
    return pde[index][0] if index.size > 0 else 0.99  # Default value if wavelength not found

def survive_photon(wavelength, emission_height, emission_angle_rad, atmos):
    """Calculate survival probability for a photon"""
    cos_emission_angle = np.cos(emission_angle_rad)  # cos of emission angle
    try:
        survival_prob, _ = atmos.prob_atm_absorbed(wavelength, emission_height, cos_emission_angle)
    except:
        print(f"Warning: Wavelength {wavelength} nm not found in atmospheric absorption data.")
        print(f"height", emission_height, "m")
        survival_prob = 0.99  # Default value if wavelength not found
    survival_prob *= pde_vs_wavelength(wavelength)  # Apply the quantum efficiency
    survival_prob *= 0.92 #Transmission through the 6 mm thick entrance window of the camera
    survival_prob *= 0.85 #Reflection losses at the mirror surface
    check_survival = np.random.uniform(0, 1)
    if check_survival < survival_prob:
        return True
    else:
        return False

def filter_photons(x_shower, y_shower, T, wavelength, zem, theta, atmos, filter=True):
    print(f"Number of photons: {len(wavelength)}")
    # Check if photon survive:
    wavelength_new = np.array(wavelength)
    zem_new = np.array(zem)
    x_shower_new = np.array(x_shower)
    y_shower_new = np.array(y_shower)
    time_new = np.array(T)
    if filter == True:
        i=0
        for j in range(len(wavelength_new)):
            if not survive_photon(wavelength_new[j], zem_new[j] / 100, theta, atmos):
                wavelength_new[j] = 0  # Set to zero if photon does not survive
                i+=1
        # Filter out photons does not survive
        mask = wavelength_new > 0
        x_shower_new = x_shower_new[mask]
        y_shower_new = y_shower_new[mask]
        wavelength_new = wavelength_new[mask]
        zem_new = zem_new[mask] 
        time_new = time_new[mask]
        print(f"Number of photons after filtering: {len(time_new)}")
        return x_shower_new, y_shower_new, time_new, wavelength_new, zem_new, mask
    else:
        # If no filtering, return original arrays
        mask = np.ones_like(wavelength_new, dtype=bool)  # All photons survive
        return x_shower, y_shower, T, wavelength, zem, mask

def convert_coordinates_to_shower(x, y, cos_x, cos_y, azimuth_rad, zenith_deg=91):
    #Rotation around z axis such that the telescope is pointing at 280 azimuth
    x_ground = (np.array(x) * np.cos(azimuth_rad) - np.array(y) * np.sin(azimuth_rad))
    y_ground = np.array(x) * np.sin(azimuth_rad) + np.array(y) * np.cos(azimuth_rad)
    x_i = (cos_x*np.cos(azimuth_rad) - cos_y*np.sin(azimuth_rad))  # normalized x direction
    y_i = (cos_x*np.sin(azimuth_rad) + cos_y*np.cos(azimuth_rad))  # normalized y direction
    z_i = np.sqrt(1 - x_i**2 - y_i**2)  # normalized z direction
    #Project the ground position to the telescope position
    t = -x_ground*np.sin(np.radians(zenith_deg))/(x_i*np.sin(np.radians(zenith_deg)) + z_i*np.cos(np.radians(zenith_deg)))
    x_tele = x_ground + t * x_i
    y_tele = y_ground + t * y_i
    z_tele = t * z_i  # z coordinate of the telescope
    #Rotation around y axis:
    x_shower = x_tele * np.cos(np.radians(zenith_deg)) - z_tele * np.sin(np.radians(zenith_deg))
    y_shower = y_tele
    z_shower = x_tele * np.sin(np.radians(zenith_deg)) + z_tele * np.cos(np.radians(zenith_deg))
    cos_x = (x_i*np.cos(np.radians(zenith_deg)) - z_i*np.sin(np.radians(zenith_deg)))
    #cos_z = (x_i*np.sin(np.radians(zenith_deg)) + z_i*np.cos(np.radians(zenith_deg)))
    cos_y = y_i  # normalized y direction
    return x_shower, y_shower, cos_x, cos_y

def get_brightest_pixel_loc(all_x_shower, all_y_shower, x_range, y_range):
    # Find all locations with maximum value
    heatmap, xedges, yedges = np.histogram2d(
                    all_x_shower, 
                    all_y_shower, 
                    bins=50, 
                    range=[x_range, y_range]
                )
    max_value = np.max(heatmap)
    max_locations = np.where(heatmap == max_value)

    # If multiple pixels have the same max value, take the first one
    brightest_row = max_locations[0][0]
    brightest_col = max_locations[1][0]

    # Convert to coordinates
    brightest_y = (xedges[brightest_col] + xedges[brightest_col + 1]) / 2
    brightest_x = (yedges[brightest_row] + yedges[brightest_row + 1]) / 2
    return brightest_x, brightest_y, max_value

def get_photon_bunches(input_file, atmos, R=150, tele_id=0, filter=False, optics=False, get_brightest_pixel=False, shower_coordinates=True, zenith_tele=91, azimuth_tele=280):
    Event = []
    with IACTFile(input_file) as f:
        events = iter(f)
        event = next(events)
        azimuth_rad = (event.header['azimuth']) #event azimuth
        zenith_rad = (event.header['zenith']) #event zenith
        if tele_id ==0: #look for the first telescope
            X = event.photon_bunches[0]['x'] #in cm
            Y = event.photon_bunches[0]['y'] #in cm
            T = event.photon_bunches[0]['time'] #in ns
            cos_X = event.photon_bunches[0]['cx'] #cosine of the x direction
            cos_Y = event.photon_bunches[0]['cy'] #cosine of the y direction
            wavelength = event.photon_bunches[0]['wavelength'] #in nm
            zem = event.photon_bunches[0]['zem']  # in cm
        if tele_id ==1: # If the first telescope is not present, use the second one
            X = event.photon_bunches[1]['x'] #in cm
            Y = event.photon_bunches[1]['y'] #in cm
            T = event.photon_bunches[1]['time'] #in ns
            cos_X = event.photon_bunches[1]['cx'] #cosine of the x direction
            cos_Y = event.photon_bunches[1]['cy'] #cosine of the y direction
            wavelength = event.photon_bunches[1]['wavelength'] #in nm
            zem = event.photon_bunches[1]['zem']  # in cm
        #Convert to shower coordinates
        if shower_coordinates==True:
            x_shower, y_shower, cos_X, cos_Y = convert_coordinates_to_shower(X, Y, cos_X, cos_Y, azimuth_rad=np.radians(-azimuth_tele), zenith_deg=zenith_tele)
        else:
            x_shower, y_shower, cos_X, cos_Y = X, Y, cos_X, cos_Y
        x_brightest, y_brightest = 0, 0
        if get_brightest_pixel==True:
            x_brightness, y_brightness, brightness_value = get_brightest_pixel_loc(x_shower, y_shower, (-R, R), (-R, R))
            x_brightest += x_brightness  # Center the shower coordinates
            y_brightest += y_brightness  # Center the shower coordinates
            print(f"Brightness pixel location in cm: ({x_brightness:.2f}, {y_brightness:.2f}) with value {brightness_value:.2f}")
            x_brightness, y_brightness, brightness_value = get_brightest_pixel_loc(x_shower, y_shower, (-R/10, R/10), (-R/10, R/10))
            x_brightest += x_brightness  # Center the shower coordinates
            y_brightest += y_brightness  # Center the shower coordinates
            print(f"Brightness pixel location in cm: ({x_brightness:.2f}, {y_brightness:.2f}) with value {brightness_value:.2f}")
        if optics == True:
            # Simulate Davies-Cotton optics
            tele_radi = 150 # Telescope radius in cm
            cam_radi = 6 # Camera radius in cm
            get_mirror_photons = False # If True, get the photons from the mirror
            plotting = False # If True, plot the optics
            x_shower, y_shower, wavelength, zem, T, _, cos_X, cos_Y = simulate_davies_cotton_optics(
                x_shower, y_shower, wavelength, zem, T, cos_X, cos_Y, zenith=np.degrees(zenith_rad),
                R=tele_radi, r_cam=cam_radi,
                get_mirror_photons=get_mirror_photons,
                plotting=plotting,
            )
            print(f"Number of photons in event after optics only: {len(x_shower)} of {len(X)} Ratio: {len(x_shower)/len(X)}")
        if filter==True:
            x_shower, y_shower, T, wavelength, zem, mask = filter_photons(x_shower, y_shower, T, wavelength, zem, zenith_rad, atmos)
            cos_X = cos_X[mask]
            cos_Y = cos_Y[mask]
            print("Number of photons after filtering only:", len(x_shower), "of", len(X), "Ratio:", len(x_shower)/len(X))
        x_shower_new = x_shower
        y_shower_new = y_shower
        time_new = T
        wavelength_new = wavelength
        zem_new = zem - 294400 # Adjust zem to the observation level
        cos_X_new = cos_X
        cos_Y_new = cos_Y
        # Append the new event to the list
        Event.append({
            'x_shower': x_shower_new,
            'y_shower': y_shower_new,
            'time': time_new,
            'wavelength': wavelength_new,
            'zem': zem_new,
            'cos_x': cos_X_new,
            'cos_y': cos_Y_new,
            'azimuth': azimuth_rad,
            'zenith': zenith_rad,
            'x_brightest': x_brightest,
            'y_brightest': y_brightest
        })
    return Event

def extract_data_from_file(input_file, atmos, R=150, tele_id=0, filter=False, optics=False, get_brightest_pixel=False, shower_coordinates=True, zenith_tele=91, azimuth_tele=280):
    x_shower_new = []
    y_shower_new = []
    wavelength_new = []
    zem_new = []
    cos_x = []
    cos_y = []
    x_brightest = []
    y_brightest = []
    events = get_photon_bunches(input_file, atmos, R=R, 
                                filter=filter, 
                                optics=optics, 
                                get_brightest_pixel=get_brightest_pixel, 
                                tele_id=tele_id,
                                shower_coordinates=shower_coordinates,
                                zenith_tele=zenith_tele,
                                azimuth_tele=azimuth_tele)
    for event in events:
        x_shower_new.extend(event['x_shower'])
        y_shower_new.extend(event['y_shower'])
        wavelength_new.extend(event['wavelength'])
        zem_new.extend(event['zem'])
        time_new = event['time']
        cos_x.extend(event['cos_x'])
        cos_y.extend(event['cos_y'])
        azimuth_event_rad = event['azimuth']
        azimuth_event_degree = np.degrees(azimuth_event_rad)
        zenith_event_rad = event['zenith']
        zenith_event_degree = np.degrees(zenith_event_rad)
        x_brightest.append(event['x_brightest'])
        y_brightest.append(event['y_brightest'])
    x_shower_new = np.array(x_shower_new)
    y_shower_new = np.array(y_shower_new)
    cos_x = np.array(cos_x)
    cos_y = np.array(cos_y)
    wavelength_new = np.array(wavelength_new)
    zem_new = np.array(zem_new)
    time_new = np.array(time_new)
    x_brightest = np.array(x_brightest)
    y_brightest = np.array(y_brightest)
    return time_new, x_shower_new, y_shower_new, wavelength_new, zem_new, cos_x, cos_y, azimuth_event_degree, zenith_event_degree, x_brightest, y_brightest

def simulate_davies_cotton_optics(X_photon, Y_photon, W_photon, Zem_photon, T_photon, cos_x, cos_y, zenith, R = 150, r_cam = 4.96, get_mirror_photons=False, plotting=False):
    # 1) compute z on the spherical dish, keep only real hits
    r2 = X_photon**2 + Y_photon**2
    ok = (r2 <= R*R)
    X = X_photon[ok]
    Y = Y_photon[ok]
    W = W_photon[ok] # Wavelength of the photons 
    Zem = Zem_photon[ok]
    T = T_photon[ok] # Time since the first interaction
    # Distance from the mirror surface to the camera plane
    Z = np.sqrt(R*R - r2[ok])
    # 2) incident unit vectors
    # assume v_i.z < 0  i.e. photon travelling towards negative z
    cos_z = np.sqrt(1.0 - cos_x[ok]**2 - cos_y[ok]**2) #cosine of the z direction
    v_i = np.vstack((cos_x[ok], cos_y[ok], -cos_z)).T # incident unit vectors, shape (M,3) Pointing towards the mirror surface
    
    # 3) normals at each hit
    # center of curvature C = (0,0,R)
    # P_hit = (X,Y,Z)
    # so n = (P_hit - C)/R
    #     = (X, Y, Z-R)/R
    R_norm  = R*2 # David-Cotton optics have a radius of curvature of 2*R for the focal point to be at the camera plane
    Z2 = np.sqrt(R_norm**2 - r2[ok]) #Z2 is the z coordinate of the normal vector that point at 2R
    #print(R_norm - Z2)
    n = np.vstack((-X, -Y, R_norm-(-Z2 + R_norm))).T  # shape (M,3) normal vectors of the mirror surface with Davies-Cotton optics
    n = n/R_norm
    #print("Check normal vector", n)

    #normalize the normal vectors
    #n = n / np.linalg.norm(n, axis=1)[:, None]  # normalize each normal vector
    # 4) reflected vectors
    # v_r = v_i - 2*(v_i·n)*n
    dot_in = np.sum(v_i * n, axis=1) # dot product of incident vector and normal
    v_r = v_i - 2*dot_in[:,None]*n # reflected vectors, shape (M,3)
    n_z = np.vstack((0, 0, -1)).T
    #print("vi",v_i, np.degrees(np.arccos((np.sum(v_i * n_z, axis=1)))))
    #print("vr",v_r, np.degrees(np.arccos((np.sum(-v_i * n, axis=1)))), np.degrees(np.arccos((np.sum(v_r * n, axis=1)))))
    # 5) intersect with camera plane z=0
    #   P_img = P_hit + t * v_r, solve Z + t*v_r_z = 0  =>  t = -Z / v_r_z
    t = Z / v_r[:,2] # time to intersection with camera plane
    #print("time",t, Z)
    t_back = -Z/ v_i[:,2] # time to back-project incident photons to the camera plane
    #print("time back",t_back)
    # keep only forward‐going intersections
    good = (t > 0) & (t_back > 0)
    cos_z = cos_z[good]
    t = t[good]
    X_hit = X[good] # keep the coordinates of the photons that hit the camera
    Y_hit = Y[good] # keep the coordinates of the photons that hit the camera
    W = W[good] # keep the wavelength of the photons that hit the camera
    Zem = Zem[good] # keep the production height of the photons that hit the camera
    T = T[good] # keep the time of the photons that hit the camera
    #reflected vectors
    v_rx  = v_r[good,0] 
    v_ry  = v_r[good,1]
    #Get the incident angles of the photons that hit the camera
    n_z = np.vstack((0, 0, -1)).T # normal vector of the camera plane
    incident_angles = np.arccos((np.sum(-v_r[good] * n_z, axis=1)))  # angle between incident vector and normal vector
    #angle between n_z and n:
    angle = np.arccos((np.sum(-n * n_z, axis=1)))
    #print("check normal angle",np.degrees(angle))
    #print("check incident angle", np.degrees(incident_angles))
    #Back‐projection of incident photons
    t_back = t_back[good]
    X_back = X[good]
    Y_back = Y[good]
    v_bx  = v_i[good,0]
    v_by  = v_i[good,1]

    # image coords
    x_img = X_hit + t * v_rx
    y_img = Y_hit + t * v_ry
    #R_hit = np.sqrt(x_img**2 + y_img**2)  # radius of the hit point on the camera
    #print("R_hit", R_hit)
    x_img_back = X_back - t_back * v_bx
    y_img_back = Y_back - t_back * v_by
    #R_back = np.sqrt(x_img_back**2 + y_img_back**2)  # radius of the back-projected point on the camera
    #print("R_back", R_back)
    # 6) select those within camera radius
    #inside = ((x_img**2 + y_img**2) <= r_cam*r_cam) 
    inside = (np.abs(x_img) <= r_cam) & (np.abs(y_img) <= r_cam)
    x_img = x_img[inside]
    y_img = y_img[inside]
    cos_z = cos_z[inside]  # Keep the cosine of the z direction of the photons that hit the camera
    W = W[inside]  # Keep the wavelength of the reflective photons that hit the camera
    Zem = Zem[inside]  # Keep the production height of the reflective photons that hit the camera
    T = T[inside]  # Keep the time of the reflective photons that hit the camera
    incident_angles = incident_angles[inside]  # Keep the incident angles of the photons that hit the camera
    y_img_back = y_img_back[inside]
    x_img_back = x_img_back[inside]
    # 6.2) Also check the back‐projection of incident photons that hit the camera
    #outside = ((x_img_back**2 + y_img_back**2) >= r_cam*r_cam)
    outside = (np.abs(x_img_back) > r_cam) | (np.abs(y_img_back) > r_cam)
    x_img = x_img[outside]
    y_img = y_img[outside]
    cos_z = cos_z[outside]  # Keep the cosine of the z direction of the incident photons that does not hit the camera
    W = W[outside]  # Keep the wavelength of incident photons that does not hit the camera
    Zem = Zem[outside]  # Keep the production height of incident photons that does not hit the camera
    T = T[outside]  # Keep the time of incident photons that does not hit the camera
    # 7) compute time of arrival at the camera
    c = 299792458  # Speed of light in m/s
    T_hit = (((Zem/100)/(c*np.cos(np.radians(zenith))*(cos_z))) + R/c)*1e9   # time it takes for the photons to arrive at the camera in nano seconds
    time = T_hit + T # Time photon arrive at the detector since the first interaction
    incident_angles = incident_angles[outside]  # Keep the incident angles of the photons that does not hit the camera
    #print(f"Photons hitting camera: {len(x_img)} of {len(X_photon)}")
    # 6.3) Plot the distribution of zenith angles of photons hitting the camera
    if plotting == True:
        #photon radius vs camera radius
        plt.figure(figsize=(10, 6))
        R_hit = np.sqrt(X_hit[inside][outside]**2 + Y_hit[inside][outside]**2)  # radius of the hit point on the camera
        r_hit = np.sqrt(x_img**2 + y_img**2)  # radius of the back-projected point on the camera
        plt.hist2d(r_hit, R_hit, bins=100, cmap='Blues', norm=LogNorm())
        plt.xlabel('Radius of Hit Point on Camera (cm)')
        plt.ylabel('Radius of Back-Projected Point on the mirror (cm)')
        plt.title('Radius of Hit Point on Mirror vs Hit Point on Camera')
        plt.show()
        theta_z = np.arccos(cos_z)  # angle between the z direction and the incident vector
        plt.figure(figsize=(10, 6))
        plt.hist(np.degrees(theta_z), bins=50, color='orange', alpha=0.5)
        plt.xlabel('Zenith Angle (degrees)')
        plt.ylabel('Number of Photons')
        plt.title('Distribution of Zenith Angles of Photons Hitting the Camera')
        plt.grid(True)
        plt.show()
        #plot incident angles distribution
        plt.figure(figsize=(10, 6))
        plt.hist(np.degrees(incident_angles), bins=50, color='blue', alpha=0.5)
        plt.xlabel('Incident Angle (degrees)')
        plt.ylabel('Number of Photons')
        plt.title('Distribution of Incident Angles of Photons Hitting the Camera')
        plt.grid(True)
        plt.show()
        #plot zenith vs camera radius
        plt.figure(figsize=(10, 6))
        plt.hist2d(np.degrees(theta_z), np.sqrt(x_img**2 + y_img**2), bins=100, cmap='Blues', norm=LogNorm())
        plt.colorbar(label='Number of Photons')
        plt.xlabel('Zenith Angle (degrees)')
        plt.ylabel('Camera Radius (cm)')
        plt.title('Zenith Angle vs Camera Radius of Photons Hitting the Camera')
        plt.show()
        #plot wavelength distribution
        plt.figure(figsize=(10, 6))
        plt.hist(W, bins=50, color='green', alpha=0.5, label='Photons Hitting Camera')
        plt.hist(W_photon, bins=50, color='blue', alpha=0.5, label='All Photons')
        plt.yscale('symlog')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Number of Photons')
        plt.title('Distribution of Wavelengths of Photons Hitting the Camera')
        plt.grid(True)
        plt.legend()
        plt.show()
        #plot Zem distribution
        plt.figure(figsize=(10, 6))
        plt.hist(Zem/1e2, bins=50, range=[0,10], color='red', alpha=0.5, label='Photons Hitting Camera')
        plt.hist(Zem_photon/1e2, range=[0,10], bins=50, color='blue', alpha=0.5, label='All Photons')
        plt.xscale('symlog')
        plt.yscale('symlog')
        plt.xlabel('Zem (m)')
        plt.ylabel('Number of Photons')
        plt.title('Distribution of Zem of Photons Hitting the Camera')
        plt.grid(True)
        plt.legend()
        plt.show()
        #Plot ZEM vs wavelength histogram 2D
        xmin = 1     # avoid zero
        xmax = np.max(Zem/1e2)
        nbins_x = 100
        xedges = np.logspace(np.log10(xmin),np.log10(xmax),nbins_x + 1)
        # Keep linear bins for W
        nbins_y = 100
        yedges = np.linspace(W.min(), W.max(), nbins_y + 1)
        plt.figure(figsize=(10, 6))
        plt.hist2d(Zem/1e2, W,bins=[xedges, yedges],cmap='inferno',norm=LogNorm())
        plt.xscale('log')
        plt.colorbar(label='Number of Photons')
        plt.xlabel('Zem (m)')
        plt.ylabel('Wavelength (nm)')
        plt.title('Zem vs Wavelength of Photons Hitting the Camera')
        plt.show()
    if get_mirror_photons==True:
        # 7) return the photons that hit the mirror
        X_hit_mirror = (X_hit[inside])[outside]
        Y_hit_mirror = (Y_hit[inside])[outside]
        cos_x_mirror = cos_x[ok][good][inside][outside]
        cos_y_mirror = cos_y[ok][good][inside][outside]
        return X_hit_mirror, Y_hit_mirror, W, Zem, time, incident_angles, cos_x_mirror, cos_y_mirror
    else:
        cos_x_img = cos_x[ok][good][inside][outside]  # Keep the cosine of the x direction of the photons that hit the camera
        cos_y_img = cos_y[ok][good][inside][outside]  # Keep the cosine of the y direction of the photons that hit the camera
        # Return the coordinates of the photons that hit the camera
        return x_img, y_img, W, Zem, time, incident_angles, cos_x_img, cos_y_img

def make_2d_histogram_consistent(x, y, x_edges, y_edges):
    """
    Create 2D histogram with consistent bin edges
    """
    if len(x) == 0 or len(y) == 0:
        # Return empty histogram with correct shape
        return np.zeros((len(y_edges)-1, len(x_edges)-1))
    
    heatmap, _, _ = np.histogram2d(x, y, bins=[x_edges, y_edges])
    return heatmap.T  # Transpose to match the orientation of the data

def Camera_Response(T_photon, x_base, y_base, w, noise=True):
    if len(T_photon) > 0:
        T_photon = T_photon - T_photon.min()  # Normalize time to start at zero
        delta_t = 10 # ns
        N = int(T_photon.max() / delta_t) + 1
        t_min = 0
        t_max = T_photon.max()
        n_window = int(np.ceil((t_max - t_min) / delta_t))
        # Pre-allocate lists
        T = [] # List to store time windows
        X = [] # List to store x coordinates
        Y = [] # List to store y coordinates
        W = [] # List to store wavelengths

        # Use numpy's digitize for efficient binning
        time_edges = np.arange(t_min, t_max + delta_t, delta_t)
        bin_indices = np.digitize(T_photon, time_edges) - 1

        # Group data by time windows
        for i in range(n_window):
            mask = (bin_indices == i)
            T.append(T_photon[mask])
            X.append(x_base[mask])
            Y.append(y_base[mask])
            W.append(w[mask])
            # Define global bin edges
        n_bins = 16
        d_cam = 9.92/2
        x_min, x_max = -d_cam, d_cam
        y_min, y_max =  -d_cam, d_cam
        # Add small padding to avoid edge effects
        x_padding = (x_max - x_min) * 0.01
        y_padding = (y_max - y_min) * 0.01

        # Create consistent bin edges for all histograms
        global_x_edges = np.linspace(x_min - x_padding, x_max + x_padding, n_bins + 1)
        global_y_edges = np.linspace(y_min - y_padding, y_max + y_padding, n_bins + 1)

        H = [] # List to store histograms
        Times = [] # List to store corresponding times
        center_times = time_edges[1:] - delta_t/2  # Fixed: should be delta_t/2, not 5
        for i in range(len(X)):
            h2d = make_2d_histogram_consistent(X[i], Y[i], global_x_edges, global_y_edges)
            # Apply your Trigger condition
            if (h2d >= 1).any(): #At least one bin must have a count >= 1
                H.append(h2d)
                Times.append(center_times[i])
        if len(H) > 0: 
            # Convert to numpy arrays for easier manipulation
            H = np.array(H) 
            Times = np.array(Times)
            # Camera response function to weight the histograms
            h = H.copy()
            for i in range(len(Times)):
                weight = camera_response_function(np.array(Times), Times[i])
                h += weight[:, np.newaxis, np.newaxis] * H[i]
            if noise == True:
                Noise = noise_function(h, lambda_fit=0.74)  # Generate noise for the histograms
                h += Noise  # Add noises to the histograms
            else:
                Noise = None
            return h, H, Times, Noise #H is no camera response, h is with camera response
        else:
            print("No histograms created. Check your data and filter conditions.")
            return None, None, None, None
    else:
        print("No photons found in the data.")
        return None, None, None, None

def get_parameters(filename, line_number):
    """
    Get E, h, z, w from a specific line number (1-based, excluding header)
    """
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    # Skip header (line 0) and get the requested line
    data_line_index = line_number  # line_number is already 1-based for data
    
    if data_line_index >= len(lines):
        raise IndexError(f"Line {line_number} not found. File has {len(lines)-1} data lines.")
    
    # Parse the line
    line = lines[data_line_index].strip()
    energy, height, zenith, weight = line.split()
    
    return float(energy), float(height), float(zenith), float(weight)

def read_particle_file(file_path, particle_list):
    #print("Reading event num:", file_path[-2]+file_path[-1])
    particles = []
    with CorsikaParticleFile(file_path, parse_blocks=False) as f:
        for event in f:
            i = 0
            #print('Start reading particles')
            for pp in event.particles:
                particle_pid = np.abs(int(pp[0] // 1000))  # Determine particle ID
                if particle_pid in particle_list:  # Muon: 5 or 6, 
                    i += 1
                    #print(i)
                    if particle_pid in [5, 6, 75, 76]:  # Muon reach the ground: 5 + or 6 -, else 85 or 86
                        mass = 0.1056583745 #GeV
                    elif particle_pid == 1: #gamma 
                        mass = 0
                    elif particle_pid in [2, 3]: #electron (2) and positron (3)
                        mass = 0.000511 #GeV
                    elif particle_pid == 7: #pi0
                        mass = 0.1349766 #GeV
                    elif particle_pid in [8,9]: #pi+ and p-
                        mass = 0.13957015 #GeV
                    else:
                        print("Unknown particle ID:", particle_pid)
                        continue
                    energy = np.sqrt(mass**2 + pp[1]**2 + pp[2]**2 + pp[3]**2)
                    direction_vector = np.array([pp[1], pp[2], pp[3]])  # Momentum components
                    direction_vector /= np.linalg.norm(direction_vector)  # Normalize to get unit vector
                    #print("Particle ID", particle_pid)
                    #print("  x", pp[4], "cm", "  y", pp[5], "cm")
                    #print("  E", energy, "GeV","  h", pp[6], "cm")
                    particles.append([particle_pid, pp[4], pp[5], pp[6], energy, direction_vector])
    #p_list = np.array([int(p[0]) for p in particles])
    #count_repeats_numpy(p_list)
    return particles

def get_muon_info(file_path, azimuth_rad, zenith_deg=91, all=True, particle_list=[5,6, 75, 76]): #5 and 6 are muons, 75 and 76 are muons with additional info
    pp = read_particle_file(file_path, particle_list)  # Read only specified particles
    x = [p[1] for p in pp if p[0] in [5, 6]]  # Extract x positions of muons on the ground
    y = [p[2] for p in pp if p[0] in [5, 6]]  # Extract y positions of muons on the ground
    h = [p[3] for p in pp if p[0] in [75, 76]]  # Extract production heights of muons if 75 or 76 else is the time since first interaction
    #t = [p[3] for p in pp if p[0] in [5, 6]]  # Extract time since first interaction for muons if 5 or 6
    E = [p[4] for p in pp if p[0] in [5, 6]]  # Extract energies of muons on the ground
    cos_x = [p[5][0] for p in pp if p[0] in [5, 6]] # Extract cosine of x direction of muons on the ground
    cos_y = [p[5][1] for p in pp if p[0] in [5, 6]] # Extract cosine of y direction of muons on the ground
    cos_z = np.sqrt(1 - np.array(cos_x)**2 - np.array(cos_y)**2)
    theta_z = np.arccos(cos_z)
    theta_z_deg = np.degrees(theta_z)
    x, y, cos_x, cos_y = convert_coordinates_to_shower(np.array(x), np.array(y), np.array(cos_x), np.array(cos_y), -azimuth_rad, zenith_deg)
    try:
        if all == False:
            mask = np.array(E) == np.max(E)  # Select the muon with the highest energy
            x = np.array(x)[mask]
            y = np.array(y)[mask]
            cos_x = np.array(cos_x)[mask]
            cos_y = np.array(cos_y)[mask]
            E = np.array(E)[mask]  # Get the highest energy value
            h = np.array(h)[mask]  # Get the production height
            theta_z_deg = theta_z_deg[mask]
    except:
        x = np.array([np.nan])
        y = np.array([np.nan])
        cos_x = np.array([np.nan])
        cos_y = np.array([np.nan])
    return x, y, cos_x, cos_y, E, h, theta_z_deg
import time

def Monte_Carlos_Avg_air_shower(gen_radi = "150", noise=True, particle = 14, N = 1000, name_tag="1E3_1E7"):
      # telescope radius in cm, available:  1500 cm
    Store_Events = []  # Store the results of each Monte Carlo iteration
    Error_ID = []
    ran_seed = int(time.time())
    for i in range(N):
        atmabs_file="/uufs/chpc.utah.edu/common/home/u1520754/corsika-78010/run/atmabs.dat"
        ran_val = np.random.randint(0, ran_seed) + i*10
        print("Random seed:", ran_val)
        np.random.seed(ran_val) # Random seed for reproducibility
        run_tag = 5 #np.random.randint(0,5)
        atmos = VAtmosAbsorption(model="CORSIKA", seed=ran_val, source_file=atmabs_file)
        id = i #np.random.randint(0, 4999)  # Randomly select an ID for each iteration
        if i > 4999:
            id = np.random.randint(0, 4999)
        E_i, h_i, zenith_i, weight_i = get_parameters(f"/uufs/chpc.utah.edu/common/home/u1520754/corsika_inputs/sampled_data_{name_tag}_5000events.txt", id+1)
        input_files = f"/uufs/chpc.utah.edu/common/home/u1520754/corsika_results/notch_r{gen_radi}.0_Trinity_Demo_{name_tag}_{run_tag}/out_sib23d-pId{particle}/OUTPUT_{1000+id}/telescope.dat"
        if 1000+id < 10000:
            particle_file = f"/uufs/chpc.utah.edu/common/home/u1520754/corsika_results/notch_r{gen_radi}.0_Trinity_Demo_{name_tag}_{run_tag}/out_sib23d-pId{particle}/OUTPUT_{1000+id}/DAT00{1000+id}"
        else:
            particle_file = f"/uufs/chpc.utah.edu/common/home/u1520754/corsika_results/notch_r{gen_radi}.0_Trinity_Demo_{name_tag}_{run_tag}/out_sib23d-pId{particle}/OUTPUT_{1000+id}/DAT0{1000+id}"
        print("Processing event ID:", id, "run tag", run_tag)
        try:
            zenith_tele = 91
            azimuth_tele = 280
            T_photon, x_base, y_base, w, Zem, cos_x, cos_y, azimuth_event_degree, zenith_event_degree, x_brightest, y_brightest = extract_data_from_file(input_files, atmos, 
                                                                                                                            R=float(gen_radi), 
                                                                                                                            filter=False, 
                                                                                                                            optics=False,
                                                                                                                            get_brightest_pixel=False, #Get the brightest pixel location
                                                                                                                            shower_coordinates=True, #Use shower coordinates
                                                                                                                            zenith_tele=zenith_tele, 
                                                                                                                            azimuth_tele=azimuth_tele
                                                                                                                            )
            x_mu, y_mu, cos_x_mu, cos_y_mu, E_mu, h_mu, zenith_mu = get_muon_info(particle_file, azimuth_rad=np.radians(azimuth_tele), zenith_deg=zenith_tele)
            print("Check zenith",zenith_event_degree, zenith_i)
            # Apply the Davies-Cotton optics simulation
            
            print("Check data: height:", h_i, h_mu, "Energy:", E_i, E_mu, "Zenith:", zenith_i, zenith_mu)
            x_offset = 0 #np.random.uniform(-150, 150)
            y_offset = 0 #np.random.uniform(-150, 150)
            x_base += x_offset
            y_base += y_offset
            x_mu += x_offset
            y_mu += y_offset
            Photon_locations = [x_base.copy(),y_base.copy()]
            Photon_directions = [cos_x.copy(), cos_y.copy()]
            Muon_locations = [x_mu.copy(),y_mu.copy()]
            Muon_directions = [cos_x_mu.copy(), cos_y_mu.copy()]
            x_base, y_base, w, Zem, T_photon, incident_angles, cos_x_optic, cos_y_optic = simulate_davies_cotton_optics(
                        x_base, y_base, w, Zem, T_photon, cos_x, cos_y, np.radians(zenith_event_degree),
                        R=150, r_cam=4.96, get_mirror_photons=False,
                        plotting=False
                    )
            #Apply the atmospheric absorption and photon detection efficiency:
            x_base, y_base, T_photon, w, Zem, mask = filter_photons(x_base, y_base, T_photon, w, Zem, np.radians(zenith_event_degree), atmos, filter=True)
            cos_x_filter = cos_x_optic[mask]  # Keep the cosine of the x direction of the photons that hit the camera
            cos_y_filter = cos_y_optic[mask]  # Keep the cosine of the y direction of the photons that hit the camera
            T_photon = np.array(T_photon)
            # Call the camera response function
            Photon_camera = [x_base.copy(), y_base.copy()]
            hist, _, cam_time, _ = Camera_Response(T_photon, x_base, y_base, w, noise=noise)
            Camera = hist if hist is not None else 0 #Events is the histogram of the camera response
            Time_cam = cam_time if cam_time is not None else 0  # Time of the camera response
            Max_signal = np.max(hist) if hist is not None else 0  # Get the maximum signal from the histogram
        except Exception as e:
            print(f"Error occurred: {e}")
            print(f"Energy: {E_i}, Production Height: {h_i}, Zenith Angle: {zenith_i}")
            Max_signal = np.nan
            Time_cam = np.nan
            Camera = np.nan
            Photon_locations = np.nan
            Photon_directions = np.nan
            Muon_locations = np.nan
            Muon_directions = np.nan
            Photon_camera = np.nan
            cos_x = np.array([])
            cos_y = np.array([])
            zenith_event_degree = zenith_i
            azimuth_event_degree = np.nan
            x_mu = np.array([np.nan])
            y_mu = np.array([np.nan])
            E_mu = np.array([np.nan])
            h_mu = np.array([np.nan])
            Error_ID.append(id)
        print(" ")
        E = [E_i, E_mu]
        h = [h_i, h_mu]
        Store_Events.append([Max_signal, E, h, weight_i, zenith_event_degree, azimuth_event_degree, Muon_locations, Muon_directions, Photon_locations, Photon_directions, Time_cam, Camera, Photon_camera])  # Store data
    Error_ID = np.array(list(set(Error_ID)))  # Unique error IDs
    #sort Error_ID
    Error_ID.sort()
    print(f"Error IDs (if any): {Error_ID}")
    return Store_Events


def Monte_Carlos_Avg_muon(gen_radi = "150", noise=True, particle = 6, N = 1000, name_tag="1E3_1E7"):
      # telescope radius in cm, available:  1500 cm
    Store_Events = []  # Store the results of each Monte Carlo iteration
    Error_ID = []
    ran_seed = int(time.time())
    for i in range(N):
        atmabs_file="/uufs/chpc.utah.edu/common/home/u1520754/corsika-78010/run/atmabs.dat"
        ran_val = np.random.randint(0, ran_seed) + i*10
        print("Random seed:", ran_val)
        np.random.seed(ran_val) # Random seed for reproducibility
        run_tag = 0 #np.random.randint(0,4)
        atmos = VAtmosAbsorption(model="CORSIKA", seed=ran_val, source_file=atmabs_file)
        id = i #np.random.randint(0, 4999)  # Randomly select an ID for each iteration
        if i > 4999:
            id = np.random.randint(0, 4999)
        E_i, h_i, zenith_i, weight_i = get_parameters(f"/uufs/chpc.utah.edu/common/home/u1520754/corsika_inputs/sampled_data_{name_tag}_5000events.txt", id+1)
        input_files = f"/uufs/chpc.utah.edu/common/home/u1520754/corsika_results/notch_r{gen_radi}.0_Trinity_Demo_{name_tag}_{run_tag}/out_sib23d-pId{particle}/OUTPUT_{1000+id}/telescope.dat"
        if 1000+id < 10000:
            particle_file = f"/uufs/chpc.utah.edu/common/home/u1520754/corsika_results/notch_r{gen_radi}.0_Trinity_Demo_{name_tag}_{run_tag}/out_sib23d-pId{particle}/OUTPUT_{1000+id}/DAT00{1000+id}"
        else:
            particle_file = f"/uufs/chpc.utah.edu/common/home/u1520754/corsika_results/notch_r{gen_radi}.0_Trinity_Demo_{name_tag}_{run_tag}/out_sib23d-pId{particle}/OUTPUT_{1000+id}/DAT0{1000+id}"
        print("Processing event ID:", id, "run tag", run_tag)
        try:
            zenith_tele = 91
            azimuth_tele = 280
            T_photon, x_base, y_base, w, Zem, cos_x, cos_y, azimuth_event_degree, zenith_event_degree, x_brightest, y_brightest = extract_data_from_file(input_files, atmos, 
                                                                                                                            R=float(gen_radi), 
                                                                                                                            filter=False, 
                                                                                                                            optics=False,
                                                                                                                            get_brightest_pixel=False, #Get the brightest pixel location
                                                                                                                            shower_coordinates=True, #Use shower coordinates
                                                                                                                            zenith_tele=zenith_tele, 
                                                                                                                            azimuth_tele=azimuth_tele
                                                                                                                            )
            x_mu, y_mu, cos_x_mu, cos_y_mu, E_mu, h_mu, zenith_mu = get_muon_info(particle_file, azimuth_rad=np.radians(azimuth_tele), zenith_deg=zenith_tele)
            print("Check zenith",zenith_event_degree, zenith_i)
            # Apply the Davies-Cotton optics simulation
            
            print("Energy:", E_i, np.max(E_mu), "Zenith:", zenith_i, np.max(zenith_mu))
            x_offset = 0 #np.random.uniform(-150, 150)
            y_offset = 0 #np.random.uniform(-150, 150)
            x_base += x_offset
            y_base += y_offset
            x_mu += x_offset
            y_mu += y_offset
            Photon_locations = [x_base.copy(),y_base.copy()]
            Photon_directions = [cos_x.copy(), cos_y.copy()]
            Muon_locations = [x_mu.copy(),y_mu.copy()]
            Muon_directions = [cos_x_mu.copy(), cos_y_mu.copy()]
            x_base, y_base, w, Zem, T_photon, incident_angles, cos_x_optic, cos_y_optic = simulate_davies_cotton_optics(
                        x_base, y_base, w, Zem, T_photon, cos_x, cos_y, np.radians(zenith_event_degree),
                        R=150, r_cam=4.96, get_mirror_photons=False,
                        plotting=False
                    )
            #Apply the atmospheric absorption and photon detection efficiency:
            x_base, y_base, T_photon, w, Zem, mask = filter_photons(x_base, y_base, T_photon, w, Zem, np.radians(zenith_event_degree), atmos, filter=True)
            cos_x_filter = cos_x_optic[mask]  # Keep the cosine of the x direction of the photons that hit the camera
            cos_y_filter = cos_y_optic[mask]  # Keep the cosine of the y direction of the photons that hit the camera
            T_photon = np.array(T_photon)
            # Call the camera response function
            Photon_camera = [x_base.copy(), y_base.copy()]
            hist, _, cam_time, _ = Camera_Response(T_photon, x_base, y_base, w, noise=noise)
            Camera = hist if hist is not None else 0 #Events is the histogram of the camera response
            Time_cam = cam_time if cam_time is not None else 0  # Time of the camera response
            Max_signal = np.max(hist) if hist is not None else 0  # Get the maximum signal from the histogram
            if Max_signal >= 20:
                print("High signal event detected!")
        except Exception as e:
            print(f"Error occurred: {e}")
            print(f"Energy: {E_i}, Production Height: {h_i}, Zenith Angle: {zenith_i}")
            Max_signal = np.nan
            Time_cam = np.nan
            Camera = np.nan
            Photon_locations = np.nan
            Photon_directions = np.nan
            Muon_locations = np.nan
            Muon_directions = np.nan
            Photon_camera = np.nan
            cos_x = np.array([])
            cos_y = np.array([])
            zenith_event_degree = zenith_i
            azimuth_event_degree = np.nan
            x_mu = np.array([np.nan])
            y_mu = np.array([np.nan])
            E_mu = np.array([np.nan])
            h_mu = np.array([np.nan])
            Error_ID.append(id)
        print(" ")
        E = [E_i, E_mu]
        h = [h_i, h_mu]
        Store_Events.append([Max_signal, E, h, weight_i, zenith_event_degree, azimuth_event_degree, Muon_locations, Muon_directions, Photon_locations, Photon_directions, Time_cam, Camera, Photon_camera])  # Store data
    Error_ID = np.array(list(set(Error_ID)))  # Unique error IDs
    #sort Error_ID
    Error_ID.sort()
    print(f"Error IDs (if any): {Error_ID}")
    return Store_Events

def main():
    parser = argparse.ArgumentParser(description='Detector Sensitivity Monte Carlo Simulation')

    parser.add_argument('--particle', type=int, default=6, 
                        help='Particle ID (5 for mu+, 6 for mu-)')
    parser.add_argument('--gen_radi', type=str, default="150", 
                        help='Generation radius in cm')
    parser.add_argument('--N', type=int, default=1000, 
                        help='Number of Monte Carlo iterations')
    parser.add_argument('--noise', action='store_true', default=True,
                        help='Include noise in simulation')
    parser.add_argument('--output_dir', type=str, 
                        default='/uufs/chpc.utah.edu/common/home/u1520754/test_jupyter/outputs',
                        help='Output directory for results')
    parser.add_argument('--job_id', type=str, default='',
                        help='Job ID for output filename')
    parser.add_argument('--name_tag', type=str, default='1E3_1E7',
                        help='Name tag for output files')

    args = parser.parse_args()
    
    # Run simulation
    print(f"Running simulation with:")
    print(f"  Particle ID: {args.particle}")
    print(f"  Generation radius: {args.gen_radi}")
    print(f"  Number of iterations: {args.N}")
    print(f"  Noise: {args.noise}")
    
    # Run the Monte Carlo simulation
    
    Events_Data = Monte_Carlos_Avg_muon(
        gen_radi=args.gen_radi, 
        noise=args.noise, 
        particle=args.particle, 
        N=args.N,
        name_tag=args.name_tag
    )

    # Create output filename
    job_suffix = f"_job{args.job_id}" if args.job_id else ""
    output_filename = f'{args.output_dir}/Test_Events_Data_pid{args.particle}_N{args.N}_{args.name_tag}.pickle'

    # Save results
    with open(output_filename, 'wb') as f:
        pickle.dump(Events_Data, f)
    
    print(f"Results saved to: {output_filename}")
    print(f"Simulation completed successfully!")
if __name__ == "__main__":
    main()