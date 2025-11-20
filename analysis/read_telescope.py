"""
=============================================================================
    corsikaIOreader is a tool to read CORSIKA eventio files
    Copyright (C) 2004, 2013, 2019 Gernot Maier and Henrike Fleischhack

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
=============================================================================

VAtmosAbsorption - atmospheric absorption of Cherenkov photons

equiv. to CORSIKA function ATABSO or kascade function atm_pass
(routines are copied from these programs)

sourcefiles for extinction values are:
    - atmabs.dat   default CORSIKA data  
    - kextint.dat  default cherenkf dat
    - us76.50km.ext MODTRAN4 calculation of Michael Daniel (from 200nm)

assume extinction values are constant above 50 km
all units for input parameters are in meter and nanometer

Example usage:
    # use extinction values from MODTRAN4
    atabso = VAtmosAbsorption("us76_new", 0)
    # set observation height to 1500m
    atabso.set_observation_level(1500.0)
    # set wavelength interval in [nm]  
    atabso.set_wavelength_interval(185.0, 700.0)
    
    # get photon wavelength from emission height 5km and with direction cosinus 0.95
    wavelength = atabso.get_wavelength(5000.0, 0.95)
"""

import math
import sys
import os
from typing import Dict, List, Tuple


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

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from eventio import IACTFile
import scipy.stats as stats
from corsikaio import CorsikaParticleFile

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
            if not survive_photon(wavelength_new[j], zem_new[j] / 100, theta, atmos) or zem_new[j]/1e5 > 20:
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
    cos_y = y_i#np.sqrt(1 - cos_x**2 - cos_z**2)  # normalized y direction
    return x_shower, y_shower, cos_x, cos_y

def get_brightness_pixel_loc(all_x_shower, all_y_shower, x_range, y_range):
    # Find all locations with maximum value
    Num_bins = 50
    heatmap, xedges, yedges = np.histogram2d(
                    all_x_shower, 
                    all_y_shower, 
                    bins=int(Num_bins), 
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

def get_photon_bunches(input_file, atmos, R=1500, tele_id=0, filter=False, optics=False, use_brightness_pixel=False, shower_coordinates=True, zenith_tele=91, azimuth_tele=280):
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
            x_shower, y_shower, cos_X, cos_Y = convert_coordinates_to_shower(X, Y, cos_X, cos_Y, azimuth_rad=-np.radians(azimuth_tele), zenith_deg=zenith_tele)
        else:
            x_shower, y_shower, cos_X, cos_Y = X, Y, cos_X, cos_Y
        if use_brightness_pixel==True:
            x_brightness, y_brightness, brightness_value = get_brightness_pixel_loc(x_shower, y_shower, (-R, R), (-R, R))
            x_shower = x_shower - x_brightness  # Center the shower coordinates
            y_shower = y_shower - y_brightness  # Center the shower coordinates
            print(f"Brightness pixel location in cm: ({x_brightness:.2f}, {y_brightness:.2f}) with value {brightness_value:.2f}")
            x_brightness, y_brightness, brightness_value = get_brightness_pixel_loc(x_shower, y_shower, (-R/10, R/10), (-R/10, R/10))
            x_shower = x_shower - x_brightness  # Center the shower coordinates
            y_shower = y_shower - y_brightness  # Center the shower coordinates
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
            'zenith': zenith_rad
        })
    return Event

def extract_data_from_file(input_file, atmos, R=1500, tele_id=0, filter=True, optics=True, use_brightness_pixel=False, shower_coordinates=True, zenith_tele=91, azimuth_tele=280):
    x_shower_new = []
    y_shower_new = []
    wavelength_new = []
    zem_new = []
    cos_x = []
    cos_y = []
    events = get_photon_bunches(input_file, atmos, R=R, 
                                filter=filter, 
                                optics=optics, 
                                use_brightness_pixel=use_brightness_pixel, 
                                shower_coordinates=shower_coordinates, 
                                tele_id=tele_id,
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
    x_shower_new = np.array(x_shower_new)
    y_shower_new = np.array(y_shower_new)
    cos_x = np.array(cos_x)
    cos_y = np.array(cos_y)
    wavelength_new = np.array(wavelength_new)
    zem_new = np.array(zem_new)
    time_new = np.array(time_new)
    return time_new, x_shower_new, y_shower_new, wavelength_new, zem_new, cos_x, cos_y, azimuth_event_degree, zenith_event_degree

def simulate_davies_cotton_optics(X_photon, Y_photon, W_photon, Zem_photon, T_photon, cos_x, cos_y, zenith, R = 1500, r_cam = 60, get_mirror_photons=True, plotting=False):
    # 1) compute z on the spherical dish, keep only real hits
    r2 = X_photon**2 + Y_photon**2
    #print("Number of original photon", len(r2))
    ok = (r2 <= R*R)
    X = X_photon[ok]
    #print("Number of photons hit the mirror:", len(X))
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
    #print("number of photons reflected:", len(X_hit))
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
    R_hit = np.sqrt(x_img**2 + y_img**2)  # radius of the hit point on the camera
    #print("R_hit", R_hit)
    x_img_back = X_back - t_back * v_bx
    y_img_back = Y_back - t_back * v_by
    R_back = np.sqrt(x_img_back**2 + y_img_back**2)  # radius of the back-projected point on the camera
    #print("R_back", R_back)
    # 6) select those within camera radius
    inside = ((x_img**2 + y_img**2) <= r_cam*r_cam) 
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
    outside = ((x_img_back**2 + y_img_back**2) >= r_cam*r_cam)
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

def count_repeats_numpy(arr):
    # Convert to numpy array
    arr = np.array(arr)
    
    # Get unique elements and counts
    unique, counts = np.unique(arr, return_counts=True)
    
    # Combine into readable format
    count_info = dict(zip(unique, counts))
    
    # Print all counts
    print("All elements:")
    for num, count in count_info.items():
        print(f"Number {num} appears {count} times")  

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
                    #print("  E", energy, "GeV","  t", pp[6], "s")
                    pp_x = pp[4]
                    pp_y = pp[5]
                    pp_t = pp[6]
                    particles.append([particle_pid, pp_x, pp_y, pp_t, energy, direction_vector])
    #p_list = np.array([int(p[0]) for p in particles])
    #count_repeats_numpy(p_list)
    return particles

def plot_photon_heatmap(input_file, particle_file, ax, event_index=0, telescope_id=None, theta=0, bins=500, tele_radi=150, scale="log"):
    """
    Plot the spatial distribution of Cherenkov photons as a heatmap with different particles.
    
    Parameters:
    -----------
    input_file : str
        Path to the CORSIKA telescope data file (.dat)
    particle_file : str
        Path to the CORSIKA particle file to get particle positions
    ax : matplotlib.axes.Axes
        The Axes object to plot on
    event_index : int
        Index of the event to analyze (0-based)
    telescope_id : int or None
        ID of specific telescope to analyze, or None for all telescopes
    theta : float
        Zenith angle in degrees
    bins : int or tuple
        Number of bins for histogram (resolution of heatmap)
    tele_radi : float
        Telescope radius in cm
    scale : str
        Color scale type: "log" or "linear"
    """
    
    # Validate scale parameter
    if scale not in ["log", "linear"]:
        print(f"Warning: scale='{scale}' not recognized. Using 'log' instead.")
        scale = "log"
    
    try:
        # First read particle data - include all particle types of interest
        particles = read_particle_file(
            particle_file, 
            [1, 2, 3, 5, 6, 7, 8, 9, 75, 76, 132]  # All particle types we want to visualize
        )
        
        # Process photon data as before
        with IACTFile(input_file) as f:
            events = iter(f)
            event = next(events)
            azimuth = (event.header['azimuth'])
            # Skip to requested event
            for _ in range(event_index):
                try:
                    event = next(events)
                except StopIteration:
                    print(f"Error: Event {event_index} not found in file {input_file}.")
                    return

            # Get telescope positions and photons
            if telescope_id is not None:
                try:
                    photons = [event.photon_bunches[telescope_id]]
                    positions = [f.telescope_positions[telescope_id]]
                except KeyError:
                    print(f"Error: Telescope {telescope_id} not found in event {event_index}.")
                    return
            else:
                photons = list(event.photon_bunches.values())
                positions = f.telescope_positions
            
            # Collect all photon positions
            all_x_shower = []
            all_y_shower = []
            
            total_photons = 0
            for pos, tel_photons in zip(positions, photons):
                n_photons = len(tel_photons)
                total_photons += n_photons

                if n_photons > 0:
                    X = tel_photons['x'] #+ pos[0]
                    Y = tel_photons['y'] #+ pos[1]
                    T = tel_photons['time']  # Time in nanoseconds
                    cos_x = tel_photons['cx']  # Cosine of the x direction
                    cos_y = tel_photons['cy']  # Cosine of the y direction
                    #print(pos)
                    Zem = tel_photons['zem']  # Adjust for telescope height
                    wavelength = tel_photons['wavelength']
                    # Project onto shower coordinate system
                    azimuth_rad = -azimuth
                    theta_rad = np.radians(theta)
                    atmabs_file="/home/teikiet/corsika-78010/run/atmabs.dat"
                    atmos = atmos = VAtmosAbsorption(model="CORSIKA", seed=12345, source_file= atmabs_file)
                    X, Y, T, wavelength, Zem, mask = filter_photons(X, Y, T, wavelength, Zem, theta_rad, atmos)
                    Zem = Zem - 294400  # Adjust Zem to the observation level (in cm)
                    x_shower = (np.array(X)*np.cos(azimuth_rad) - np.array(Y)*np.sin(azimuth_rad))*np.cos(theta_rad)
                    y_shower = np.array(X)*np.sin(azimuth_rad) + np.array(Y)*np.cos(azimuth_rad)
                    # cetering the shower coordinates
                    #x_shower = x_shower - tele_radi*np.tan(theta_rad)*np.cos(theta_rad)*np.cos(-azimuth_rad)  
                    #y_shower = y_shower - tele_radi*np.tan(theta_rad)*np.cos(theta_rad)*np.sin(-azimuth_rad) 
                    all_x_shower.extend(x_shower/100) #convert to meters
                    all_y_shower.extend(y_shower/100) #convert to meters
                    

            # Create the heatmap using histogram2d
            if total_photons > 0:
                # Define the bins
                x_range = (-tele_radi*1.2/100, tele_radi*1.2/100) # Convert to meters
                y_range = (-tele_radi*1.2/100, tele_radi*1.2/100) # Convert to meters
                # Create 2D histogram
                heatmap, xedges, yedges = np.histogram2d(
                    all_x_shower, 
                    all_y_shower, 
                    bins=bins, 
                    range=[x_range, y_range]
                )
                
                # Configure normalization and colorbar based on scale
                if scale == "log":
                    # For log scale, handle zeros by setting minimum value
                    heatmap_display = np.where(heatmap == 0, 0, heatmap)  # Keep zeros as zeros for masking
                    vmin = max(1, heatmap[heatmap > 0].min()) if heatmap[heatmap > 0].size > 0 else 1
                    vmax = heatmap.max()
                    
                    norm = LogNorm(vmin=vmin, vmax=vmax)
                    colorbar_label = 'Log₁₀(Photon Density)'
                    
                    # Create the heatmap plot with LogNorm
                    im = ax.imshow(
                        heatmap.T,  # Transpose for correct orientation
                        origin='lower',
                        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                        aspect='auto',
                        cmap='inferno',  # Good colormaps: 'viridis', 'plasma', 'inferno'
                        norm=norm
                    )
                    
                else:  # linear scale
                    vmin = 0
                    vmax = heatmap.max()
                    norm = None  # Use default linear normalization
                    colorbar_label = 'Photon Density'
                    
                    # Create the heatmap plot with linear scale
                    im = ax.imshow(
                        heatmap.T,  # Transpose for correct orientation
                        origin='lower',
                        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                        aspect='auto',
                        cmap='inferno',  # Good colormaps: 'viridis', 'plasma', 'inferno'
                        vmin=vmin,
                        vmax=vmax
                    )
                
                # Add colorbar with appropriate configuration
                cbar = plt.colorbar(im, ax=ax)
                cbar.set_label(colorbar_label, rotation=270, labelpad=15, fontsize=12, color='white')
                
                # Configure colorbar ticks based on scale
                if scale == "log":
                    # Use logarithmic ticks
                    log_ticks = np.logspace(np.log10(vmin), np.log10(vmax), 5)
                    cbar.set_ticks(log_ticks)
                    # Format as scientific notation for log scale
                    cbar.set_ticklabels([f'{int(np.log10(t))}' for t in log_ticks])
                else:
                    # Use linear ticks - matplotlib handles this automatically
                    pass
                
                plt.setp(cbar.ax.get_yticklabels(), color='white', fontsize=17)

                # Add contour lines to highlight structure (optional)
                #contour_levels = np.logspace(0, np.log10(heatmap.max()), 5)
                #ax.contour(
                #    heatmap.T, 
                #    levels=contour_levels,
                #    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                #    colors='white',
                #    alpha=0.3,
                #    linewidths=0.5
                #)
                
                # Define particle visualization settings
                particle_settings = {
                    1: {'marker': 'o', 'color': 'white', 'label': 'Gamma', 'size_scale': 100, 'alpha': 0.2},  # gamma
                    2: {'marker': '^', 'color': 'white', 'label': 'Electron', 'size_scale': 100, 'alpha': 0.5},  # electron
                    3: {'marker': '^', 'color': 'yellow', 'label': 'Positron', 'size_scale': 100, 'alpha': 0.5},  # positron
                    5: {'marker': 's', 'color': 'white', 'label': 'μ+', 'size_scale': 100, 'alpha': 0.5},  # mu+
                    6: {'marker': 's', 'color': 'yellow', 'label': 'μ-', 'size_scale': 100, 'alpha': 0.5},  # mu-
                    7: {'marker': '*', 'color': 'white', 'label': 'π0', 'size_scale': 100, 'alpha': 0.5},  # pi0
                    8: {'marker': 'h', 'color': 'white', 'label': 'π+', 'size_scale': 100, 'alpha': 0.5},  # pi+
                    9: {'marker': 'h', 'color': 'yellow', 'label': 'π-', 'size_scale': 100, 'alpha': 0.5},  # pi-
                    75: {'marker': 'D', 'color': 'white', 'label': 'μ+ (decay)', 'size_scale': 100, 'alpha': 0.5},  # mu+ decay
                    76: {'marker': 'D', 'color': 'yellow', 'label': 'μ- (decay)', 'size_scale': 100, 'alpha': 0.5},  # mu- decay
                    132: {'marker': 'o', 'color': 'cyan', 'label': 'tau-', 'size_scale': 100, 'alpha': 0.5}  # tau-
                }
                
                # Group particles by type
                particle_groups = {}
                for p in particles:
                    pid = p[0]
                    if pid not in particle_groups:
                        particle_groups[pid] = []
                    particle_groups[pid].append(p)
                
                # Plot each particle type with proper transformation
                legend_handles = []
                particle_counts = {}
                #Get the telescope position
                tele_x = positions[0][0]/100  # Convert to meters
                tele_y = positions[0][1]/100  # Convert to meters
                for pid, particle_list in particle_groups.items():
                    if pid in particle_settings:
                        settings = particle_settings[pid]
                        
                        # Prepare arrays for this particle type
                        p_x = []
                        p_y = []
                        p_energy = []
                        p_time = []
                        # Apply coordinate transformation
                        for p in particle_list:
                            # Extract raw coordinates
                            raw_x = p[1]/100 - tele_x # x position #convert to meters
                            raw_y = p[2]/100 - tele_y # y position #convert to meters
                            
                            # Adjust for telescope height
                            raw_x = raw_x - tele_radi/100 * np.tan(np.radians(theta))*np.cos(-azimuth_rad)  
                            raw_y = raw_y - tele_radi/100 * np.tan(np.radians(theta))*np.sin(-azimuth_rad) 
                             
                            energy = p[4]  # energy GeV
                            raw_time = p[3] * 1e-9  # time in seconds (convert from nanoseconds)
                            # Apply the same transformation as for photons
                            azimuth_rad = -azimuth
                            theta_rad = np.radians(theta)
                            # Transform coordinates to shower coordinates
                            x_transformed = (raw_x*np.cos(azimuth_rad) - raw_y*np.sin(azimuth_rad))*np.cos(theta_rad)
                            y_transformed = raw_x*np.sin(azimuth_rad) + raw_y*np.cos(azimuth_rad)
                            # Check if particle is within plot range
                            if (x_range[0] <= x_transformed <= x_range[1]) and (y_range[0] <= y_transformed <= y_range[1]):
                                p_x.append(x_transformed)
                                p_y.append(y_transformed)
                                p_energy.append(energy)
                                p_time.append(raw_time)
                        
                        # Count visible particles
                        particle_counts[pid] = len(p_x)
                        
                        # Plot particles if any are visible
                        if p_x:
                            # Calculate sizes based on energy
                            #sizes = np.array(p_energy) * settings['size_scale']
                            if len(p_x) > 1000:
                                sizes = settings['size_scale']/10
                                alphas = settings['alpha']/2
                            elif len(p_x) < 5:
                                sizes = settings['size_scale']
                                alphas = settings['alpha']
                                print(f"{settings['label']} Position in cm",raw_x*100, raw_y*100, np.array(p_energy),"GeV", p_time, "s")
                            elif len(p_x) > 100:
                                sizes = settings['size_scale']/5
                                alphas = settings['alpha']/1.5
                            else:
                                sizes = settings['size_scale']
                                alphas = settings['alpha']
                            # Plot particle positions
                            scatter = ax.scatter(
                                p_x, p_y,
                                s=sizes,
                                c=settings['color'],
                                marker=settings['marker'],
                                alpha=alphas,
                                label=f"{settings['label']} ({len(p_x)})",
                                edgecolors='black' ,
                                linewidths=1
                            )
                            legend_handles.append(scatter)
                
                # Add legend only if particles were plotted
                if legend_handles:
                    ax.legend(
                        handles=legend_handles,
                        loc='upper right',
                        fontsize=8,
                        framealpha=0.7,
                        facecolor='black',
                        edgecolor='white',
                        labelcolor='white'
                    )
                
                # Format the plot
                ax.set_facecolor('black')
                ax.set_aspect('equal')
                ax.set_xlim(x_range)
                ax.set_ylim(y_range)
                
                # Add labels
                ax.set_xlabel('X position (m)', fontsize=12, color='white')
                ax.set_ylabel('Y position (m)', fontsize=12, color='white')
                
                # Update title to include total particles and scale info
                total_particles = sum(particle_counts.values()) if particle_counts else 0
                ax.set_title(f'Photons: {total_photons:,}  |  Particles: {total_particles:,} ({scale.capitalize()} Scale)', 
                             fontsize=14, color='white')

                # Set ticks and grid properties
                ax.tick_params(axis='both', colors='white', labelsize=10)
                ax.grid(True, alpha=0.2, color='white', linestyle='--')
                
                # Set spines color
                for spine in ax.spines.values():
                    spine.set_edgecolor('white')
            else:
                ax.text(0.5, 0.5, 'No photons detected', 
                        horizontalalignment='center',
                        verticalalignment='center',
                        transform=ax.transAxes,
                        color='white', fontsize=14)
                ax.set_facecolor('black')
                
    except FileNotFoundError:
        print(f"Error: File '{input_file}' or '{particle_file}' not found.")
    except Exception as e:
        print(f"Error processing files: {str(e)}")
        import traceback
        traceback.print_exc()

def plot_4_events_with_ground_particles(ID, particle, zenith, E_mag, production_height, tele_radi, num_bins=500, scale="log"):
    """
    Plot 4 events with ground particles using specified color scale.
    
    Parameters:
    -----------
    ID : list
        List of event IDs to plot
    particle : int
        Particle type ID
    zenith : float
        Zenith angle in degrees
    E_mag : int
        Energy magnitude (power of 10)
    production_height : float
        Production height in cm
    tele_radi : float
        Telescope radius
    num_bins : int
        Number of bins for histogram
    scale : str
        Color scale type: "log" or "linear"
    """
    
    # Validate scale parameter
    if scale not in ["log", "linear"]:
        print(f"Warning: scale='{scale}' not recognized. Using 'log' instead.")
        scale = "log"
    
    # Create a 2x2 grid of subplots
    input_files = [
        f"/home/teikiet/notch_h{production_height}_r{tele_radi}./out_sib23d-pId{particle}-En1.00E{E_mag}-{zenith}./OUTPUT_1.00E{E_mag}_{particle}_{zenith}._{100+i}/telescope.dat"
        for i in ID
    ]
    
    particle_files = [
        f"/home/teikiet/notch_h{production_height}_r{tele_radi}./out_sib23d-pId{particle}-En1.00E{E_mag}-{zenith}./OUTPUT_1.00E{E_mag}_{particle}_{zenith}._{100+i}/DAT000{100+i}"
        for i in ID
    ]
    
    fig, axs = plt.subplots(2, 2, figsize=(14, 12), facecolor='black')
    
    # Set figure background to black for better visualization
    fig.patch.set_facecolor('black')
    
    # Determine particle name
    if particle == 6:
        particle_name = "μ⁻"
    elif particle == 5:
        particle_name = "μ⁺"
    elif particle == 132:
        particle_name = "τ⁻"
    else:
        particle_name = f"p{particle}"
        
    # Set title with white text on black background - include scale info
    title = f'Cherenkov Photon Density Map with Ground Particles ({scale.capitalize()} Scale)\nPrimary: {particle_name}, θ = {zenith}°, E = 10^{E_mag} GeV, h = {production_height} cm, r = {int(float(tele_radi)/100)} m'
    fig.suptitle(title, fontsize=16, color='white', y=0.98)

    # Flatten the 2D array of axes for easy iteration
    axs = axs.flatten()

    # Plot each file in a separate subplot
    for i, (input_file, particle_file, ax) in enumerate(zip(input_files, particle_files, axs)):
        print(f"Plotting event {ID[i]}:")
        print(f"with {scale} scale: {input_file} with particles from {particle_file}")
        plot_photon_heatmap(
            input_file, 
            particle_file, 
            ax, 
            event_index=0, 
            telescope_id=None, 
            theta=zenith, 
            bins=num_bins, 
            tele_radi=int(float(tele_radi)),
            scale=scale  # Pass the scale parameter
        )
        ax.set_title(f'Event {ID[i]}', fontsize=14, color='white')

    # Adjust layout with more space for title
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    
    # Add description at the bottom of the figure
    fig.text(
        0.5, 0.01, 
        "",
        ha='center', 
        color='white', 
        fontsize=12
    )
    
    # Save the figure with high resolution - include scale in filename
    #plt.savefig(
    #    f'photon_heatmap_with_particles_{scale}_z{zenith}_1E{E_mag}_p{particle}_h{production_height}_r{tele_radi}.png', 
    #    dpi=300, 
    #    bbox_inches='tight', 
    #    facecolor='black'
    #)
    #plt.show()




"""Example usage:
x_shower_new = []
y_shower_new = []
wavelength_new = []
zem_new = []
for i in range(len(input_files)):
    events = get_photon_bunches([input_files[i]])
    for event in events:
        x_shower_new.extend(event['x_shower'])
        y_shower_new.extend(event['y_shower'])
        wavelength_new.extend(event['wavelength'])
        zem_new.extend(event['zem'])
import matplotlib.pyplot as plt
import numpy as np
#log scale for the plots
def plot_results(x_shower_new, y_shower_new, wavelength_new, zem_new):
    # Set the style for dark background
    plt.style.use('dark_background')
    x_shower_new = np.array(x_shower_new)
    y_shower_new = np.array(y_shower_new)
    wavelength_new = np.array(wavelength_new)
    zem_new = np.array(zem_new)

    # Create 2x2 subplot layout
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.patch.set_facecolor('black')  # Set figure background to black

    # Plot 1: Production height histogram
    axes[0, 0].hist(zem/1e5, bins=int(len(zem)/10), histtype='step', label='zem', color='cyan')
    axes[0, 0].hist(zem_new/1e5, bins=int(len(zem)/10), histtype='step', label='zem_new', color='orange')
    axes[0, 0].set_xlabel('Production height (km)', color='white')
    axes[0, 0].set_yscale('log')
    axes[0, 0].set_xscale('log')
    axes[0, 0].set_ylabel('Photon Counts', color='white')
    axes[0, 0].legend()
    axes[0, 0].set_title('Production Height Distribution', color='white')
    axes[0, 0].tick_params(colors='white')

    # Plot 2: Wavelength histogram
    axes[0, 1].hist(wavelength, bins=100, histtype='step', label='wavelength', color='cyan')
    axes[0, 1].hist(wavelength_new, bins=100, histtype='step', label='wavelength_new', color='orange')
    axes[0, 1].set_yscale('log')
    axes[0, 1].set_xscale('log')
    axes[0, 1].legend()
    axes[0, 1].set_xlabel('Wavelength (nm)', color='white')
    axes[0, 1].set_ylabel('Photon Counts', color='white')
    axes[0, 1].set_title('Wavelength Distribution', color='white')
    axes[0, 1].tick_params(colors='white')

    from matplotlib.colors import LogNorm

    # Plot 3: 2D histogram (R vs wavelength) normalized by 2π*R
    R = np.sqrt(x_shower_new**2 + y_shower_new**2)

    # Create custom 2D histogram
    Bins = 100
    R_bins = np.logspace(np.log10(R[R>0].min()), np.log10(R.max()), Bins+1)  # 50 bins in log space
    wavelength_bins = np.linspace(wavelength_new.min(), wavelength_new.max(), Bins+1)

    # Create 2D histogram
    counts, R_edges, wl_edges = np.histogram2d(R, wavelength_new, bins=[R_bins, wavelength_bins])

    # Calculate bin centers
    R_centers = (R_edges[:-1] + R_edges[1:]) / 2
    wl_centers = (wl_edges[:-1] + wl_edges[1:]) / 2

    # Normalize by 2π*R for each R bin
    R_mesh, wl_mesh = np.meshgrid(R_centers, wl_centers, indexing='ij')
    normalized_counts = counts / (2 * np.pi * R_mesh)

    # Replace zeros and infinities for log scale
    normalized_counts[normalized_counts == 0] = np.nan
    normalized_counts[~np.isfinite(normalized_counts)] = np.nan

    # Plot 
    im3 = axes[1, 0].pcolormesh(R_edges, wl_edges, normalized_counts.T, 
                            cmap='rainbow', norm=LogNorm())
    cbar3 = fig.colorbar(im3, ax=axes[1, 0])
    cbar3.set_label('Counts/(2π*R) (log scale)', rotation=270, labelpad=15, color='white')
    cbar3.ax.yaxis.set_tick_params(color='white')
    cbar3.ax.yaxis.set_ticklabels(cbar3.ax.yaxis.get_ticklabels(), color='white')
    axes[1, 0].set_xlabel('R (cm)', color='white')
    axes[1, 0].set_ylabel('Wavelength (nm)', color='white')
    axes[1, 0].set_xscale('log')
    axes[1, 0].set_title('2D Histogram: R vs Wavelength (normalized)', color='white')
    axes[1, 0].tick_params(colors='white')

    # Plot 4: 2D histogram (X vs Y)
    x_range = [x_shower_new.min(), x_shower_new.max()]
    y_range = [y_shower_new.min(), y_shower_new.max()]
    bins = 200  # You can adjust this

    # Create 2D histogram
    heatmap, xedges, yedges = np.histogram2d(
        x_shower_new, 
        y_shower_new, 
        bins=bins, 
        range=[x_range, y_range]
    )

    # Create the heatmap plot - using LogNorm for better visualization
    im4 = axes[1, 1].imshow(
        heatmap.T,  # Transpose for correct orientation
        origin='lower',
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
        aspect='auto',
        cmap='rainbow',  # Changed to match your previous plot
        norm=LogNorm(vmin=1, vmax=heatmap.max())  # Log scale for better visibility
    )

    # Add colorbar
    cbar4 = fig.colorbar(im4, ax=axes[1, 1])
    cbar4.set_label('Counts (log scale)', rotation=270, labelpad=15, color='white')

    # Option 2: Logarithmic scale ticks (if using LogNorm)
    vmin, vmax = im4.norm.vmin, im4.norm.vmax
    # Use 10^n format for better readability with logarithmic scale
    log_ticks = np.logspace(np.log10(vmin), np.log10(vmax), 5)
    cbar4.set_ticks(log_ticks)
    # Format as 10^n notation
    cbar4.set_ticklabels([f'$10^{{{np.round(np.log10(t),1)}}}$' for t in log_ticks])
    cbar4.ax.yaxis.set_tick_params(color='white')
    cbar4.ax.yaxis.set_ticklabels(cbar4.ax.yaxis.get_ticklabels(), color='white')

    axes[1, 1].set_xlabel('X (cm)', color='white')
    axes[1, 1].set_ylabel('Y (cm)', color='white')
    axes[1, 1].set_title('2D Histogram: X vs Y Coordinates', color='white')
    axes[1, 1].tick_params(colors='white')

    # Adjust layout to prevent overlap
    plt.tight_layout()
    plt.show()

plot_results(x_shower_new, y_shower_new, wavelength_new, zem_new)
"""



import plotly.graph_objects as go
import plotly.io as pio
import tempfile
import os
def save_plotly_animation_to_mp4(fig, filename="animation.mp4", fps=10):
    """
    Save a Plotly figure with frames as an MP4 video.
    
    Parameters:
    fig: Plotly figure with frames
    filename: output MP4 filename
    fps: frames per second for the video
    """
    import subprocess
    import shutil
    
    # Check if ffmpeg is available
    if not shutil.which('ffmpeg'):
        print("Error: ffmpeg is not installed. Please install ffmpeg to create MP4 videos.")
        return False
    
    # Create temporary directory for frames
    with tempfile.TemporaryDirectory() as temp_dir:
        frame_files = []
        
        # Export each frame as PNG
        for i, frame in enumerate(fig.frames):
            # Create a figure with just this frame's data
            temp_fig = go.Figure(data=frame.data)
            temp_fig.update_layout(fig.layout)
            
            # Remove slider and animation controls for cleaner frames
            temp_fig.update_layout(sliders=None, updatemenus=None)
            
            # Save frame as PNG
            frame_path = os.path.join(temp_dir, f"frame_{i:04d}.png")
            pio.write_image(temp_fig, frame_path, width=700, height=700)
            frame_files.append(frame_path)
            print(f"Exported frame {i+1}/{len(fig.frames)}")
        
        # Use ffmpeg to create MP4 from PNG frames
        ffmpeg_cmd = [
            'ffmpeg',
            '-y',  # overwrite output file
            '-r', str(fps),  # input frame rate
            '-i', os.path.join(temp_dir, 'frame_%04d.png'),  # input pattern
            '-c:v', 'libx264',  # video codec
            '-pix_fmt', 'yuv420p',  # pixel format for compatibility
            '-r', str(fps),  # output frame rate
            filename
        ]
        
        try:
            subprocess.run(ffmpeg_cmd, check=True, capture_output=True)
            print(f"Successfully created MP4: {filename}")
            return True
        except subprocess.CalledProcessError as e:
            print(f"Error creating MP4: {e}")
            return False

