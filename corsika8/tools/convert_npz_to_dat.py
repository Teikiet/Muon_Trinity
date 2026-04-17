#!/usr/bin/env python3
"""
Convert CHASM abstable.npz files to text format for CORSIKA 8 Cherenkov absorption.

Usage:
    python convert_npz_to_dat.py input.npz output.dat
    
The input NPZ file should contain:
    - 'ecoeff': 2D array of optical depth τ(λ, h), shape (n_wavelengths, n_heights)
    - 'wavelength': 1D array of wavelengths in nm, shape (n_wavelengths,)
    - 'height': 1D array of heights in m, shape (n_heights,)

The output text file format is:
    Line 1: n_wavelengths n_heights
    Line 2: wavelength_0 wavelength_1 ... wavelength_(n-1)  [nm]
    Line 3: height_0 height_1 ... height_(m-1)  [m]
    Lines 4+: ecoeff values, one row per wavelength, columns are heights

Example:
    python convert_npz_to_dat.py abstable_MODTRAN_new.npz abstable_MODTRAN_new.dat
"""

import numpy as np
import sys
import os

def convert_npz_to_dat(input_file: str, output_file: str) -> None:
    """Convert NPZ absorption table to text format."""
    
    # Load the NPZ file
    print(f"Loading {input_file}...")
    data = np.load(input_file)
    
    # Extract arrays
    ecoeff = data['ecoeff']  # Shape: (n_wavelengths, n_heights)
    wavelengths = data['wavelength']  # Shape: (n_wavelengths,)
    heights = data['height']  # Shape: (n_heights,)
    
    n_wavelengths = len(wavelengths)
    n_heights = len(heights)
    
    print(f"  wavelength range: [{wavelengths.min():.1f}, {wavelengths.max():.1f}] nm ({n_wavelengths} values)")
    print(f"  height range: [{heights.min():.1f}, {heights.max():.1f}] m ({n_heights} values)")
    print(f"  ecoeff shape: {ecoeff.shape}")
    print(f"  ecoeff range: [{ecoeff.min():.6g}, {ecoeff.max():.6g}]")
    
    # Verify shapes match
    if ecoeff.shape != (n_wavelengths, n_heights):
        raise ValueError(f"ecoeff shape {ecoeff.shape} doesn't match expected ({n_wavelengths}, {n_heights})")
    
    # Write output file
    print(f"Writing {output_file}...")
    with open(output_file, 'w') as f:
        # Line 1: dimensions
        f.write(f"{n_wavelengths} {n_heights}\n")
        
        # Line 2: wavelengths
        f.write(" ".join(f"{w:.6g}" for w in wavelengths) + "\n")
        
        # Line 3: heights
        f.write(" ".join(f"{h:.6g}" for h in heights) + "\n")
        
        # Lines 4+: ecoeff values (one row per wavelength)
        for i in range(n_wavelengths):
            f.write(" ".join(f"{e:.6g}" for e in ecoeff[i]) + "\n")
    
    # Verify file size
    file_size = os.path.getsize(output_file)
    print(f"  Output file size: {file_size} bytes")
    print("Done!")


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        print("Error: Expected 2 arguments: input.npz output.dat")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(input_file):
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)
    
    convert_npz_to_dat(input_file, output_file)


if __name__ == "__main__":
    main()
