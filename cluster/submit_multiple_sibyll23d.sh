#!/bin/bash

# process_zenith_angles.sh - A script to run CORSIKA simulations with multiple zenith angles

# Define the zenith angles to process
zenith_angles=(0. 30. 60. 88.)

# Original script path
SUBMIT_SCRIPT="/uufs/chpc.utah.edu/common/home/u1520754/cluster/submit_single_sibyll23d.sh"

# Optional: Define particle types and energies
# If you want to override the defaults in the submit script
particle_types=(6)    # 6 = muon-, 5 = muon+, 133 = tau neutrino, 132 is tau-, 14 = proton
energy_values=(1.00E4) #in GeV
production_height=(20.E5) # production height in cm
telescope_radius=(1500.) # telescope radius in cm origin: 150. cm, test: 15000. cm
# Generate a random azimuth angle between -180 and 180 degrees
azimuth=(280)

# Process each zenith angle
echo "Starting batch processing for multiple zenith angles..."
for zenith in "${zenith_angles[@]}"; do
    echo "==================================="
    echo "Processing zenith angle: $zenith°"
    echo "==================================="
    
    # Method 1: Modify the script using sed (temporary copy)
    TMP_SCRIPT=$(mktemp)
    cp "$SUBMIT_SCRIPT" "$TMP_SCRIPT"
    
    # Update zenith value in the temporary script
    sed -i "s/^zenith=theta*/zenith=$zenith/" "$TMP_SCRIPT"
    # Update particle types and energies if provided
    sed -i "s/^particles=(.*)/particles=(${particle_types[@]})/" "$TMP_SCRIPT"
    # Update energies if provided
    sed -i "s/^energies=(.*)/energies=(${energy_values[@]})/" "$TMP_SCRIPT"
    # Update production height if provided
    sed -i "s/^production_height=.*/production_height=${production_height}/" "$TMP_SCRIPT"
    # Update telescope radius if provided
    sed -i "s/^telescope_radius=.*/telescope_radius=${telescope_radius}/" "$TMP_SCRIPT"
    # Update azimuth angle in the temporary script
    sed -i "s/^azimuth=.*/azimuth=$azimuth/" "$TMP_SCRIPT"
    # Execute the modified script
    bash "$TMP_SCRIPT"
    
    # Clean up
    rm "$TMP_SCRIPT"
    
    echo "Completed processing for zenith angle: $zenith°"
    echo "Waiting 0 seconds before next angle..."
    sleep 0
done

echo "Batch processing complete for all zenith angles."
