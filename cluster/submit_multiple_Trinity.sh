#!/bin/bash

# process_production_heights.sh - A script to run CORSIKA simulations with multiple production heights

# Define the production heights to process (in cm)
production_heights=(5.E5) #5.E5 10.E5 15.E5 20.E5 25.E5 30.E5 35.E5 40.E5 

# Original script path
SUBMIT_SCRIPT="/uufs/chpc.utah.edu/common/home/u1520754/cluster/submit_single_Trinity.sh"

# Optional: Define particle types and energies
# If you want to override the defaults in the submit script
particle_types=(6)    # 6 = muon-, 5 = muon+, 133 = tau neutrino, 132 is tau-, 14 = proton
energy_values=(1.00E6) #in GeV
zenith=(90.)           # Fixed zenith angle in degrees
telescope_radius=(1500.) # telescope radius in cm origin: 150. cm, test: 15000. cm
azimuth_min=(280.)
azimuth_max=(280.)

# Process each production height
echo "Starting batch processing for multiple production heights..."
for height in "${production_heights[@]}"; do
    echo "==================================="
    echo "Processing production height: $height cm"
    echo "==================================="
    
    # Method 1: Modify the script using sed (temporary copy)
    TMP_SCRIPT=$(mktemp)
    cp "$SUBMIT_SCRIPT" "$TMP_SCRIPT"
    
    # Update production height value in the temporary script
    sed -i "s/^production_height=.*/production_height=$height/" "$TMP_SCRIPT"
    # Update zenith angle (fixed)
    sed -i "s/^zenith_tele=.*/zenith_tele=$zenith/" "$TMP_SCRIPT"
    # Update particle types and energies if provided
    sed -i "s/^particles=(.*)/particles=(${particle_types[@]})/" "$TMP_SCRIPT"
    # Update energies if provided
    sed -i "s/^energies=(.*)/energies=(${energy_values[@]})/" "$TMP_SCRIPT"
    # Update telescope radius if provided
    sed -i "s/^telescope_radius=.*/telescope_radius=${telescope_radius}/" "$TMP_SCRIPT"
    # Update azimuth angle in the temporary script
    sed -i "s/^azimuth_tele=.*/azimuth_tele=$azimuth/" "$TMP_SCRIPT"

    # Execute the modified script
    bash "$TMP_SCRIPT"
    
    # Clean up
    rm "$TMP_SCRIPT"
    
    echo "Completed processing for production height: $height cm"
    echo "Waiting 0 seconds before next height..."
    sleep 0
done

echo "Batch processing complete for all production heights."
