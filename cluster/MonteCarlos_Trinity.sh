#!/bin/bash

timeOut=0    
name_tag="1E2"   
data_file="/uufs/chpc.utah.edu/common/home/u1520754/corsika_inputs/sampled_muon_parameters_${name_tag}.txt"

# Check if data file exists
if [[ ! -f "$data_file" ]]; then
    echo "Error: Data file $data_file not found!"
    exit 1
fi

# Read data into array using temp file (most compatible)
temp_file=$(mktemp)
tail -n +2 "$data_file" > "$temp_file"

lines=()
while IFS= read -r line; do
    lines+=("$line")
done < "$temp_file"

rm "$temp_file"  # Clean up

# Fixed parameters
particles=(6)  # -mu
mass=(0.105658) #mass in GeV/c^2
telescope_radius=150.0  # Or read from a fourth column if you want to vary this
events_file=1   

# Calculate base coordinates (keep your existing logic)
base_shift_x=2705.36
base_shift_y=11522.87

work_direc="/uufs/chpc.utah.edu/common/home/u1520754"
home_direc="/uufs/chpc.utah.edu/common/home/u1520754"

echo "Processing ${#lines[@]} parameter combinations..."

# Loop through each line of sampled parameters
for i in "${!lines[@]}"; do
    line="${lines[$i]}"
    
    # Parse the three columns
    read -r energy production_height zenith_event <<< "$line"
    
    # Convert scientific notation if needed
    #energy=$(printf "%.2e" "$energy")
    #production_height=$(printf "%.2e" "$production_height")
    #zenith_event=$(printf "%.2e" "$zenith_event")
    echo "Processing sample before $((i+1)): E=${energy} GeV, h=${production_height} cm, zenith=${zenith_event}°"
    if (( $(echo "$zenith_event > 88" | bc -l) )); then
    zenith_event=88.0  # Also fix the assignment
    fi
    echo "Processing sample after $((i+1)): E=${energy} GeV, h=${production_height} cm, zenith=${zenith_event}°"
    particle=${particles[0]}  # Using single particle type
    
    # Calculate momentum
    P_long=$(awk -v E="$energy" -v m="$mass" 'BEGIN{
        if (E <= m) exit 2;
        printf("%.7E", sqrt(E*E - m*m));
    }') || { echo "E=$energy invalid (E<=m)"; continue; }
    
    # Calculate telescope positioning
    zenith_tele=$zenith_event
    theta_rad=$(echo "scale=10; $zenith_tele * 4*a(1) / 180" | bc -l)
    
    # Use a fixed telescope radius or make it variable
    RANDOM=$((1000+$i))
    # Calculate positions
    azimuth_min=250
    azimuth_max=300
    azimuth_tele=$(echo "scale=1; $azimuth_min + ($azimuth_max - $azimuth_min) * $RANDOM / 32767" | bc -l)
    azimuth_rad=$(echo "scale=10; $azimuth_tele * 4*a(1) / 180" | bc -l)

    x=$(echo "scale=6; $telescope_radius * s($theta_rad)/c($theta_rad) * c($azimuth_rad)" | bc -l)
    y=$(echo "scale=6; $telescope_radius * s($theta_rad)/c($theta_rad) * s($azimuth_rad)" | bc -l)
    x=$(echo "-1 * $x" | bc -l)
    y=$(echo "-1 * $y" | bc -l)
    
    x1=$(echo "$x - $base_shift_x" | bc -l)
    y1=$(echo "$y + $base_shift_y" | bc -l)
    

    # Generate unique identifiers
    evtnr=$((1000+$i))  # Unique event number
    zenith_min=$zenith_event
    zenith_max=$zenith_event
    azimuth_min=$azimuth_tele
    azimuth_max=$azimuth_tele
    # Create output directories
    out_direc=$work_direc/corsika_results/notch_r${telescope_radius}_Trinity_MC_${name_tag}/out_sib23d-pId${particle}
    mkdir -p "$out_direc"

    dat_direc="$out_direc/OUTPUT_${evtnr}"

    # Clean and create data directory
    if [ -d "$dat_direc" ]; then
        rm -rf "$dat_direc"
    fi
    mkdir -p "$dat_direc"
    
    # Create STACKIN file
    STACKIN_FILE="$dat_direc/STACKIN_$evtnr"
    cp "$home_direc/corsika_inputs/muon_stackin" "$STACKIN_FILE"
    
    sed -i "s#TTPRMPARTT#$particle#g" "$STACKIN_FILE"
    sed -i "s#E_tot#$energy#g" "$STACKIN_FILE"
    sed -i "s#P_long#$P_long#g" "$STACKIN_FILE"
    
    # Create INPUT file
    INPUT_FILE="$dat_direc/INPUTS_$evtnr"
    cp "$home_direc/corsika_inputs/input-sib23d-params_muon" "$INPUT_FILE"
    
    # Generate seeds based on sample index
    seed1="${evtnr}1"
    seed2="${evtnr}2"
    seed3="${evtnr}3"
    
    # Replace all placeholders
    sed -i "s#TTPRMPARTT#$particle#g" "$INPUT_FILE"
    sed -i "s#TTERANGEMINTT#$energy#g" "$INPUT_FILE"
    sed -i "s#TTERANGEMAXTT#$energy#g" "$INPUT_FILE"
    sed -i "s#TTTHETAPMINTT#$zenith_min#g" "$INPUT_FILE"
    sed -i "s#TTTHETAPMAXTT#$zenith_max#g" "$INPUT_FILE"
    sed -i "s#PHIPMIN#$azimuth_min#g" "$INPUT_FILE"
    sed -i "s#PHIPMAX#$azimuth_max#g" "$INPUT_FILE"
    sed -i "s#TTRUNNRTT#$evtnr#g" "$INPUT_FILE"
    sed -i "s#TTSEED1TT#$seed1#g" "$INPUT_FILE"
    sed -i "s#TTSEED2TT#$seed2#g" "$INPUT_FILE"
    sed -i "s#TTSEED3TT#$seed3#g" "$INPUT_FILE"
    sed -i "s#OUTPUTDIR#${dat_direc}/#g" "$INPUT_FILE"
    sed -i "s#OUTFILEDIR#${dat_direc}/telescope.dat#g" "$INPUT_FILE"
    sed -i "s#PROHEIGHT#$production_height#g" "$INPUT_FILE"
    sed -i "s#TELERADI#$telescope_radius#g" "$INPUT_FILE"
    sed -i "s#XSHIFT#$x#g" "$INPUT_FILE"
    sed -i "s#YSHIFT#$y#g" "$INPUT_FILE"
    sed -i "s#Xone#$x1#g" "$INPUT_FILE"
    sed -i "s#Yone#$y1#g" "$INPUT_FILE"
    sed -i "s#STACKINFILE#$STACKIN_FILE#g" "$INPUT_FILE"
    
    # Submit job
    submit="$home_direc/cluster/run_corsika_sibyll23d.slurm ${INPUT_FILE} ${dat_direc}/log_out_$evtnr"
    error_treatment="--output=${dat_direc}/log_$evtnr.out --error=${dat_direc}/log_$evtnr.err"
    
    sbatch $error_treatment $submit
    echo "Submitted job $((i+1))/${#lines[@]}: $evtnr"
    
    sleep "$timeOut"
done

echo "Submitted ${#lines[@]} jobs total!"
