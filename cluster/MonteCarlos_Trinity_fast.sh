#!/bin/bash

timeOut=0    
name_tag="1E3_1E7"   
run_tag="6"
data_file="/uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/corsika_inputs/sampled_data_${name_tag}_5000events.txt"
# Check if data file exists
if [[ ! -f "$data_file" ]]; then
    echo "Error: Data file $data_file not found!"
    exit 1
fi
echo "Running tag: $name_tag with run tag: $run_tag"
# Read data into array using temp file (most compatible)
temp_file=$(mktemp)
tail -n +2 "$data_file" > "$temp_file"

lines=()
while IFS= read -r line; do
    [[ -n "$line" ]] && lines+=("$line")
done < "$temp_file"

rm "$temp_file"  # Clean up

# Pre-calculate fixed parameters once
particles=6
mass=0.105658
telescope_radius=1500.0
work_direc="/uufs/chpc.utah.edu/common/home/u1520754"
home_direc="/uufs/chpc.utah.edu/common/home/u1520754"
azimuth_min=275
azimuth_max=285

# Pre-calculate constants
pi4=$(echo "4*a(1)" | bc -l)  # Calculate π once
deg_to_rad=$(echo "scale=15; $pi4/180" | bc -l)  # π/180 conversion factor

echo "Processing ${#lines[@]} parameter combinations..."

# Batch directory creation - create all main directories first
out_base="$work_direc/corsika_results/notch_r${telescope_radius}_Trinity_Demo_${name_tag}_${run_tag}/out_sib23d-pId${particles}"
mkdir -p "$out_base"

# Pre-read template files once
stackin_template=$(cat "$home_direc/corsika_inputs/muon_stackin")
input_template=$(cat "$home_direc/corsika_inputs/input-sib23d-params_muon")

# Define specific event IDs to process (0-based indexing)
event_ids=(4997 4947 3595 3660 3676 671 4958 2233 2916 3277 3515 3658 3669 4301 4945 2235 2940 4110 4243 4246 4351 4354 3311 3595 3626 3654 3664 4862 4918)
echo "Processing ${#event_ids[@]} specific events..."

# Process in chunks to avoid overwhelming the scheduler  
chunk_size=50
total_events=${#event_ids[@]}

for ((chunk_start=0; chunk_start<total_events; chunk_start+=chunk_size)); do
    chunk_end=$((chunk_start + chunk_size - 1))
    if [ $chunk_end -ge $total_events ]; then
        chunk_end=$((total_events - 1))
    fi
    
    echo "Processing chunk: events $((chunk_start+1))-$((chunk_end+1))"
    
    # Process chunk
    for ((j=chunk_start; j<=chunk_end; j++)); do
        i=${event_ids[$j]}  # Get the actual line index
        
        # Validate index
        if [ $i -ge ${#lines[@]} ]; then
            echo "Warning: Event ID $i exceeds available lines (${#lines[@]}), skipping..."
            continue
        fi
        line="${lines[$i]}"
        
        # Parse the 4 columns
        read -r energy production_height zenith_event weight <<< "$line"

        # Zenith angle check - simplified
        zenith_event=$(awk -v z="$zenith_event" 'BEGIN{print (z > 89.9) ? 89.9 : z}')
        
        evtnr=$((1000+$i))
        
        # Calculate momentum - simplified
        P_long=$(awk -v E="$energy" -v m="$mass" 'BEGIN{
            if (E <= m) exit 2;
            printf("%.7E", sqrt(E*E - m*m));
        }') || { echo "E=$energy invalid (E<=m)"; continue; }
        
        # All trigonometric calculations in one awk call
        read -r theta_rad azimuth_tele azimuth_rad x y <<< $(awk -v run_tag="$run_tag" -v zenith="$zenith_event" -v amin="$azimuth_min" -v amax="$azimuth_max" -v seed="$((1000+i))" -v tr="$telescope_radius" -v deg2rad="$deg_to_rad" 'BEGIN{
            seed0 = seed "0" run_tag;
            srand(seed0);
            theta_rad = zenith * deg2rad;
            azimuth_tele = amin + (amax - amin) * rand();
            azimuth_rad = azimuth_tele * deg2rad;
            
            tan_theta = sin(theta_rad)/cos(theta_rad);
            x = -tr * tan_theta * cos(azimuth_rad);
            y = -tr * tan_theta * sin(azimuth_rad);
            
            printf("%.10f %.1f %.10f %.6f %.6f", theta_rad, azimuth_tele, azimuth_rad, x, y);
        }')
        
        
        # Create directories
        dat_direc="$out_base/OUTPUT_${evtnr}"
        [ -d "$dat_direc" ] && rm -rf "$dat_direc"
        if [ -d "$dat_direc" ]; then
            rm -f "$dat_direc/DAT${evtnr}"*
            rm -f "$dat_direc/telescope.dat"
            rm -f "$$out_base/STACKIN_${evtnr}"
            rm -f "$dat_direc/INPUTS_${evtnr}"
        fi
        mkdir -p "$dat_direc"
        
        # Generate all files using string operations instead of multiple sed calls
        STACKIN_FILE="$out_base/STACKIN_$evtnr"
        INPUT_FILE="$dat_direc/INPUTS_$evtnr"
        
        # Create STACKIN file with single substitution
        {
            echo "$stackin_template" | awk -v particle="$particles" -v energy="$energy" -v plong="$P_long" '
            {
                gsub(/TTPRMPARTT/, particle);
                gsub(/E_tot/, energy);
                gsub(/P_long/, plong);
                print;
            }'
        } > "$STACKIN_FILE"
        
        # Create INPUT file with single substitution
        seed1="${evtnr}1${run_tag}"
        seed2="${evtnr}2${run_tag}"
        seed3="${evtnr}3${run_tag}"
        seed4="${evtnr}4${run_tag}"
        {
            echo "$input_template" | awk -v particle="$particles" -v energy="$energy" -v zenith="$zenith_event" -v azimin="$azimuth_tele" -v azimax="$azimuth_tele" -v evtnr="$evtnr" -v seed1="$seed1" -v seed2="$seed2" -v seed3="$seed3" -v seed4="$seed4" -v datdir="$dat_direc" -v proheight="$production_height" -v teleradius="$telescope_radius" -v xshift="$x" -v yshift="$y" -v stackfile="$STACKIN_FILE" '
            {
                gsub(/TTPRMPARTT/, particle);
                gsub(/TTERANGEMINTT/, energy);
                gsub(/TTERANGEMAXTT/, energy);
                gsub(/TTTHETAPMINTT/, zenith);
                gsub(/TTTHETAPMAXTT/, zenith);
                gsub(/PHIPMIN/, azimin);
                gsub(/PHIPMAX/, azimax);
                gsub(/TTRUNNRTT/, evtnr);
                gsub(/TTSEED1TT/, seed1);
                gsub(/TTSEED2TT/, seed2);
                gsub(/TTSEED3TT/, seed3);
                gsub(/TTSEED4TT/, seed4);
                gsub(/OUTPUTDIR/, datdir "/");
                gsub(/OUTFILEDIR/, datdir "/telescope.dat");
                gsub(/PROHEIGHT/, proheight);
                gsub(/TELERADI/, teleradius);
                gsub(/XSHIFT/, xshift);
                gsub(/YSHIFT/, yshift);
                gsub(/STACKINFILE/, stackfile);
                print;
            }'
        } > "$INPUT_FILE"
        
        # Submit job
        sbatch --output="${dat_direc}/log_${evtnr}.out" --error="${dat_direc}/log_${evtnr}.err" "$home_direc/cluster/run_corsika_sibyll23d.slurm" "$INPUT_FILE" "${dat_direc}/log_out_$evtnr"
        
        echo "Submitted job $((j+1))/${#event_ids[@]}: evtnr=$evtnr (line=$i)"
    done
    
    # Small delay between chunks to avoid overwhelming scheduler
    if [ $((chunk_end+1)) -lt $total_events ]; then
        echo "Pausing 0 second between chunks..."
        sleep 0
    fi
done

echo "Submitted ${#event_ids[@]} jobs total!"
