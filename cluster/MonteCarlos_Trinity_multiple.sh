#!/bin/bash

timeOut=0    
name_tag="1E3_1E7"   
run_tag="0"
sample_size=5000
event_ids=($(seq 0 4999))
events_per_job=20
data_file="/uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/corsika_inputs/sampled_muon_parameters_${name_tag}.txt"
echo "Using data file: $data_file sampling $sample_size events for tag $run_tag"

# Check if data file exists
if [[ ! -f "$data_file" ]]; then
    echo "Error: Data file $data_file not found!"
    exit 1
fi

# Define the sampled data file path
sampled_data_file="/uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/corsika_inputs/sampled_data_${name_tag}_${sample_size}events.txt"
# Check if sampled file already exists
if [[ -f "$sampled_data_file" ]]; then
    echo "Sample file already exists: $sampled_data_file"
    echo "Skipping file creation and using existing sample."
    
    # Verify the existing file has the expected number of lines (header + sample_size)
    existing_lines=$(($(wc -l < "$sampled_data_file") - 1))  # Subtract header
    echo "Existing file has $existing_lines data lines (expected: $sample_size)"
    
    if [[ $existing_lines -ne $sample_size ]]; then
        echo "WARNING: Line count mismatch! Expected $sample_size, found $existing_lines"
        echo "You may want to delete the file and regenerate: rm '$sampled_data_file'"
    fi

else
    echo "Sample file not found. Creating new sample..."
    
    # Create temporary file with better error handling
    temp_sample=$(mktemp) || { echo "Failed to create temp file"; exit 1; }
    echo "Debug: Using temp file: $temp_sample"

    # Get random sample of lines (skip header)
    echo "Sampling $sample_size lines from $data_file..."
    tail -n +2 "$data_file" | shuf -n "$sample_size" > "$temp_sample"

    # Check if sampling worked
    sample_lines=$(wc -l < "$temp_sample")
    echo "Debug: Sampled $sample_lines lines to temp file"

    if [[ $sample_lines -eq 0 ]]; then
        echo "ERROR: No lines were sampled!"
        echo "Original file info:"
        wc -l "$data_file"
        head -3 "$data_file"
        rm "$temp_sample"
        exit 1
    fi

    echo "Saving sampled data to: $sampled_data_file"
    
    # Create directory if needed
    mkdir -p "$(dirname "$sampled_data_file")"

    # Add header to the saved file
    head -n 1 "$data_file" > "$sampled_data_file"

    # Append the sampled data
    cat "$temp_sample" >> "$sampled_data_file"
    
    # Clean up temp file
    rm "$temp_sample"
    
    echo "Successfully created sample file with $sample_lines events"
fi

# Now load the sample data (either existing or newly created)
temp_file=$(mktemp)
tail -n +2 "$sampled_data_file" > "$temp_file"

lines=()
while IFS= read -r line; do
    [[ -n "$line" ]] && lines+=("$line")
done < "$temp_file"

rm "$temp_file"  # Clean up

# Safety check
if [[ ${#lines[@]} -eq 0 ]]; then
    echo "ERROR: No lines loaded into array!"
    exit 1
fi

# Pre-calculate fixed parameters once
particles=5 # + muon (5), - muon (6)
mass=0.105658
telescope_radius=150.0
work_direc="/uufs/chpc.utah.edu/common/home/u1520754"
home_direc="/uufs/chpc.utah.edu/common/home/u1520754"
azimuth_min=275
azimuth_max=285

# Pre-calculate constants
pi4=$(echo "4*a(1)" | bc -l)
deg_to_rad=$(echo "scale=15; $pi4/180" | bc -l)

# Batch directory creation
out_base="$work_direc/corsika_results/notch_r${telescope_radius}_Trinity_Demo_${name_tag}_${run_tag}/out_sib23d-pId${particles}"
mkdir -p "$out_base"

# Pre-read template files once
stackin_template=$(cat "$home_direc/corsika_inputs/muon_stackin")
input_template=$(cat "$home_direc/corsika_inputs/input-sib23d-params_muon")

# Define specific event IDs to process

total_events=${#event_ids[@]}
total_jobs=$(((total_events + events_per_job - 1) / events_per_job))

echo "Submitting $total_jobs jobs with $events_per_job events each..."

# Create sequential job runner script
runner_script="$work_direc/sequential_runner.sh"
cat > "$runner_script" << 'EOF'
#!/bin/bash

# Arguments: start_idx end_idx job_id
start_idx=$1
end_idx=$2
job_id=$3

# Source all the variables and data
EOF

# Append all variables to the runner script
cat >> "$runner_script" << EOF
particles=$particles
mass=$mass
telescope_radius=$telescope_radius
work_direc="$work_direc"
home_direc="$home_direc"
azimuth_min=$azimuth_min
azimuth_max=$azimuth_max
deg_to_rad="$deg_to_rad"
out_base="$out_base"
name_tag="$name_tag"
run_tag="$run_tag"

# Recreate the lines array
lines=(
EOF

# Add all lines to the script
for line in "${lines[@]}"; do
    echo "\"$line\"" >> "$runner_script"
done

# Close the lines array first (keep single-quoted so nothing gets expanded)
cat >> "$runner_script" << 'EOF'
)
EOF

# Inject event_ids from the parent script (needs expansion now)
cat >> "$runner_script" << EOF
# Recreate event_ids array from parent
event_ids=(${event_ids[@]})
EOF

# Resume single-quoted heredoc for the rest of the child script
cat >> "$runner_script" << 'EOF'
# Read templates
stackin_template=$(cat "$home_direc/corsika_inputs/muon_stackin")
input_template=$(cat "$home_direc/corsika_inputs/input-sib23d-params_muon")

echo "Job $job_id: Processing events $start_idx to $end_idx sequentially..."

# Process events sequentially
for ((j=start_idx; j<=end_idx; j++)); do
    i=${event_ids[$j]}

    # Get line index
    line="${lines[$i]}"

    # Parse the 4 columns
    read -r energy production_height zenith_event weight <<< "$line"

    # Zenith angle check
    zenith_event=$(awk -v z="$zenith_event" 'BEGIN{print (z > 89.9) ? 89.9 : z}')

    evtnr=$((1000+$i))

    echo "Job $job_id: Processing event $((j-start_idx+1))/$((end_idx-start_idx+1)): evtnr=$evtnr"
    
    # Calculate momentum
    P_long=$(awk -v E="$energy" -v m="$mass" 'BEGIN{
        if (E <= m) exit 2;
        printf("%.7E", sqrt(E*E - m*m));
    }') || { echo "E=$energy invalid (E<=m)"; continue; }

    # All trigonometric calculations
    read -r theta_rad azimuth_tele azimuth_rad x y <<< $(awk -v zenith="$zenith_event" -v amin="$azimuth_min" -v amax="$azimuth_max" -v seed="$evtnr" -v tr="$telescope_radius" -v deg2rad="$deg_to_rad" 'BEGIN{
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

    # Final coordinate calculations

    # Create directories
    dat_direc="$out_base/OUTPUT_${evtnr}"
    [ -d "$dat_direc" ] && rm -rf "$dat_direc"
    mkdir -p "$dat_direc"

    # Generate files
    STACKIN_FILE="$out_base/STACKIN_$evtnr"
    INPUT_FILE="$dat_direc/INPUTS_$evtnr"
    # Create STACKIN file
    echo "$stackin_template" | awk -v particle="$particles" -v energy="$energy" -v plong="$P_long" '{
        gsub(/TTPRMPARTT/, particle);
        gsub(/E_tot/, energy);
        gsub(/P_long/, plong);
        print;
    }' > "$STACKIN_FILE"

    # Create INPUT file
    seed1="${evtnr}1${run_tag}"
    seed2="${evtnr}2${run_tag}"
    seed3="${evtnr}3${run_tag}"
    seed4="${evtnr}4${run_tag}"

    echo "$input_template" | awk -v particle="$particles" -v energy="$energy" -v zenith="$zenith_event" -v azimin="$azimuth_tele" -v azimax="$azimuth_tele" -v evtnr="$evtnr" -v seed1="$seed1" -v seed2="$seed2" -v seed3="$seed3" -v seed4="$seed4" -v datdir="$dat_direc" -v proheight="$production_height" -v teleradius="$telescope_radius" -v xshift="$x" -v yshift="$y" -v stackfile="$STACKIN_FILE" '{
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
    }' > "$INPUT_FILE"

    # Run corsika and wait for completion
    module load gcc
    echo "Input file: $INPUT_FILE"
    echo "Output file: ${dat_direc}/log_out_$evtnr"
    cd /uufs/chpc.utah.edu/common/home/u1520754/corsika-78010/run/
    ./corsika78010Linux_SIBYLL_urqmd < "$INPUT_FILE" > "${dat_direc}/log_out_$evtnr"

    echo "Job $job_id: Completed event $((j-start_idx+1))/$((end_idx-start_idx+1)): evtnr=$evtnr"
done

echo "Job $job_id: All events completed!"
EOF

chmod +x "$runner_script"
##srun "$home_direc/cluster/run_corsika_sibyll23d.slurm" "$INPUT_FILE" "${dat_direc}/log_out_$evtnr"
# Submit jobs
for ((job=0; job<total_jobs; job++)); do
    start_idx=$((job * events_per_job))
    end_idx=$((start_idx + events_per_job - 1))
    
    # Adjust last job if necessary
    if [ $end_idx -ge $total_events ]; then
        end_idx=$((total_events - 1))
    fi
    
    log_dir="$out_base/job_logs"
    mkdir -p "$log_dir"
    
    sbatch --account=owner-guest \
       --partition=kingspeak-shared-guest \
       --time=40:00:00 \
       --mem=16G \
       --output="${log_dir}/job_${job}.out" \
       --error="${log_dir}/job_${job}.err" \
       --job-name="corsika_job_${job}" \
       "$runner_script" $start_idx $end_idx $job
    
    echo "Submitted job $((job+1))/$total_jobs: events $start_idx-$end_idx"
done

echo "Submitted $total_jobs jobs total!"
