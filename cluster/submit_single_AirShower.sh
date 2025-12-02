#!/bin/bash

timeOut=0       # Timeout to wait between submits 

i_values=($(seq 1 10)) # Job numbers
echo "Number of i_values: ${#i_values[@]}"
particles=(14)  # corresponding particles for each i 6 is -mu, 5 is +mu, 14 proton, 132 -tau
energies=(1.00E6)  # corresponding energies for each i
telescope_radius=(15000.0) # telescope radius in cm
events_file=1   # Number of events per file

# Define the RANGES for random generation
zenith_min=85    # minimum zenith angle in degrees
zenith_max=90.0   # maximum zenith angle in degrees
azimuth_min=275   # minimum azimuth angle in degrees
azimuth_max=285 # maximum azimuth angle in degrees

work_direc="/uufs/chpc.utah.edu/common/home/u1520754"
home_direc="/uufs/chpc.utah.edu/common/home/u1520754"

# Check if particles and energies have either 1 entry or match i_values in size
num_particles=${#particles[@]}
num_energies=${#energies[@]}
num_i_values=${#i_values[@]}

if [[ $num_particles -ne 1 && $num_particles -ne $num_i_values ]]; then
    echo "Error: particles must have either one entry or the same number of entries as i_values."
    exit 1
fi

if [[ $num_energies -ne 1 && $num_energies -ne $num_i_values ]]; then
    echo "Error: energies must have either one entry or the same number of entries as i_values."
    exit 1
fi

for index in "${!i_values[@]}"; do
  i=${i_values[$index]}  # get the current job number

  # Generate RANDOM zenith and azimuth for THIS event
  zenith_tele=$(awk -v min=$zenith_min -v max=$zenith_max 'BEGIN{srand(); print min+rand()*(max-min)}')
  azimuth_tele=$(awk -v min=$azimuth_min -v max=$azimuth_max 'BEGIN{srand(); print min+rand()*(max-min)}')
  
  echo "Event $i: Random Zenith = $zenith_tele deg, Azimuth = $azimuth_tele deg"

  # Convert to radians for calculation
  theta_rad=$(echo "scale=10; $zenith_tele * 4*a(1) / 180" | bc -l)
  azimuth_rad=$(echo "scale=10; $azimuth_tele * 4*a(1) / 180" | bc -l)

  # Calculate X and Y shift using correct bc syntax
  x=$(echo "scale=6; $telescope_radius * s($theta_rad)/c($theta_rad) * c(1 * $azimuth_rad)" | bc -l)
  y=$(echo "scale=6; $telescope_radius * s($theta_rad)/c($theta_rad) * s(1 * $azimuth_rad)" | bc -l)
  x=$(echo "-1 * $x" | bc -l)
  y=$(echo "-1 * $y" | bc -l)

  echo "Calculated X = $x cm, Y = $y cm"

  # If particles array has only one entry, use it for all; otherwise, use the indexed value
  if [[ $num_particles -eq 1 ]]; then
      particle=${particles[0]}
  else
      particle=${particles[$index]}
  fi

  # If energies array has only one entry, use it for all; otherwise, use the indexed value
  if [[ $num_energies -eq 1 ]]; then
      energy=${energies[0]}
  else
      energy=${energies[$index]}
  fi

  evtnr=$((1000+$i))
  out_direc=$work_direc/corsika_results/notch_r${telescope_radius}_Trinity_MC_proton/out_sib23d-pId${particle}
  if [[ ! -d $out_direc ]]; then
    mkdir -p "$out_direc"
  fi
  dat_direc="$out_direc/OUTPUT_${evtnr}"

  # Check if the directory exists
  if [ -d "$dat_direc" ]; then
      echo "The directory $dat_direc exists. Deleting it..."
      rm -rf "$dat_direc"  # Remove the directory and all its contents
  fi

  # Create the directory again
  mkdir -p "$dat_direc"
  echo "The directory $dat_direc has been created."

  INPUT_FILE="$dat_direc/INPUTS_$evtnr"
  cp "$home_direc/Muon_Trinity/corsika_inputs/input-sib23d-params_proton" "$INPUT_FILE"

  day=$(date +"%d")
  seed1="${i}10"
  seed2="${i}20"
  seed3="${i}30"
  seed4="${i}40"

  # Safely replace placeholders - use the SAME value for min and max to fix the angle
  sed -i "s#TTPRMPARTT#$particle#g" "$INPUT_FILE"
  sed -i "s#TTERANGEMINTT#$energy#g" "$INPUT_FILE"
  sed -i "s#TTERANGEMAXTT#$energy#g" "$INPUT_FILE"
  sed -i "s#TTTHETAPMINTT#$zenith_tele#g" "$INPUT_FILE"
  sed -i "s#TTTHETAPMAXTT#$zenith_tele#g" "$INPUT_FILE"
  sed -i "s#PHIPMIN#$azimuth_tele#g" "$INPUT_FILE"
  sed -i "s#PHIPMAX#$azimuth_tele#g" "$INPUT_FILE"
  sed -i "s#TTRUNNRTT#$evtnr#g" "$INPUT_FILE"
  sed -i "s#TTSEED1TT#$seed1#g" "$INPUT_FILE"
  sed -i "s#TTSEED2TT#$seed2#g" "$INPUT_FILE"
  sed -i "s#TTSEED3TT#$seed3#g" "$INPUT_FILE"
  sed -i "s#TTSEED4TT#$seed4#g" "$INPUT_FILE"
  sed -i "s#OUTPUTDIR#${dat_direc}/#g" "$INPUT_FILE"
  sed -i "s#OUTFILEDIR#${dat_direc}/telescope.dat#g" "$INPUT_FILE"
  sed -i "s#TELERADI#$telescope_radius#g" "$INPUT_FILE"
  sed -i "s#XSHIFT#$x#g" "$INPUT_FILE"
  sed -i "s#YSHIFT#$y#g" "$INPUT_FILE"

  submit="$home_direc/Muon_Trinity/cluster/run_corsika_sibyll23d.slurm ${INPUT_FILE} ${dat_direc}/log_out_$evtnr"

  error_treatment="--output=${dat_direc}/log_$evtnr.out --error=${dat_direc}/log_$evtnr.err"

  sbatch $error_treatment $submit

  echo "sbatch $error_treatment $submit"

  sleep "$timeOut"
done
