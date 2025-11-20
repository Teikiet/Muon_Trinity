#!/bin/bash

timeOut=2       # Timeout to wait between submits 
#FIXHEI  PROHEIGHT  0              Fixes the height (in cm) of the first interaction of hadronic primaries
#THETAP  TTTHETAPMINTT  TTTHETAPMAXTT                range zenith angle (deg)
#FIXHEI  PROHEIGHT  0              Fixes the height (in cm) of the first interaction of hadronic primaries
i_values=(1) #($(seq 1 10)) # Job numbers 1 2 3 4 5 6 7 8 9 10
particles=(6)  # corresponding particles for each i 6 is -mu, 5 is +mu, 14 proton, 132 -tau
mass=(0.105658) #mass in GeV/c^2
energies=(1.00E6)  # corresponding energies for each i
production_height=(5.E5) # production height in cm
zenith_tele=88.
telescope_radius=(1500.) # telescope radius in cm
events_file=1   # Number of events per file

for E in "${energies[@]}"; do
  P_long=$(awk -v E="$E" -v m="$mass" 'BEGIN{
    if (E <= m) exit 2;
    printf("%.7E", sqrt(E*E - m*m));
  }') || { echo "E=$E invalid (E<=m)"; continue; }

  # Save/print as requested:
  echo "P_long=$P_long"
done

azimuth_tele=280.
theta_rad=$(echo "scale=10; $zenith_tele * 4*a(1) / 180" | bc -l)
azimuth_rad=$(echo "scale=10; $azimuth_tele * 4*a(1) / 180" | bc -l)
zenith_min=$zenith_tele
zenith_max=$zenith_tele  # zenith angle range in degrees
azimuth_min=$azimuth_tele  # azimuth angle in degrees
azimuth_max=$azimuth_tele  # azimuth angle in degrees
# Telescope radius
#telescope_radius=1500.0

# Calculate X and Y shift using correct bc syntax
x=$(echo "scale=6; $telescope_radius * s($theta_rad)/c($theta_rad) * c(1 * $azimuth_rad)" | bc -l)
y=$(echo "scale=6; $telescope_radius * s($theta_rad)/c($theta_rad) * s(1 * $azimuth_rad)" | bc -l)
x=$(echo "-1 * $x" | bc -l)
y=$(echo "-1 * $y" | bc -l)
#Shift the coordinates to the origin: 6: 2705.36,11522.87 ; 5: 26958.12,114816.8 ; 4: 265352.31,1131024 ; 3: 2211070,9428740
base_shift_x=2705.36
base_shift_y=11522.87
echo "Base shift X = $base_shift_x cm"
echo "Base shift Y = $base_shift_y cm"
x1=$(echo "$x - $base_shift_x" | bc -l) 
y1=$(echo "$y + $base_shift_y" | bc -l)
echo "Calculated X = $x cm, Y = $y cm"
echo "Calculated X1 = $x1 cm, Y1 = $y1 cm"

# Convert to radians using bc
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

  evtnr=$((100+$i))
  out_direc=$work_direc/corsika_results/notch_h${production_height}_r${telescope_radius}_Trinity/out_sib23d-pId${particle}-En${energy}-$zenith_tele
  if [[ ! -d $out_direc ]]; then
    mkdir -p "$out_direc"
  fi
  dat_direc="$out_direc/OUTPUT_${energy}_${particle}_${zenith_tele}_${evtnr}"

  # Check if the directory exists
  if [ -d "$dat_direc" ]; then
      echo "The directory $dat_direc exists. Deleting it..."
      rm -rf "$dat_direc"  # Remove the directory and all its contents
  fi

  # Create the directory again
  mkdir -p "$dat_direc"

  echo "The directory $dat_direc has been created."
  STACKIN_FILE="$out_direc/STACKIN_$evtnr"
  cp "$home_direc/corsika_inputs/muon_stackin" "$STACKIN_FILE"

  sed -i "s#TTPRMPARTT#$particle#g" "$STACKIN_FILE"
  sed -i "s#E_tot#$energy#g" "$STACKIN_FILE"
  sed -i "s#P_long#$P_long#g" "$STACKIN_FILE"

  INPUT_FILE="$dat_direc/INPUTS_$evtnr"
  cp "$home_direc/corsika_inputs/input-sib23d-params_muon" "$INPUT_FILE" #input-sib23d-params_muon

  day=$(date +"%d")
  seed1="${i}10"
  seed2="${i}20"
  seed3="${i}30"

  # Safely replace placeholders
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

  submit="$home_direc/cluster/run_corsika_sibyll23d.slurm ${INPUT_FILE} ${dat_direc}/log_out_$evtnr"

  error_treatment="--output=${dat_direc}/log_$evtnr.out --error=${dat_direc}/log_$evtnr.err"

  sbatch $error_treatment $submit

  echo "sbatch $error_treatment $submit"

  sleep "$timeOut"
done
