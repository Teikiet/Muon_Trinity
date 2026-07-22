#!/bin/bash

timeOut=2       # Timeout to wait between submits 
#FIXHEI  PROHEIGHT  0              Fixes the height (in cm) of the first interaction of hadronic primaries
#THETAP  TTTHETAPMINTT  TTTHETAPMAXTT                range zenith angle (deg)
#FIXHEI  PROHEIGHT  0              Fixes the height (in cm) of the first interaction of hadronic primaries
i_values=(1 2 3 4 5) # Job numbers 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
particles=(6)  # corresponding particles for each i 6 is -mu, 5 is +mu, 14 proton, 132 -tau
energies=(1.00E6)  # corresponding energies for each i
production_height=(20.E5) # production height in cm
telescope_radius=(1500.) # telescope radius in cm
events_file=1   # Number of events per file

zenith_tele=88.
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
#Shift the coordinates to the origin
# power-law parameters from your fit
Ax=1726813325.693995
Bx=-0.961981573607791
Ay=7402197702.982274
By=-0.9637319068917075
# your primary energy in GeV (or whatever units you used in the fit)
energy=${energies[0]}     # <-- you already have this as “energies=(1.00E3)”, extract it if needed
base_shift_x=$(awk -v A=$Ax -v B=$Bx -v E=$energy \
             'BEGIN{ printf "%.10f", A * (E**B) }')
base_shift_y=$(awk -v A=$Ay -v B=$By -v E=$energy \
             'BEGIN{ printf "%.10f", A * (E**B) }')
echo "Base shift X = $base_shift_x cm"
echo "Base shift Y = $base_shift_y cm"
x1=$(echo "$x - $base_shift_x" | bc -l) 
y1=$(echo "$y + $base_shift_y" | bc -l)
echo "Calculated X = $x cm, Y = $y cm"
echo "Calculated X1 = $x1 cm, Y1 = $y1 cm"