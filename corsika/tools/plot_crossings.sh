#!/bin/sh

# (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
#
# See file AUTHORS for a list of contributors.
#
# This software is distributed under the terms of the 3-clause BSD license.
# See file LICENSE for a full version of the license.

# with this script you can plot an animation of output of TrackWriter

crossings_dat=$1
if [ -z "$crossings_dat" ]; then
  echo "usage: $0 <crossings.dat> [output.gif]" >&2
  exit 1
fi

output=$2
if [ -z "$output" ]; then
  output="$crossings_dat.gif"
fi

cat <<EOF | gnuplot
set term gif animate size 600,600
set output "$output"

set xlabel "x / m"
set ylabel "y / m"
set zlabel "z / m"
set title "CORSIKA 8 preliminary"

do for [t=0:360:1] {
	set view 60, t
        splot "$crossings_dat" u 2:3:4 w points t ""
}
EOF

exit $?
