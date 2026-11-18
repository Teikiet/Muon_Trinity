#!/bin/bash

# (c) Copyright 2018 CORSIKA Project, corsika8@kit.edu
#
# See file AUTHORS for a list of contributors.
#
# This software is distributed under the terms of the 3-clause BSD license.
# See file LICENSE for a full version of the license.

# with this script you can plot an animation of output of TrackWriter

tracks_dat="$1"
if [ -z "$tracks_dat" ]; then
  echo "usage: $0 <tracks.dat> [output.gif]" >&2
  exit 1
fi

directory=$(dirname "$tracks_dat")
hadrons_dat="$directory/hadrons.dat"
muons_dat="$directory/muons.dat"

if [ ! -e "$tracks_dat" ]; then
  echo "$tracks_dat does not exist" >&2
  exit 2
fi

if [ "$muons_dat" -ot "$tracks_dat" ] || [ "$hadrons_dat" -ot "$tracks_dat" ]; then
  echo "splitting $tracks_dat into $muons_dat and $hadrons_dat."
  cat "$tracks_dat" | egrep '^\s+-*13\s' > "$muons_dat"
  cat "$tracks_dat" | egrep -v '^\s+-*13\s' > "$hadrons_dat"
fi

output="$2"
if [ -z "$output" ]; then
  output="$tracks_dat.gif"
fi

echo "creating $output..."

cat <<EOF | gnuplot
set term gif animate size 900,900
set output "$output"

#set zrange [0:40e3]
#set xrange [-10:10]
#set yrange [-10:10]
set xlabel "x / m"
set ylabel "y / m"
set zlabel "z / m"
set title "CORSIKA 8 preliminary"

do for [t=0:359:1] {
# for separate file per angle:
#	set output sprintf("%03d_$output", t)

	set view 80, t
        splot "$muons_dat" u 3:4:5:6:7:8 w vectors nohead lt rgb "red" t "", "$hadrons_dat" u 3:4:5:6:7:8 w vectors nohead  lc rgb "black" t ""
}
EOF

exit $?
