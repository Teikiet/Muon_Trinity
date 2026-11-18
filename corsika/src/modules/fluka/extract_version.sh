#!/bin/bash

# (c) Copyright 2024 CORSIKA Project, corsika8@kit.edu
#
# See file AUTHORS for a list of contributors.
#
# This software is distributed under the terms of the 3-clause BSD license.
# See file LICENSE for a full version of the license.

# This script extracts the FLUKA version from libflukahp.a.
# This is just a hack.

flukalib=`realpath $1`

if [ -d "$flukalib" ]; then
    echo "\"$flukalib\" is a directory. Please specify full path to libflukahp.a?" 1>&2
    exit 1
fi

if [ ! -f "$flukalib" ]; then
    echo "\"$flukalib\" is not a regular file. Did you specify the full path to libflukahp.a?" 1>&2
    exit 1
fi

if [ ! -r "$flukalib" ]; then
    echo "\"$flukalib\" not readable" 1>&2
    exit 1
fi

# check if FLUKA has the required symbols. If not, it is an imcompatible
# version (e.g. CERN FLUKA or too old)
if ! ar t "$flukalib" | grep ndmhep.o >/dev/null; then
    echo "The provided libflukahp.a is incompatible." 1>&2
    exit 1
fi

MYTMPDIR="$(mktemp -d)"
trap 'rm -rf -- "$MYTMPDIR"' EXIT # delete temp dir on exit
cd $MYTMPDIR

ar x ${flukalib} flukam.o || exit 2
objcopy -O binary --only-section=.rodata.str1.8 flukam.o /dev/stdout | perl -nle 'm/(FLUKA20\d\d Version \d.\d)/; print $1' || exit 3
