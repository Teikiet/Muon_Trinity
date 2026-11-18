#!/bin/sh

# (c) Copyright 2023 CORSIKA Project, corsika8@kit.edu
#
# See file AUTHORS for a list of contributors.
#
# This software is distributed under the terms of the 3-clause BSD license.
# See file LICENSE for a full version of the license.

# This script strips off flrndm() and flrnlp() from libflukahp.a so that
# we can provide our own implementation.

flukalibOrig=`realpath $1`
target="$2"

if [ -d "$flukalibOrig" ]; then
    echo "\"$flukalibOrig\" is a directory. Please specify full path to libflukahp.a?" 1>&2
    exit 1
fi

if [ ! -f "$flukalibOrig" ]; then
    echo "\"$flukalibOrig\" is not a regular file. Did you specify the full path to libflukahp.a?" 1>&2
    exit 1
fi

if [ ! -r "$flukalibOrig" ]; then
    echo "\"$flukalibOrig\" not readable" 1>&2
    exit 1
fi

# check if FLUKA has the required symbols. If not, it is an imcompatible
# version (e.g. CERN FLUKA or too old)
if ! ar t "$flukalibOrig" | grep ndmhep.o >/dev/null; then
    echo "The provided libflukahp.a is incompatible." 1>&2
    exit 1
fi

cp "${flukalibOrig}" "${target}/libflukahp-norndm.a" && \
ar -d "${target}/libflukahp-norndm.a" flrndm.o flrnlp.o
