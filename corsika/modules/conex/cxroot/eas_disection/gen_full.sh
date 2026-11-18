#!/bin/bash

echo "This runs conex to produce an output file containing the highest energy interactions"
echo

# rm -f *.log *.root *~

export CONEX_ROOT=${PWD}/param_highE
export seed=90
export energy=17

../bin/conex2r -x full -n 50 -s ${seed} -e ${energy} -E ${energy} -p 100 -m 5 -L 1 -F full-interactions.${seed}.root >& full.${seed}.log

