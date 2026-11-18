#!/bin/bash

echo "This runs conex to re-process the highest energy interactions of one air-shower individually"
echo

export CONEX_ROOT=${PWD}/param_fullMC
export seed=70
export energy=19

export energy=17
export i=1
../bin/conex2r -x fullHad -n 10 -e ${energy} -E ${energy} -p 100 -m 5 -L -${i} -F fullHad-interactions.${energy}.root >& fullHad.${energy}.log


for i in `seq 1 100`; do

    prefix="int_${i}"
#    ../bin/conex2r -x ${prefix} -n 1 -s ${seed} -e ${energy} -E ${energy} -p 100 -m 5 -L -${i} -F full-interactions.${seed}.root >& ${prefix}.${seed}.log

done

