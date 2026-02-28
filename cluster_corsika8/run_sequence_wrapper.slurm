#!/bin/bash
# =============================================================================
# run_sequence_wrapper.slurm
#
# Reads a batch file and runs each simulation sequentially using the
# worker script (run_corsika8_trinity_chain.slurm) as a bash function.
# =============================================================================

#SBATCH --account=owner-guest
#SBATCH --partition=kingspeak-guest
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G

set -uo pipefail
# NOTE: not using -e so one failed sim doesn't kill the whole sequence

BATCH_FILE="$1"
WORKER_SCRIPT="$2"

if [ ! -f "${BATCH_FILE}" ]; then
    echo "ERROR: Batch file not found: ${BATCH_FILE}"
    exit 1
fi

if [ ! -f "${WORKER_SCRIPT}" ]; then
    echo "ERROR: Worker script not found: ${WORKER_SCRIPT}"
    exit 1
fi

TOTAL_SIMS=$(wc -l < "${BATCH_FILE}")
echo "============================================================"
echo "Sequence wrapper starting"
echo "  Batch file   : ${BATCH_FILE}"
echo "  Worker script: ${WORKER_SCRIPT}"
echo "  Simulations  : ${TOTAL_SIMS}"
echo "  SLURM Job ID : ${SLURM_JOB_ID}"
echo "  Node         : $(hostname)"
echo "  Start time   : $(date)"
echo "============================================================"

SIM_NUM=0
PASS=0
FAIL=0

while IFS='|' read -r ARG1 ARG2 ARG3 ARG4 ARG5 ARG6 ARG7 ARG8 ARG9 ARG10 ARG11 ARG12 ARG13 ARG14 ARG15; do

    SIM_NUM=$((SIM_NUM + 1))
    echo ""
    echo "============================================================"
    echo "Starting simulation ${SIM_NUM}/${TOTAL_SIMS} at $(date)"
    echo "============================================================"

    # Run the worker script as a bash subprocess
    # Pass the same arguments the slurm script expects
    bash "${WORKER_SCRIPT}" \
        "${ARG1}" "${ARG2}" "${ARG3}" "${ARG4}" "${ARG5}" \
        "${ARG6}" "${ARG7}" "${ARG8}" "${ARG9}" "${ARG10}" \
        "${ARG11}" "${ARG12}" "${ARG13}" "${ARG14}" "${ARG15}"

    EXIT_CODE=$?

    if [ ${EXIT_CODE} -eq 0 ]; then
        PASS=$((PASS + 1))
        echo "Simulation ${SIM_NUM}/${TOTAL_SIMS} PASSED (exit ${EXIT_CODE})"
    else
        FAIL=$((FAIL + 1))
        echo "Simulation ${SIM_NUM}/${TOTAL_SIMS} FAILED (exit ${EXIT_CODE})"
    fi

    echo "Finished at $(date)"
    echo ""

done < "${BATCH_FILE}"

echo "============================================================"
echo "Sequence wrapper finished at $(date)"
echo "  Total : ${TOTAL_SIMS}"
echo "  Passed: ${PASS}"
echo "  Failed: ${FAIL}"
echo "============================================================"

if [ ${FAIL} -gt 0 ]; then
    exit 1
fi
exit 0
