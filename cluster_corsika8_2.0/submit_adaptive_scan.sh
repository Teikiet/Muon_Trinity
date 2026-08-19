#!/bin/bash

set -u

SCRIPT_DIR="$HOME/Muon_Trinity/cluster_corsika8_2.0"
INPUT_SCRIPT="${SCRIPT_DIR}/Trinity_adaptive_input.py"
CONTROLLER_SCRIPT="${SCRIPT_DIR}/adaptive_controller.py"

ENERGY=""
N_INIT=""
N_MAX=""
MIN_TRIGGERS=""
H_START=""
H_MAX=""
DELTA_H=""
H_TERMINATE_ABOVE=""
H_BATCH=""
ZEN_START=""
ZEN_MAX=""
DELTA_ZEN=""
AZ=""
TRIGGER_PE=""
MAX_QUEUED=""
SAVE_BATCHES=""
H_VALUES=""
RESUME="false"
RERUN_FAILED="false"
DRY_RUN="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --energy)
            if [[ $# -gt 1 && ! "$2" =~ ^-- ]]; then
                ENERGY="$2"
                shift 2
            else
                echo "ERROR: --energy requires a value"
                exit 1
            fi
            ;;
        --energy=*)
            ENERGY="${1#*=}"
            shift
            ;;
        --n-init)
            N_INIT="$2"
            shift 2
            ;;
        --n-init=*)
            N_INIT="${1#*=}"
            shift
            ;;
        --n-max)
            N_MAX="$2"
            shift 2
            ;;
        --n-max=*)
            N_MAX="${1#*=}"
            shift
            ;;
        --min-triggers)
            MIN_TRIGGERS="$2"
            shift 2
            ;;
        --min-triggers=*)
            MIN_TRIGGERS="${1#*=}"
            shift
            ;;
        --h-start)
            H_START="$2"
            shift 2
            ;;
        --h-start=*)
            H_START="${1#*=}"
            shift
            ;;
        --h-max)
            H_MAX="$2"
            shift 2
            ;;
        --h-max=*)
            H_MAX="${1#*=}"
            shift
            ;;
        --delta-h)
            DELTA_H="$2"
            shift 2
            ;;
        --delta-h=*)
            DELTA_H="${1#*=}"
            shift
            ;;
        --h-terminate-above)
            H_TERMINATE_ABOVE="$2"
            shift 2
            ;;
        --h-terminate-above=*)
            H_TERMINATE_ABOVE="${1#*=}"
            shift
            ;;
        --h-batch)
            H_BATCH="$2"
            shift 2
            ;;
        --h-batch=*)
            H_BATCH="${1#*=}"
            shift
            ;;
        --zen-start)
            ZEN_START="$2"
            shift 2
            ;;
        --zen-start=*)
            ZEN_START="${1#*=}"
            shift
            ;;
        --zen-max)
            ZEN_MAX="$2"
            shift 2
            ;;
        --zen-max=*)
            ZEN_MAX="${1#*=}"
            shift
            ;;
        --h-values)
            H_VALUES="$2"
            shift 2
            ;;
        --h-values=*)
            H_VALUES="${1#*=}"
            shift
            ;;
        --delta-zen)
            DELTA_ZEN="$2"
            shift 2
            ;;
        --delta-zen=*)
            DELTA_ZEN="${1#*=}"
            shift
            ;;
        --az)
            AZ="$2"
            shift 2
            ;;
        --az=*)
            AZ="${1#*=}"
            shift
            ;;
        --trigger-pe)
            TRIGGER_PE="$2"
            shift 2
            ;;
        --trigger-pe=*)
            TRIGGER_PE="${1#*=}"
            shift
            ;;
        --max-queued)
            MAX_QUEUED="$2"
            shift 2
            ;;
        --max-queued=*)
            MAX_QUEUED="${1#*=}"
            shift
            ;;
        --save-batches)
            SAVE_BATCHES="$2"
            shift 2
            ;;
        --save-batches=*)
            SAVE_BATCHES="${1#*=}"
            shift
            ;;
        --resume)
            RESUME="true"
            shift
            ;;
        --rerun-failed)
            RERUN_FAILED="true"
            shift
            ;;
        --dry-run)
            DRY_RUN="true"
            shift
            ;;
        --help|-h)
            echo "Usage: bash submit_adaptive_scan.sh --energy E [--n-init 1] [--n-max 64] [--min-triggers 1] [--h-start 3000] [--h-max 50000] [--delta-h 200] [--h-terminate-above 20000] [--h-batch 1] [--zen-start 90.0] [--zen-max 92.0] [--delta-zen 0.3] [--az 270.0] [--trigger-pe 20] [--max-queued 900] [--save-batches 10] [--h-values 3000,4000,5000] [--resume] [--rerun-failed] [--dry-run]"
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option '$1'"
            exit 1
            ;;
    esac
done

if [[ -z "${ENERGY}" ]]; then
    echo "ERROR: --energy is required"
    exit 1
fi

INPUT_JSON=$(mktemp /tmp/trinity_adaptive_input_XXXXXX.json)
trap 'rm -f "${INPUT_JSON}"' EXIT
python3 "${INPUT_SCRIPT}" --output "${INPUT_JSON}"

read -r PDG TEL_RADIUS TEL_Y ENERGY_CSV DEFAULTS_JSON <<< "$(python3 - "${INPUT_JSON}" <<'PYEOF'
import json
import sys

with open(sys.argv[1], "r", encoding="utf-8") as handle:
    data = json.load(handle)

pdg = data.get("pdg", 13)
tel_radius = data.get("tel_radius", 5)
tel_y = data.get("tel_y", 0)
energies = data.get("energy_strings", [])
defaults = data.get("adaptive_defaults", {})

print(f"{pdg} {tel_radius} {tel_y} {','.join(str(v) for v in energies)} {json.dumps(defaults)}")
PYEOF
)"

IFS="," read -ra ENERGY_STRS <<< "${ENERGY_CSV}"
FOUND_ENERGY="false"
for ENERGY_STR in "${ENERGY_STRS[@]}"; do
    if [[ "${ENERGY_STR}" == "${ENERGY}" ]]; then
        FOUND_ENERGY="true"
        break
    fi
done

if [[ "${FOUND_ENERGY}" != "true" ]]; then
    echo "ERROR: requested energy '${ENERGY}' is not in the adaptive energy list"
    exit 1
fi

export DEFAULTS_JSON

if [[ -z "${N_INIT}" ]]; then N_INIT=$(python3 -c 'import json, os; print(json.loads(os.environ["DEFAULTS_JSON"]) ["n_init"])'); fi
if [[ -z "${N_MAX}" ]]; then N_MAX=$(python3 -c 'import json, os; print(json.loads(os.environ["DEFAULTS_JSON"]) ["n_max"])'); fi
if [[ -z "${MIN_TRIGGERS}" ]]; then MIN_TRIGGERS=$(python3 -c 'import json, os; print(json.loads(os.environ["DEFAULTS_JSON"]).get("min_triggers", 1))'); fi
if [[ -z "${H_START}" ]]; then H_START=$(python3 -c 'import json, os; print(json.loads(os.environ["DEFAULTS_JSON"]) ["h_start"])'); fi
if [[ -z "${H_MAX}" ]]; then H_MAX=$(python3 -c 'import json, os; print(json.loads(os.environ["DEFAULTS_JSON"]) ["h_max"])'); fi
if [[ -z "${DELTA_H}" ]]; then DELTA_H=$(python3 -c 'import json, os; print(json.loads(os.environ["DEFAULTS_JSON"]) ["delta_h"])'); fi
if [[ -z "${H_TERMINATE_ABOVE}" ]]; then H_TERMINATE_ABOVE=$(python3 -c 'import json, os; print(json.loads(os.environ["DEFAULTS_JSON"]) ["h_terminate_above"])'); fi
if [[ -z "${H_BATCH}" ]]; then H_BATCH=$(python3 -c 'import json, os; print(json.loads(os.environ["DEFAULTS_JSON"]) ["h_batch"])'); fi
if [[ -z "${ZEN_START}" ]]; then ZEN_START=$(python3 -c 'import json, os; print(json.loads(os.environ["DEFAULTS_JSON"]) ["zen_start"])'); fi
if [[ -z "${ZEN_MAX}" ]]; then ZEN_MAX=$(python3 -c 'import json, os; print(json.loads(os.environ["DEFAULTS_JSON"]).get("zen_max", 92.0))'); fi
if [[ -z "${DELTA_ZEN}" ]]; then DELTA_ZEN=$(python3 -c 'import json, os; print(json.loads(os.environ["DEFAULTS_JSON"]) ["delta_zen"])'); fi
if [[ -z "${AZ}" ]]; then AZ=$(python3 -c 'import json, os; print(json.loads(os.environ["DEFAULTS_JSON"]) ["az"])'); fi
if [[ -z "${TRIGGER_PE}" ]]; then TRIGGER_PE=$(python3 -c 'import json, os; print(json.loads(os.environ["DEFAULTS_JSON"]) ["trigger_pe"])'); fi
if [[ -z "${MAX_QUEUED}" ]]; then MAX_QUEUED=$(python3 -c 'import json, os; print(json.loads(os.environ["DEFAULTS_JSON"]) ["max_queued"])'); fi
if [[ -z "${SAVE_BATCHES}" ]]; then SAVE_BATCHES=$(python3 -c 'import json, os; print(json.loads(os.environ["DEFAULTS_JSON"]).get("save_batches", 1))'); fi
if [[ -z "${H_VALUES}" ]]; then H_VALUES=$(python3 -c 'import json, os; print(",".join(str(v) for v in json.loads(os.environ["DEFAULTS_JSON"]).get("h_values", [])))'); fi

CONTROLLER_ARGS=(
    "${CONTROLLER_SCRIPT}"
    --energy "${ENERGY}"
    --n-init "${N_INIT}"
    --n-max "${N_MAX}"
    --min-triggers "${MIN_TRIGGERS}"
    --h-start "${H_START}"
    --h-max "${H_MAX}"
    --delta-h "${DELTA_H}"
    --h-terminate-above "${H_TERMINATE_ABOVE}"
    --h-batch "${H_BATCH}"
    --zen-start "${ZEN_START}"
    --zen-max "${ZEN_MAX}"
    --delta-zen "${DELTA_ZEN}"
    --az "${AZ}"
    --trigger-pe "${TRIGGER_PE}"
    --max-queued "${MAX_QUEUED}"
    --save-batches "${SAVE_BATCHES}"
    --pid "${PDG}"
    --radius "${TEL_RADIUS}"
    --tel-y "${TEL_Y}"
    --obs-level 2944
    --hadron-model SIBYLL-2.3d
)

if [[ -n "${H_VALUES}" ]]; then
    CONTROLLER_ARGS+=(--h-values "${H_VALUES}")
fi

if [[ "${RESUME}" == "true" ]]; then
    CONTROLLER_ARGS+=(--resume)
fi
if [[ "${RERUN_FAILED}" == "true" ]]; then
    CONTROLLER_ARGS+=(--rerun-failed)
fi
if [[ "${DRY_RUN}" == "true" ]]; then
    CONTROLLER_ARGS+=(--dry-run)
fi

LOG_DIR="$HOME/csv_logs/adaptive_E${ENERGY}"
mkdir -p "${LOG_DIR}"

SBATCH_CMD=(
    sbatch
    --parsable
    --account=owner-guest
    --partition=kingspeak-guest
    --time=0:15:00
    --mem=2G
    --job-name="adaptive_controller_E${ENERGY}"
    --output="${LOG_DIR}/cycle_%j.out"
    --error="${LOG_DIR}/cycle_%j.err"
    --wrap
    "source ${HOME}/miniconda3/etc/profile.d/conda.sh; conda activate jupyter_env; python3 ${CONTROLLER_ARGS[*]}"
)

if [[ "${DRY_RUN}" == "true" ]]; then
    printf 'DRY-RUN: '
    printf '%q ' "${SBATCH_CMD[@]}"
    printf '\n'
    exit 0
fi

"${SBATCH_CMD[@]}"