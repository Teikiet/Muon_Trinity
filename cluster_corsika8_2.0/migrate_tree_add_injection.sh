#!/bin/bash
# =============================================================================
# migrate_tree_add_injection.sh
#
# One-off migration for the injection-offset tree change.
#
# The run tree level changed from
#     h${INJ_HEIGHT}/
# to
#     h${INJ_HEIGHT}_ix${INJ_X}_iy${INJ_Y}/
# so every pre-existing run directory is now invisible to the updated path
# builders (save_CARE2csv_chunk_tree.py, adaptive_controller.py).
#
# This script renames existing  h<number>/  directories in place, adding the
# suffix "_ix0.0_iy0.0" (i.e. declaring them as injection-free runs).
#
# It only touches directories whose name matches  h<number>  exactly:
#   h112750      -> h112750_ix0.0_iy0.0
#   h2944        -> h2944_ix0.0_iy0.0
# Directories that already carry an _ix suffix are left alone.
#
# Usage:
#   bash migrate_tree_add_injection.sh <BASE_PATH>            # DRY RUN (default)
#   bash migrate_tree_add_injection.sh <BASE_PATH> --apply    # actually rename
#
# Example:
#   bash migrate_tree_add_injection.sh \
#       /scratch/general/vast/u1520754/muon_sim_chain_tree_adaptive
# =============================================================================

set -euo pipefail

BASE_PATH="${1:-}"
APPLY="${2:-}"

if [ -z "${BASE_PATH}" ]; then
    echo "ERROR: missing <BASE_PATH>."
    echo "Usage: bash $0 <BASE_PATH> [--apply]"
    exit 1
fi

if [ ! -d "${BASE_PATH}" ]; then
    echo "ERROR: not a directory: ${BASE_PATH}"
    exit 1
fi

DO_APPLY=false
if [ "${APPLY}" = "--apply" ]; then
    DO_APPLY=true
elif [ -n "${APPLY}" ]; then
    echo "ERROR: unknown option '${APPLY}'. Use --apply or nothing (dry run)."
    exit 1
fi

SUFFIX="_ix0.0_iy0.0"

echo "============================================================"
echo "Tree migration: add ${SUFFIX} to bare h<number> directories"
echo "BASE_PATH : ${BASE_PATH}"
echo "Mode      : $([ "${DO_APPLY}" = true ] && echo APPLY || echo 'DRY RUN (no changes)')"
echo "============================================================"
echo ""

RENAMED=0
SKIPPED=0
FAILED=0

# Find directories named exactly h<digits> or h<digits>.<digits>, anywhere in the tree.
while IFS= read -r DIR; do
    [ -d "${DIR}" ] || continue
    BASE="$(basename "${DIR}")"

    # Only bare h<number> levels (no existing _ix/_iy suffix)
    if ! [[ "${BASE}" =~ ^h[0-9]+([.][0-9]+)?$ ]]; then
        continue
    fi

    TARGET="${DIR}${SUFFIX}"

    if [ -e "${TARGET}" ]; then
        echo "SKIP  (target exists): ${DIR} -> ${TARGET}"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    if ${DO_APPLY}; then
        if mv "${DIR}" "${TARGET}"; then
            echo "RENAMED: ${DIR} -> ${TARGET}"
            RENAMED=$((RENAMED + 1))
        else
            echo "FAILED : ${DIR}"
            FAILED=$((FAILED + 1))
        fi
    else
        echo "WOULD RENAME: ${DIR} -> ${TARGET}"
        RENAMED=$((RENAMED + 1))
    fi
done < <(find "${BASE_PATH}" -type d -regextype posix-extended \
             -regex '.*/h[0-9]+(\.[0-9]+)?' 2>/dev/null)

echo ""
echo "============================================================"
echo "$([ "${DO_APPLY}" = true ] && echo Renamed || echo 'Would rename'): ${RENAMED}"
echo "Skipped : ${SKIPPED}"
echo "Failed  : ${FAILED}"
if [ "${DO_APPLY}" != true ]; then
    echo ""
    echo "This was a DRY RUN. Re-run with --apply to perform the renames."
fi
echo "============================================================"
