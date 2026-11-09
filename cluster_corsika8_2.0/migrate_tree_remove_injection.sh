#!/bin/bash
# =============================================================================
# migrate_tree_remove_injection.sh
#
# One-off migration for the injection-offset tree change (REVERSE of
# migrate_tree_add_injection.sh).
#
# The run tree level changed from
#     h${INJ_HEIGHT}_ix${INJ_X}_iy${INJ_Y}/
# back to
#     h${INJ_HEIGHT}/
# because the injection point/direction is now PREDICTED by the magnetic-field
# backtracker rather than scanned, so it is no longer part of the tree path.
#
# This script renames existing  h<number>_ix<...>_iy<...>/  directories in
# place, stripping the "_ix*_iy*" suffix:
#   h5000_ix0.0_iy0.0   -> h5000
#   h112750_ix0.0_iy0.0 -> h112750
# Bare  h<number>  directories (already migrated, or never had the suffix) are
# left alone.
#
# Usage:
#   bash migrate_tree_remove_injection.sh <BASE_PATH>            # DRY RUN
#   bash migrate_tree_remove_injection.sh <BASE_PATH> --apply    # rename
#
# Example:
#   bash migrate_tree_remove_injection.sh \
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

echo "============================================================"
echo "Tree migration: strip _ix*_iy* from h<number>_ix*_iy* directories"
echo "BASE_PATH : ${BASE_PATH}"
echo "Mode      : $([ "${DO_APPLY}" = true ] && echo APPLY || echo 'DRY RUN (no changes)')"
echo "============================================================"
echo ""

RENAMED=0
SKIPPED=0
FAILED=0

# Find directories named h<number>_ix<...> anywhere in the tree.
while IFS= read -r DIR; do
    [ -d "${DIR}" ] || continue
    BASE="$(basename "${DIR}")"

    # Strip a trailing _ix... suffix (and everything after it).
    NEW_BASE="$(printf '%s' "${BASE}" | sed -E 's/_ix.*$//')"

    if [ "${NEW_BASE}" = "${BASE}" ]; then
        continue
    fi

    PARENT="$(dirname "${DIR}")"
    TARGET="${PARENT}/${NEW_BASE}"

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
             -regex '.*/h[0-9]+_ix[^/]*' 2>/dev/null)

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
