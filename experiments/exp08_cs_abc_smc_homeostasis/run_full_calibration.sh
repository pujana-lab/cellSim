#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

BINARY="${PROJECT_ROOT}/build/run_bootstrapping"
OUTPUT_DIR="${SCRIPT_DIR}/data"

echo "=== Launching Full Production ABC-SMC Calibration (1000 cells, p in [0.047, 0.057], nu = 10*p) ==="
mkdir -p "${OUTPUT_DIR}"

python3 "${SCRIPT_DIR}/abc_smc.py" \
    --particles 200 \
    --generations 15 \
    --n-sim 200 \
    --n-cells 1000 \
    --max-healthy-cells 1000 \
    --workers 10 \
    --binary "${BINARY}" \
    --output-dir "${OUTPUT_DIR}"

echo "=== Full ABC-SMC calibration complete. Results saved in ${OUTPUT_DIR} ==="
