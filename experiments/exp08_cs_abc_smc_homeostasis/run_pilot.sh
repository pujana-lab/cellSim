#!/usr/bin/env bash
set -euo pipefail

# Script to run a 4-generation proof-of-concept pilot run of exp08 ABC-SMC calibration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

BINARY="${PROJECT_ROOT}/build/run_bootstrapping"
OUTPUT_DIR="${SCRIPT_DIR}/results_pilot"

echo "=== Running exp08 ABC-SMC Pilot Calibration (4 generations, 20 particles) ==="
python3 "${SCRIPT_DIR}/abc_smc.py" \
    --particles 20 \
    --generations 4 \
    --n-sim 50 \
    --n-cells 1000 \
    --max-healthy-cells 1000 \
    --workers 4 \
    --binary "${BINARY}" \
    --output-dir "${OUTPUT_DIR}"

echo "=== Pilot calibration complete. Results saved in ${OUTPUT_DIR} ==="
