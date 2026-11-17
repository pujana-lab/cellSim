#!/usr/bin/env bash
# run_pipeline.sh — Automated, parameterized ABC-SMC calibration and figure generation
# Canonical 4-parameter model (beta, low_delta, high_delta, neoplastic_div_rate)
# Creates a dedicated, timestamped run directory to preserve full experimental history.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Default settings (fast, compact calibration)
PARTICLES=50
GENERATIONS=8
N_SIM=200
WORKERS=10
SEED=42
TAG=""
RESUME_FROM=""
HEALTHY_MIN=0.15
HEALTHY_MAX=0.35
N_CELLS=500
TP53_RATE=0.003

print_help() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]

Runs the ABC-SMC calibration pipeline with conditioned beta_prior,
evaluates cell stage dynamics, and generates publication figures.

Options:
  -p, --particles INT       Particles per SMC generation (default: 50)
  -g, --generations INT     Number of SMC generations (default: 8)
  -s, --n-sim INT           Monte Carlo runs per particle evaluation (default: 200)
  -w, --workers INT         Worker threads (default: 10)
  --seed INT                Random seed (default: 42)
  --tag STRING              Optional label to append to run folder name
  --resume-from PATH        Resume from a previous run directory or data directory
  --healthy-min FLOAT       Minimum retention of baseline cells at age 80 (default: 0.15)
  --healthy-max FLOAT       Maximum retention of baseline cells at age 80 (default: 0.35)
  --n-cells INT             Starting cells and carrying capacity (default: 500)
  --tp53-rate FLOAT         Somatic TP53 mutation rate (default: 0.003)
  -h, --help                Show this help message
EOF
    exit 0
}

# Parse CLI arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--particles)
            PARTICLES="$2"
            shift 2
            ;;
        -g|--generations)
            GENERATIONS="$2"
            shift 2
            ;;
        -s|--n-sim)
            N_SIM="$2"
            shift 2
            ;;
        -w|--workers)
            WORKERS="$2"
            shift 2
            ;;
        --seed)
            SEED="$2"
            shift 2
            ;;
        --tag)
            TAG="_$2"
            shift 2
            ;;
        --resume-from)
            RESUME_FROM="$2"
            shift 2
            ;;
        --healthy-min)
            HEALTHY_MIN="$2"
            shift 2
            ;;
        --healthy-max)
            HEALTHY_MAX="$2"
            shift 2
            ;;
        --n-cells)
            N_CELLS="$2"
            shift 2
            ;;
        --tp53-rate)
            TP53_RATE="$2"
            shift 2
            ;;
        -h|--help)
            print_help
            ;;
        *)
            echo "Unknown option: $1" >&2
            print_help
            ;;
    esac
done

TIMESTAMP="$(date +'%Y%m%d_%H%M%S')"
RUN_NAME="run_${TIMESTAMP}${TAG}"
RUN_DIR="${SCRIPT_DIR}/runs/${RUN_NAME}"

DATA_DIR="${RUN_DIR}/data"
FIGURES_DIR="${RUN_DIR}/figures"

mkdir -p "${DATA_DIR}" "${FIGURES_DIR}"

BOOTSTRAP_BIN="${PROJECT_ROOT}/build/run_bootstrapping"
STAGE_BIN="${PROJECT_ROOT}/build/run_stage_trajectories"

for bin in "${BOOTSTRAP_BIN}" "${STAGE_BIN}"; do
    if [[ ! -x "${bin}" ]]; then
        echo "Error: executable '${bin}' not found. Please run 'cmake --build build' first." >&2
        exit 1
    fi
done

# Save run metadata
cat << EOF > "${RUN_DIR}/run_config.json"
{
  "timestamp": "${TIMESTAMP}",
  "run_name": "${RUN_NAME}",
  "particles": ${PARTICLES},
  "generations": ${GENERATIONS},
  "n_sim": ${N_SIM},
  "workers": ${WORKERS},
  "seed": ${SEED},
  "n_cells": ${N_CELLS},
  "tp53_rate": ${TP53_RATE},
  "healthy_min": ${HEALTHY_MIN},
  "healthy_max": ${HEALTHY_MAX}
}
EOF

echo "========================================================================"
echo " Starting Calibration Run: ${RUN_NAME}"
echo " Directory: ${RUN_DIR}"
echo " Parameters: particles=${PARTICLES}, generations=${GENERATIONS}, n_sim=${N_SIM}, n_cells=${N_CELLS}"
echo " Conditioned beta retention at T=80: [${HEALTHY_MIN}, ${HEALTHY_MAX}]"
echo "========================================================================"

# Step 1: Run ABC-SMC
echo ""
echo "--- [Step 1/3] Running ABC-SMC with Conditioned beta_prior ---"
EXTRA_ARGS=()
if [[ -n "${RESUME_FROM}" ]]; then
    EXTRA_ARGS+=("--resume-from" "${RESUME_FROM}")
fi

python3 "${SCRIPT_DIR}/abc_smc.py" \
    --particles "${PARTICLES}" \
    --generations "${GENERATIONS}" \
    --n-sim "${N_SIM}" \
    --n-cells "${N_CELLS}" \
    --max-healthy-cells "${N_CELLS}" \
    --tp53-rate "${TP53_RATE}" \
    --healthy-retention-min "${HEALTHY_MIN}" \
    --healthy-retention-max "${HEALTHY_MAX}" \
    --workers "${WORKERS}" \
    --seed "${SEED}" \
    --binary "${BOOTSTRAP_BIN}" \
    --output-dir "${DATA_DIR}" \
    "${EXTRA_ARGS[@]}"

# Identify final generation file
LAST_GEN_IDX=$(printf "%02d" $((GENERATIONS - 1)))
POST_CSV="${DATA_DIR}/generation_${LAST_GEN_IDX}.csv"

if [[ ! -f "${POST_CSV}" ]]; then
    POST_CSV="$(ls -1 "${DATA_DIR}"/generation_*.csv 2>/dev/null | sort | tail -n 1)"
fi

echo "Selected final posterior file: ${POST_CSV}"

# Step 2: Simulate Stage Trajectories
echo ""
echo "--- [Step 2/3] Simulating Tissue Stage Trajectories ---"
DIST_COL=$(python3 -c "import pandas as pd; df=pd.read_csv('${POST_CSV}', nrows=1); print('distance_weighted_sse' if 'distance_weighted_sse' in df.columns else 'distance')")
MAP_P=$(python3 -c "import pandas as pd; df=pd.read_csv('${POST_CSV}'); print(df.loc[df['${DIST_COL}'].idxmin()]['division_rate'])")
MAP_BETA=$(python3 -c "import pandas as pd; df=pd.read_csv('${POST_CSV}'); print(df.loc[df['${DIST_COL}'].idxmin()]['brca1_rate'])")
MAP_TP53=$(python3 -c "import pandas as pd; df=pd.read_csv('${POST_CSV}'); row=df.loc[df['${DIST_COL}'].idxmin()]; print(row['tp53_rate'] if 'tp53_rate' in row else '${TP53_RATE}')")
MAP_LOW_DELTA=$(python3 -c "import pandas as pd; df=pd.read_csv('${POST_CSV}'); print(df.loc[df['${DIST_COL}'].idxmin()]['low_delta'])")
MAP_HIGH_DELTA=$(python3 -c "import pandas as pd; df=pd.read_csv('${POST_CSV}'); print(df.loc[df['${DIST_COL}'].idxmin()]['high_delta'])")
MAP_NU=$(python3 -c "import pandas as pd; df=pd.read_csv('${POST_CSV}'); print(df.loc[df['${DIST_COL}'].idxmin()]['neoplastic_div_rate'])")

echo "MAP posterior parameters:"
echo "  p (basal division)     = ${MAP_P}"
echo "  beta (BRCA1 rate)      = ${MAP_BETA}"
echo "  tp53_rate              = ${MAP_TP53}"
echo "  low_delta              = ${MAP_LOW_DELTA}"
echo "  high_delta             = ${MAP_HIGH_DELTA}"
echo "  nu (neoplastic div)    = ${MAP_NU}"

STAGE_CSV="${DATA_DIR}/stage_trajectories_100tissues.csv"
"${STAGE_BIN}" \
    --n_tissues 100 \
    --n_cells "${N_CELLS}" \
    --max_cells "${N_CELLS}" \
    --max_t 80 \
    --division_rate "${MAP_P}" \
    --brca1_rate "${MAP_BETA}" \
    --tp53_rate "${MAP_TP53}" \
    --delta_low "${MAP_LOW_DELTA}" \
    --delta_high "${MAP_HIGH_DELTA}" \
    --neoplastic_div_rate "${MAP_NU}" \
    --output "${STAGE_CSV}" \
    --seed "${SEED}"

# Step 3: Render Figures
echo ""
echo "--- [Step 3/3] Rendering Figures ---"
python3 "${SCRIPT_DIR}/plot_figure3.py" \
    --data-dir "${DATA_DIR}" \
    --output-dir "${FIGURES_DIR}"

python3 "${SCRIPT_DIR}/plot_sfig4_homeostasis.py" \
    --data-dir "${DATA_DIR}" \
    --output-dir "${FIGURES_DIR}" \
    --variant all

python3 "${SCRIPT_DIR}/plot_cell_stage_dynamics.py" \
    --csv "${STAGE_CSV}" \
    --out-dir "${FIGURES_DIR}"

# Create symlink to latest run
ln -sfn "${RUN_NAME}" "${SCRIPT_DIR}/runs/latest"

echo ""
echo "========================================================================"
echo " Run completed successfully!"
echo " Results and plots saved to: ${RUN_DIR}"
echo " Figures available at: ${FIGURES_DIR}"
echo "========================================================================"
