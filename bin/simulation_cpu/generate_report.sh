#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

RESULT_FOLDER="./results"
mkdir -p "${RESULT_FOLDER}"

LOGFILE="${RESULT_FOLDER}/logfile"

# Tulis ke logfile dan tetap tampil ke stdout untuk debugging SSH/API
exec > >(tee -a "${LOGFILE}") 2>&1

echo "=================================================="
echo "Generate report started at $(date)"
echo "SCRIPT_DIR=${SCRIPT_DIR}"
echo "PWD=$(pwd)"
echo "=================================================="

trap 'echo "ERROR: generate_report.sh failed at line ${LINENO} with exit code $?"' ERR

# Aktifkan venv yang sama seperti exec_bash.sh
# (only active this in the production machine)
#VENV_PATH="../../../../../../venv/myenv/bin/activate"
#if [ -f "${VENV_PATH}" ]; then
#  # shellcheck disable=SC1091
#  source "${VENV_PATH}"
#  echo "Venv activated from ${VENV_PATH}"
#else
#  echo "ERROR: venv activate script not found: ${VENV_PATH}"
#  exit 1
#fi

echo "Using python: $(which python || true)"
python --version || true

PARAM_FILE="param.txt"
if [ ! -f "${PARAM_FILE}" ]; then
  echo "ERROR: ${PARAM_FILE} not found in $(pwd)"
  exit 1
fi

extract_scalar_param() {
  local key="$1"
  grep "^${key}" "${PARAM_FILE}" \
    | head -n 1 \
    | cut -d'=' -f2- \
    | sed 's|//.*||' \
    | xargs
}

extract_path_param() {
  local key="$1"
  grep "^${key}" "${PARAM_FILE}" \
    | head -n 1 \
    | cut -d'=' -f2- \
    | sed 's|//.*||' \
    | xargs
}

DRUG_NAME="$(extract_scalar_param "drug_name")"
USER_NAME="$(extract_scalar_param "user_name")"
DRUG_CONCENTRATIONS="$(extract_scalar_param "drug_concentrations")"
CELL_MODEL="$(extract_scalar_param "cell_model")"
HILL_FILE="$(extract_path_param "hill_file")"

echo "DRUG_NAME=${DRUG_NAME}"
echo "USER_NAME=${USER_NAME}"
echo "DRUG_CONCENTRATIONS=${DRUG_CONCENTRATIONS}"
echo "CELL_MODEL=${CELL_MODEL}"
echo "HILL_FILE=${HILL_FILE}"

if [ -z "${DRUG_NAME}" ] || [ -z "${USER_NAME}" ] || [ -z "${CELL_MODEL}" ] || [ -z "${HILL_FILE}" ]; then
  echo "ERROR: one or more required parameters are empty"
  exit 1
fi

if [ ! -f "${HILL_FILE}" ]; then
  echo "ERROR: HILL_FILE not found: ${HILL_FILE}"
  exit 1
fi

SAMPLE_SIZE=$(( $(wc -l < "${HILL_FILE}") - 1 ))
if [ "${SAMPLE_SIZE}" -le 0 ]; then
  echo "ERROR: invalid SAMPLE_SIZE=${SAMPLE_SIZE} from HILL_FILE=${HILL_FILE}"
  exit 1
fi

echo "SAMPLE_SIZE=${SAMPLE_SIZE}"

RESULT_DRUG_PATH="${RESULT_FOLDER}"
MEDIAN_CI_FILE="${RESULT_FOLDER}/${DRUG_NAME}-median-ci.csv"

PLOT_TIME_SERIES_SCRIPT="../scripts/plot_time_series.py"
PLOT_FEATURES_SCRIPT="../scripts/plot_features.py"
COMBINE_AVG_SCRIPT="../scripts/combine_and_average_features.py"
COMPUTE_MEAN_CI_SCRIPT="../scripts/compute_mean_ci.py"
PERPLEXITY_SCRIPT="./exec_perplexity_report_generator.sh"

for required_file in \
  "${PLOT_TIME_SERIES_SCRIPT}" \
  "${PLOT_FEATURES_SCRIPT}" \
  "${COMBINE_AVG_SCRIPT}" \
  "${COMPUTE_MEAN_CI_SCRIPT}" \
  "${PERPLEXITY_SCRIPT}"
do
  if [ ! -f "${required_file}" ]; then
    echo "ERROR: required file not found: ${required_file}"
    exit 1
  fi
done

chmod +x "${PERPLEXITY_SCRIPT}" || true

# Bersihkan artifact report lama supaya validasi PDF tidak salah ambil file lama
# (removed the tex file removal since we need the LaTeX files for processing)
find "${RESULT_FOLDER}" -maxdepth 1 -type f \( \
  -name "*.pdf" -o \
  -name "*.aux" -o \
  -name "*.out" -o \
  -name "*.toc" -o \
  -name "*.log" -o \
\) -print -delete || true

mkdir -p responses

PROGRESS=0
echo "Report Progress: ${PROGRESS}%"

echo "[1/4] Generate plots from time-series result"
python3 -u "${PLOT_TIME_SERIES_SCRIPT}" \
  "${RESULT_DRUG_PATH}" \
  "${DRUG_NAME}" \
  "${SAMPLE_SIZE}" \
  "${CELL_MODEL}"
PROGRESS=25
echo "Report Progress: ${PROGRESS}%"

echo "[2/4] Unifying / plotting feature data"
python3 -u "${PLOT_FEATURES_SCRIPT}" \
  "${RESULT_DRUG_PATH}" \
  "${DRUG_NAME}" \
  "${USER_NAME}" \
  "${CELL_MODEL}"
PROGRESS=50
echo "Report Progress: ${PROGRESS}%"

echo "[3/4] Calculating average feature for each concentration"
python3 -u "${COMBINE_AVG_SCRIPT}" --base_dir "${RESULT_DRUG_PATH}"

AVG_FILE="./results/sample-averages-across-concentrations.csv"
if [ ! -f "${AVG_FILE}" ]; then
  echo "ERROR: average feature file not found: ${AVG_FILE}"
  exit 1
fi

echo "[3/4] Calculating median and 95% confidence interval"
python3 -u "${COMPUTE_MEAN_CI_SCRIPT}" \
  "${AVG_FILE}" \
  "${MEDIAN_CI_FILE}" \
  --n-boot "${SAMPLE_SIZE}"
PROGRESS=75
echo "Report Progress: ${PROGRESS}%"

echo "[4/4] Generate PDF report"
bash "${PERPLEXITY_SCRIPT}"
PROGRESS=100
echo "Report Progress: ${PROGRESS}%"

PDF_FILE="$(find "${RESULT_FOLDER}" -maxdepth 1 -type f -name '*.pdf' -printf '%T@ %p\n' 2>/dev/null | sort -nr | head -n 1 | cut -d' ' -f2- || true)"

if [ -z "${PDF_FILE}" ] || [ ! -f "${PDF_FILE}" ]; then
  echo "ERROR: report generation finished but no PDF found in ${RESULT_FOLDER}"
  exit 1
fi

echo "PDF generated: ${PDF_FILE}"
echo "Generate report finished at $(date)"

deactivate || true
exit 0
