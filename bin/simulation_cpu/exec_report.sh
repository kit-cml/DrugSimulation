#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

RESULT_FOLDER="./results"
mkdir -p "${RESULT_FOLDER}"

LOGFILE="${RESULT_FOLDER}/logfile"

echo "==================================================" >> "${LOGFILE}" 2>&1
echo "Exec generate started at $(date)" >> "${LOGFILE}" 2>&1
echo "SCRIPT_DIR=${SCRIPT_DIR}" >> "${LOGFILE}" 2>&1
echo "PWD=$(pwd)" >> "${LOGFILE}" 2>&1
echo "==================================================" >> "${LOGFILE}" 2>&1

VENV_PATH="../../../../../../venv/myenv/bin/activate"

if [ -f "${VENV_PATH}" ]; then
  # shellcheck disable=SC1091
  source "${VENV_PATH}"
  echo "Venv activated from ${VENV_PATH}" >> "${LOGFILE}" 2>&1
else
  echo "ERROR: venv activate script not found: ${VENV_PATH}" >> "${LOGFILE}" 2>&1
  exit 1
fi

chmod +x "./generate_report.sh" || true

bash "./generate_report.sh" >> "${LOGFILE}" 2>&1
EXIT_CODE=$?

if [ "${EXIT_CODE}" -ne 0 ]; then
  echo "ERROR: generate_report.sh failed with exit code ${EXIT_CODE}" >> "${LOGFILE}" 2>&1
  deactivate || true
  exit "${EXIT_CODE}"
fi

PDF_FILE="$(find "${RESULT_FOLDER}" -maxdepth 1 -type f -name '*.pdf' | head -n 1 || true)"

if [ -z "${PDF_FILE}" ]; then
  echo "ERROR: generate_report.sh finished but no PDF found in ${RESULT_FOLDER}" >> "${LOGFILE}" 2>&1
  deactivate || true
  exit 1
fi

echo "PDF generated: ${PDF_FILE}" >> "${LOGFILE}" 2>&1
echo "Exec generate finished at $(date)" >> "${LOGFILE}" 2>&1

deactivate || true
exit 0
