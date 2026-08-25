#!/usr/bin/env bash
set -euo pipefail

# Tunables (can be overridden by environment variables)
VARIANT_FLAG="${VARIANT_FLAG:--DCIPAORDV1}"
SIM_DIR="${SIM_DIR:-bin/simulation_cpu}"
NP="${NP:-10}"
RUN_SIM="${RUN_SIM:-1}"
PARAM_FILE="${PARAM_FILE:-param.txt}"
PACING_OVERRIDE="${PACING_OVERRIDE:-}"
PACING_WRITE_OVERRIDE="${PACING_WRITE_OVERRIDE:-}"
SOLVER_OVERRIDE="${SOLVER_OVERRIDE:-}"
TIME_STEP_MIN_OVERRIDE="${TIME_STEP_MIN_OVERRIDE:-}"
HILL_FILE_OVERRIDE="${HILL_FILE_OVERRIDE:-}"
HERG_FILE_OVERRIDE="${HERG_FILE_OVERRIDE:-}"
DRUG_NAME_OVERRIDE="${DRUG_NAME_OVERRIDE:-}"
CONC_OVERRIDE="${CONC_OVERRIDE:-}"
# Extra compiler flags appended after the Makefile's own BASE_CXXFLAGS
# (-std=c++11 -g -O0). Empty by default so plain/manual builds keep that
# debug build unchanged; callers that want an optimized binary (e.g. the
# benchmark scripts) set this to something like "-O2 -DNDEBUG -g0" -- since
# gcc/clang use the last -O/-g flag seen on the command line, this cleanly
# overrides the debug build's -O0 -g.
BUILD_OPT_FLAGS="${BUILD_OPT_FLAGS:-}"

variant_suffix_from_flag() {
  case "$1" in
    "")
      echo "ORd-static"
      ;;
    "-DCIPAORDV1")
      echo "CiPAORdv1.0"
      ;;
    "-DTOR_ORD")
      echo "ToR-ORd"
      ;;
    "-DTOR_ORD_DYNCL")
      echo "ToR-ORd-dynCl"
      ;;
    "-DGRANDI")
      echo "Grandi"
      ;;
    "-DORD_STATIC_LAND")
      echo "ORd-static_Land"
      ;;
    "-DCIPAORDV1_LAND")
      echo "CiPAORdv1.0_Land"
      ;;
    "-DTOR_ORD_LAND")
      echo "ToR-ORd_Land"
      ;;
    *)
      echo ""
      ;;
  esac
}

BIN_SUFFIX="$(variant_suffix_from_flag "$VARIANT_FLAG")"
if [[ -z "$BIN_SUFFIX" ]]; then
  echo "Unsupported VARIANT_FLAG: $VARIANT_FLAG"
  exit 2
fi

echo "[1/3] Building libCML static library in Linux container..."
cd /workspace/libCML
make clean || true
make build/libcml.a

echo "[2/3] Building DrugSimulation binary for variant: $BIN_SUFFIX"
cd /workspace/DrugSimulation
make clean || true
make all \
  EXTRA_CXXFLAGS+="-std=c++11 $VARIANT_FLAG $BUILD_OPT_FLAGS"

if [[ "$RUN_SIM" != "1" ]]; then
  echo "RUN_SIM=$RUN_SIM, skipping simulation run."
  exit 0
fi

echo "[3/3] Running simulation with mpiexec in $SIM_DIR"
cd "/workspace/DrugSimulation/$SIM_DIR"

PARAM_TO_USE="$PARAM_FILE"
TEMP_PARAM_FILE=""

if [[ -n "$PACING_OVERRIDE" || -n "$PACING_WRITE_OVERRIDE" || -n "$SOLVER_OVERRIDE" || -n "$TIME_STEP_MIN_OVERRIDE" || -n "$HILL_FILE_OVERRIDE" || -n "$HERG_FILE_OVERRIDE" || -n "$DRUG_NAME_OVERRIDE" || -n "$CONC_OVERRIDE" ]]; then
  # PID-suffixed so concurrent invocations (e.g. multiple `docker compose exec`
  # calls running different drugs in parallel against this same bind-mounted
  # SIM_DIR) never race on the same file -- a fixed name here previously let
  # one invocation's sed edits clobber another's mid-write.
  TEMP_PARAM_FILE=".docker-test-param.$$.txt"
  cp "$PARAM_FILE" "$TEMP_PARAM_FILE"

  if [[ -n "$PACING_OVERRIDE" ]]; then
    sed -E -i "s#^(number_pacing[[:space:]]*=[[:space:]]*)[0-9]+#\1${PACING_OVERRIDE}#" "$TEMP_PARAM_FILE"
  fi

  if [[ -n "$PACING_WRITE_OVERRIDE" ]]; then
    sed -E -i "s#^(number_pacing_write[[:space:]]*=[[:space:]]*)[0-9]+#\1${PACING_WRITE_OVERRIDE}#" "$TEMP_PARAM_FILE"
  elif [[ -n "$PACING_OVERRIDE" ]]; then
    sed -E -i "s#^(number_pacing_write[[:space:]]*=[[:space:]]*)[0-9]+#\1${PACING_OVERRIDE}#" "$TEMP_PARAM_FILE"
  fi

  if [[ -n "$SOLVER_OVERRIDE" ]]; then
    sed -E -i "s#^(solver_type[[:space:]]*=[[:space:]]*).*	*//?#\1${SOLVER_OVERRIDE} // #" "$TEMP_PARAM_FILE"
  fi

  if [[ -n "$TIME_STEP_MIN_OVERRIDE" ]]; then
    # Euler uses this as a fixed integration step for the whole run (time_step_max is ignored
    # by the Euler branch), so this needs to be small (e.g. 0.005) to stay numerically stable.
    sed -E -i "s#^(time_step_min[[:space:]]*=[[:space:]]*)[0-9.]+#\1${TIME_STEP_MIN_OVERRIDE}#" "$TEMP_PARAM_FILE"
  fi

  if [[ -n "$HILL_FILE_OVERRIDE" ]]; then
    sed -E -i "s#^(hill_file[[:space:]]*=[[:space:]]*)\S+#\1${HILL_FILE_OVERRIDE}#" "$TEMP_PARAM_FILE"
  fi

  if [[ -n "$HERG_FILE_OVERRIDE" ]]; then
    sed -E -i "s#^(herg_file[[:space:]]*=[[:space:]]*)\S+#\1${HERG_FILE_OVERRIDE}#" "$TEMP_PARAM_FILE"
  fi

  if [[ -n "$DRUG_NAME_OVERRIDE" ]]; then
    sed -E -i "s#^(drug_name[[:space:]]*=[[:space:]]*)\S+#\1${DRUG_NAME_OVERRIDE}#" "$TEMP_PARAM_FILE"
  fi

  if [[ -n "$CONC_OVERRIDE" ]]; then
    sed -E -i "s#^(drug_concentrations[[:space:]]*=[[:space:]]*)\S+#\1${CONC_OVERRIDE}#" "$TEMP_PARAM_FILE"
  fi

  PARAM_TO_USE="$TEMP_PARAM_FILE"
  trap 'rm -f "$TEMP_PARAM_FILE"' EXIT
fi

DRUG_CONCENTRATIONS_RAW="$(grep '^drug_concentrations' "$PARAM_TO_USE" | cut -d'=' -f2 | sed 's#//.*##' | xargs || true)"
mkdir -p results/0.00
if [[ -n "$DRUG_CONCENTRATIONS_RAW" ]]; then
  IFS=',' read -r -a concs <<< "$DRUG_CONCENTRATIONS_RAW"
  for conc in "${concs[@]}"; do
    conc_trimmed="$(echo "$conc" | xargs)"
    if [[ -n "$conc_trimmed" ]]; then
      conc_dir="$(printf "%.2f" "$conc_trimmed")"
      mkdir -p "results/$conc_dir"
    fi
  done
fi

echo "Using parameter file: $PARAM_TO_USE"
mpiexec --allow-run-as-root -np "$NP" "../drugsim_${BIN_SUFFIX}" -input_deck "$PARAM_TO_USE" \
  | tee results/logfile.docker

echo "Simulation finished. Logs: $SIM_DIR/results/logfile.docker"
