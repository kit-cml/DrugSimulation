#!/usr/bin/env bash
# Hardware-performance benchmark harness for DrugSimulation.
#
# Sweeps the drugsim container across drugs x solvers x MPI process counts,
# holding every other hyperparameter fixed, and records wall-clock elapsed
# time for each run to a CSV so runs are comparable across hardware. Every
# row is tagged with host + container hardware info (CPU, RAM, arch, whether
# Docker is emulating a foreign architecture) so CSVs from different machines
# can be compared or concatenated without losing track of what ran where.
#
# Run from the repository root (or anywhere; it cd's to the repo root itself):
#   bash docker/linux/benchmark.sh
#
# Everything is configured via environment variables (same style as
# build_and_run.sh / run_only.sh). See `bash docker/linux/benchmark.sh -h`.
set -u -o pipefail

usage() {
  cat <<'EOF'
Usage: [ENV=value ...] bash docker/linux/benchmark.sh

Sweeps DRUGS x SOLVERS x NP_LIST, running each combination REPEATS times via
`docker compose exec ... run_only.sh`, and appends one row per run to a CSV
under OUT_DIR. All other simulation hyperparameters are held fixed across the
whole sweep so timings are comparable.

Environment variables (defaults shown):
  DRUGS              (auto-discovered from chantest_hill/*)  comma/space separated
  SOLVERS             CVode,Euler                             comma/space separated
  NP_LIST              1                                       comma/space separated MPI process counts
  REPEATS               1                                       repeats per combination
  VARIANT_FLAG        -DCIPAORDV1                             cell-model build variant (fixed for the whole sweep)
  PACING                5                                       number_pacing override
  PACING_WRITE          5                                       number_pacing_write override
  TIME_STEP_MIN         0.001                                   time_step_min override
  CONC                (unset -> use param.txt's own concentrations for every drug)
  HILL_SAMPLES          30                                      -> chantest_hill/<drug>/IC50_samples<HILL_SAMPLES>.csv
  HERG_SAMPLES          10                                      -> chantest_herg/<drug>/boot_pars<HERG_SAMPLES>.csv
  SIM_DIR              bin/simulation_cpu
  SKIP_BUILD             0                                       set 1 to reuse an already-built binary
  FAIL_FAST              0                                       set 1 to stop the sweep on first failure
  OUT_DIR              benchmark_results/<timestamp>            where results.csv/summary.txt/logs/ are written
  DRY_RUN                0                                       set 1 to print the commands without running them

Hardware/environment info (host OS/arch/CPU/RAM, container OS/arch/CPU/RAM,
whether Docker is running inside a VM, whether the container is CPU-emulated)
is auto-detected once per sweep, printed up front, written to
OUT_DIR/environment.txt, and appended as extra columns on every results.csv
row -- no configuration needed.

Example - sweep all drugs across both solvers at 1 and 4 MPI ranks, 3 repeats:
  NP_LIST="1,4" REPEATS=3 bash docker/linux/benchmark.sh

Example - just two drugs, Euler only, for a quick sanity check:
  DRUGS="quinidine,bepridil" SOLVERS="Euler" bash docker/linux/benchmark.sh
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

REPO_ROOT="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "$REPO_ROOT"

SIM_DIR="${SIM_DIR:-bin/simulation_cpu}"
VARIANT_FLAG="${VARIANT_FLAG:--DCIPAORDV1}"
SOLVERS="${SOLVERS:-CVode,Euler}"
NP_LIST="${NP_LIST:-1}"
REPEATS="${REPEATS:-1}"
PACING="${PACING:-5}"
PACING_WRITE="${PACING_WRITE:-5}"
TIME_STEP_MIN="${TIME_STEP_MIN:-0.001}"
CONC="${CONC:-}"
HILL_SAMPLES="${HILL_SAMPLES:-30}"
HERG_SAMPLES="${HERG_SAMPLES:-10}"
SKIP_BUILD="${SKIP_BUILD:-0}"
FAIL_FAST="${FAIL_FAST:-0}"
DRY_RUN="${DRY_RUN:-0}"
OUT_DIR="${OUT_DIR:-benchmark_results/$(date +%Y%m%d_%H%M%S)}"

# --- Hardware / environment detection --------------------------------------
# Captured once per sweep (hardware doesn't change mid-run) and stamped onto
# every CSV row, so results from different machines stay self-describing.
norm_arch() {
  case "$1" in
    arm64 | aarch64) echo "arm64" ;;
    x86_64 | amd64) echo "x86_64" ;;
    *) echo "$1" ;;
  esac
}

gather_host_info() {
  HOST_OS="$(uname -s)"
  HOST_ARCH="$(uname -m)"
  case "$HOST_OS" in
    Darwin)
      HOST_CPU_MODEL="$(sysctl -n machdep.cpu.brand_string 2>/dev/null)"
      HOST_CPU_CORES="$(sysctl -n hw.logicalcpu 2>/dev/null)"
      local ram_bytes
      ram_bytes="$(sysctl -n hw.memsize 2>/dev/null)"
      HOST_RAM_GB="$(awk -v b="${ram_bytes:-0}" 'BEGIN { printf "%.1f", b / 1024 / 1024 / 1024 }')"
      ;;
    Linux)
      HOST_CPU_MODEL="$(awk -F: '/model name/ { gsub(/^ /, "", $2); print $2; exit }' /proc/cpuinfo 2>/dev/null)"
      HOST_CPU_CORES="$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null)"
      local ram_kb
      ram_kb="$(awk '/MemTotal/ { print $2 }' /proc/meminfo 2>/dev/null)"
      HOST_RAM_GB="$(awk -v kb="${ram_kb:-0}" 'BEGIN { printf "%.1f", kb / 1024 / 1024 }')"
      ;;
    *)
      HOST_CPU_MODEL="unknown"
      HOST_CPU_CORES="unknown"
      HOST_RAM_GB="unknown"
      ;;
  esac
  HOST_CPU_MODEL="${HOST_CPU_MODEL:-unknown}"
  HOST_CPU_MODEL="${HOST_CPU_MODEL//,/}" # keep CSV columns intact
  HOST_CPU_CORES="${HOST_CPU_CORES:-unknown}"
  HOST_RAM_GB="${HOST_RAM_GB:-unknown}"
}

# Requires the drugsim-linux container to already be running.
gather_container_info() {
  DOCKER_OS="$(docker info --format '{{.OperatingSystem}}' 2>/dev/null)"
  if [[ "$DOCKER_OS" == *"Docker Desktop"* ]]; then
    DOCKER_RUNTIME="docker-desktop-vm"
  else
    DOCKER_RUNTIME="linux-native"
  fi

  CONTAINER_ARCH="$(docker compose exec -T drugsim-linux uname -m 2>/dev/null | tr -d '\r\n')"
  CONTAINER_CPU_CORES="$(docker compose exec -T drugsim-linux nproc 2>/dev/null | tr -d '\r\n')"
  local mem_kb
  mem_kb="$(docker compose exec -T drugsim-linux awk '/MemTotal/ { print $2 }' /proc/meminfo 2>/dev/null | tr -d '\r\n')"
  CONTAINER_RAM_GB="$(awk -v kb="${mem_kb:-0}" 'BEGIN { printf "%.1f", kb / 1024 / 1024 }')"

  if [[ -n "$CONTAINER_ARCH" && "$(norm_arch "$HOST_ARCH")" != "$(norm_arch "$CONTAINER_ARCH")" ]]; then
    EMULATED="yes"
  else
    EMULATED="no"
  fi

  CONTAINER_ARCH="${CONTAINER_ARCH:-unknown}"
  CONTAINER_CPU_CORES="${CONTAINER_CPU_CORES:-unknown}"
  CONTAINER_RAM_GB="${CONTAINER_RAM_GB:-unknown}"
}

# Placeholders so results.csv columns line up even in DRY_RUN mode, where the
# container is never contacted.
DOCKER_RUNTIME="n/a (dry run)"
CONTAINER_ARCH="n/a"
CONTAINER_CPU_CORES="n/a"
CONTAINER_RAM_GB="n/a"
EMULATED="n/a"

gather_host_info

if [[ -z "${DRUGS:-}" ]]; then
  DRUGS="$(find "$SIM_DIR/chantest_hill" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort | tr '\n' ',')"
fi

# Accept comma- and/or space-separated lists.
split_list() { echo "$1" | tr ',' ' '; }
read -r -a DRUG_ARR <<< "$(split_list "$DRUGS")"
read -r -a SOLVER_ARR <<< "$(split_list "$SOLVERS")"
read -r -a NP_ARR <<< "$(split_list "$NP_LIST")"

if [[ ${#DRUG_ARR[@]} -eq 0 || ${#SOLVER_ARR[@]} -eq 0 || ${#NP_ARR[@]} -eq 0 ]]; then
  echo "DRUGS, SOLVERS and NP_LIST must each resolve to at least one value." >&2
  exit 2
fi

mkdir -p "$OUT_DIR/logs"
CSV_FILE="$OUT_DIR/results.csv"
echo "timestamp_utc,variant_flag,drug,solver,np,repeat,pacing,pacing_write,time_step_min,concentrations,hill_file,herg_file,elapsed_seconds,exit_code,log_file,host_os,host_arch,host_cpu_model,host_cpu_cores,host_ram_gb,docker_runtime,container_arch,container_cpu_cores,container_ram_gb,emulated" > "$CSV_FILE"

echo "Output directory: $OUT_DIR"
echo "Drugs:   ${DRUG_ARR[*]}"
echo "Solvers: ${SOLVER_ARR[*]}"
echo "NP:      ${NP_ARR[*]}"
echo "Repeats: $REPEATS"
echo

if [[ "$DRY_RUN" != "1" ]]; then
  if ! docker compose ps --status running --services 2>/dev/null | grep -qx "drugsim-linux"; then
    echo "drugsim-linux container is not running, starting it (docker compose up -d)..."
    docker compose up -d
  fi

  gather_container_info

  if [[ "$SKIP_BUILD" != "1" ]]; then
    echo "Building binary for VARIANT_FLAG=$VARIANT_FLAG once before the sweep (RUN_SIM=0)..."
    docker compose exec -T drugsim-linux env \
      "VARIANT_FLAG=${VARIANT_FLAG}" "RUN_SIM=0" \
      bash docker/linux/build_and_run.sh
    build_status=$?
    if [[ $build_status -ne 0 ]]; then
      echo "Build failed (exit $build_status); aborting sweep." >&2
      exit $build_status
    fi
    echo
  fi
fi

echo "Hardware / environment:"
echo "  Host:      ${HOST_OS} ${HOST_ARCH} | ${HOST_CPU_MODEL} | ${HOST_CPU_CORES} logical cores | ${HOST_RAM_GB} GB RAM"
echo "  Container: ${DOCKER_RUNTIME} | arch=${CONTAINER_ARCH} (emulated=${EMULATED}) | ${CONTAINER_CPU_CORES} cores visible | ${CONTAINER_RAM_GB} GB RAM visible"
if [[ "$EMULATED" == "yes" ]]; then
  echo "  WARNING: container arch differs from host arch -- this run is CPU-emulated (e.g. Rosetta/QEMU)."
  echo "           Timings will NOT reflect native hardware performance; do not compare them against native runs."
fi
echo

{
  echo "host_os=${HOST_OS}"
  echo "host_arch=${HOST_ARCH}"
  echo "host_cpu_model=${HOST_CPU_MODEL}"
  echo "host_cpu_cores=${HOST_CPU_CORES}"
  echo "host_ram_gb=${HOST_RAM_GB}"
  echo "docker_runtime=${DOCKER_RUNTIME}"
  echo "container_arch=${CONTAINER_ARCH}"
  echo "container_cpu_cores=${CONTAINER_CPU_CORES}"
  echo "container_ram_gb=${CONTAINER_RAM_GB}"
  echo "emulated=${EMULATED}"
} > "$OUT_DIR/environment.txt"

total=$(( ${#DRUG_ARR[@]} * ${#SOLVER_ARR[@]} * ${#NP_ARR[@]} * REPEATS ))
count=0
failures=0

for drug in "${DRUG_ARR[@]}"; do
  hill_file="./chantest_hill/${drug}/IC50_samples${HILL_SAMPLES}.csv"
  herg_file="./chantest_herg/${drug}/boot_pars${HERG_SAMPLES}.csv"

  for solver in "${SOLVER_ARR[@]}"; do
    for np in "${NP_ARR[@]}"; do
      for ((rep = 1; rep <= REPEATS; rep++)); do
        count=$((count + 1))
        label="drug=${drug} solver=${solver} np=${np} rep=${rep}"
        printf '[%d/%d] %s\n' "$count" "$total" "$label"

        run_id="${drug}_${solver}_np${np}_rep${rep}"
        log_file="$OUT_DIR/logs/${run_id}.log"

        cmd=(docker compose exec -T drugsim-linux env
          "VARIANT_FLAG=${VARIANT_FLAG}"
          "NP=${np}"
          "PACING_OVERRIDE=${PACING}"
          "PACING_WRITE_OVERRIDE=${PACING_WRITE}"
          "TIME_STEP_MIN_OVERRIDE=${TIME_STEP_MIN}"
          "SOLVER_OVERRIDE=${solver}"
          "DRUG_NAME_OVERRIDE=${drug}"
          "HILL_FILE_OVERRIDE=${hill_file}"
          "HERG_FILE_OVERRIDE=${herg_file}")
        [[ -n "$CONC" ]] && cmd+=("CONC_OVERRIDE=${CONC}")
        cmd+=(bash docker/linux/run_only.sh)

        if [[ "$DRY_RUN" == "1" ]]; then
          printf '  DRY RUN: %s\n' "${cmd[*]}"
          continue
        fi

        timestamp="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        TIMEFORMAT='%R'
        time_file="$(mktemp)"
        { time "${cmd[@]}" >"$log_file" 2>&1; } 2>"$time_file"
        exit_code=$?
        elapsed="$(cat "$time_file")"
        rm -f "$time_file"

        conc_used="${CONC:-<param.txt default>}"
        echo "${timestamp},${VARIANT_FLAG},${drug},${solver},${np},${rep},${PACING},${PACING_WRITE},${TIME_STEP_MIN},${conc_used},${hill_file},${herg_file},${elapsed},${exit_code},${log_file},${HOST_OS},${HOST_ARCH},${HOST_CPU_MODEL},${HOST_CPU_CORES},${HOST_RAM_GB},${DOCKER_RUNTIME},${CONTAINER_ARCH},${CONTAINER_CPU_CORES},${CONTAINER_RAM_GB},${EMULATED}" >> "$CSV_FILE"

        if [[ $exit_code -ne 0 ]]; then
          failures=$((failures + 1))
          echo "  FAILED (exit $exit_code) after ${elapsed}s -- see $log_file"
          if [[ "$FAIL_FAST" == "1" ]]; then
            echo "FAIL_FAST=1, stopping sweep." >&2
            exit $exit_code
          fi
        else
          echo "  OK in ${elapsed}s"
        fi
      done
    done
  done
done

if [[ "$DRY_RUN" == "1" ]]; then
  echo "Dry run complete: $total combinations would run."
  exit 0
fi

echo
echo "Sweep complete: $((count - failures))/$total runs succeeded ($failures failed)."
echo "Raw results: $CSV_FILE"
echo "Environment: $OUT_DIR/environment.txt"

summary_file="$OUT_DIR/summary.txt"
awk -F',' '
  NR == 1 { next }
  {
    key = $3 "|" $4 "|" $5
    if ($14 == "0") {
      n[key]++
      sum[key] += $13
      sumsq[key] += $13 * $13
    } else {
      fail[key]++
    }
    seen[key] = 1
  }
  END {
    printf "%-22s %-8s %-4s %-6s %-12s %-12s %-6s\n", "drug", "solver", "np", "n_ok", "mean_s", "stdev_s", "fail"
    for (k in seen) {
      split(k, a, "|")
      cnt = n[k] + 0
      mean = (cnt > 0) ? sum[k] / cnt : 0
      var = (cnt > 1) ? (sumsq[k] - cnt * mean * mean) / (cnt - 1) : 0
      sd = (var > 0) ? sqrt(var) : 0
      printf "%-22s %-8s %-4s %-6d %-12.3f %-12.3f %-6d\n", a[1], a[2], a[3], cnt, mean, sd, (fail[k] + 0)
    }
  }
' "$CSV_FILE" | sort | tee "$summary_file"

echo
echo "Summary: $summary_file"
