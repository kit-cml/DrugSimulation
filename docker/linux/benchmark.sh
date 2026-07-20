#!/usr/bin/env bash
# Hardware-performance benchmark harness for DrugSimulation.
#
# Sweeps the drugsim container across drugs x solvers x MPI process counts,
# holding every other hyperparameter fixed, and records wall-clock elapsed
# time for each run to a CSV so runs are comparable across hardware.
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
echo "timestamp_utc,variant_flag,drug,solver,np,repeat,pacing,pacing_write,time_step_min,concentrations,hill_file,herg_file,elapsed_seconds,exit_code,log_file" > "$CSV_FILE"

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
        echo "${timestamp},${VARIANT_FLAG},${drug},${solver},${np},${rep},${PACING},${PACING_WRITE},${TIME_STEP_MIN},${conc_used},${hill_file},${herg_file},${elapsed},${exit_code},${log_file}" >> "$CSV_FILE"

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
