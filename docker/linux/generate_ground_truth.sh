#!/usr/bin/env bash
# Ground-truth generator for DrugSimulation.
#
# Unlike benchmark.sh (which shrinks the workload to measure hardware speed),
# this runs the FULL, real simulation -- real number_pacing, full hill/herg
# population, real per-drug concentrations -- using the slow-but-trusted
# Euler solver, for every drug. It gets done faster not by cutting the
# workload down, but by running several drugs' simulations CONCURRENTLY
# (parallel across drugs) instead of one drug at a time, since MPI ranks
# within a single drug run (NP) hit diminishing returns long before the
# container's CPU budget is used up.
#
# Run from the repository root (or anywhere; it cd's to the repo root itself):
#   bash docker/linux/generate_ground_truth.sh
set -u -o pipefail

usage() {
  cat <<'EOF'
Usage: [ENV=value ...] bash docker/linux/generate_ground_truth.sh

Runs the full (real number_pacing, full hill/herg population, real per-drug
concentrations) Euler simulation for every drug, PARALLEL_DRUGS at a time
concurrently, so the whole batch finishes faster than running drugs one by
one. Each drug still gets NP MPI ranks internally, same as any other run.

Environment variables (defaults shown):
  DRUGS              (auto-discovered from chantest_hill/*)  comma/space separated
  VARIANT_FLAG        -DCIPAORDV1                             cell-model build variant
  SOLVER               Euler                                   solver_type override (ground truth = trusted Euler)
  TIME_STEP_MIN         0.005                                   fixed Euler step (must stay small for stability)
  PACING              (unset -> param.txt's own number_pacing, e.g. 1000: a REAL run, not a smoke test)
  PACING_WRITE        (unset -> param.txt's own number_pacing_write)
  NP                    2                                       MPI ranks PER DRUG
  PARALLEL_DRUGS      (auto: container cores / NP, min 1)      how many drugs run concurrently
  CONC_LEVELS           all                                     which CMAX_TABLE columns per drug (default: full 1x-4x dose range)
  CMAX_TABLE           <SIM_DIR>/chantest_hill/CiPA_drug_concentrations.csv   per-drug Cmax table
  SAMPLES               all                                     hill/herg rows per drug; "all" = full HILL_SAMPLES/HERG_SAMPLES file
  HILL_SAMPLES          30                                      -> chantest_hill/<drug>/IC50_samples<HILL_SAMPLES>.csv
  HERG_SAMPLES          10                                      -> chantest_herg/<drug>/boot_pars<HERG_SAMPLES>.csv
  SIM_DIR              bin/simulation_cpu
  SKIP_BUILD             0                                       set 1 to reuse an already-built binary
  FAIL_FAST              0                                       set 1 to stop launching further batches after any failure
  OUT_DIR              ground_truth_results/<timestamp>         where results.csv/logs/ are written
  DRY_RUN                0                                       set 1 to print the per-drug commands without running them

Real simulation output (time-series/feature CSVs) still lands in the shared
bin/simulation_cpu/results/<conc>/ tree, same as any other run -- those
filenames are drug-name-prefixed so concurrent drugs don't clobber each
other's data. The one shared, non-prefixed file is results/logfile.docker,
which concurrent runs will interleave; ignore it and use OUT_DIR/logs/ instead
(this script captures each drug's own output independently via its own
docker compose exec stream, unaffected by that in-container file).

Example - defaults (all drugs, real pacing, Euler, full population):
  bash docker/linux/generate_ground_truth.sh

Example - bump per-drug parallelism to 4 ranks, 3 drugs at a time:
  NP=4 PARALLEL_DRUGS=3 bash docker/linux/generate_ground_truth.sh

Example - quick plumbing smoke test (NOT real ground truth):
  PACING=5 PACING_WRITE=5 DRUGS="quinidine,bepridil" bash docker/linux/generate_ground_truth.sh
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
SOLVER="${SOLVER:-Euler}"
TIME_STEP_MIN="${TIME_STEP_MIN:-0.005}"
PACING="${PACING:-}"
PACING_WRITE="${PACING_WRITE:-}"
NP="${NP:-2}"
PARALLEL_DRUGS="${PARALLEL_DRUGS:-}"
CONC_LEVELS="${CONC_LEVELS:-all}"
CMAX_TABLE="${CMAX_TABLE:-$SIM_DIR/chantest_hill/CiPA_drug_concentrations.csv}"
SAMPLES="${SAMPLES:-all}"
HILL_SAMPLES="${HILL_SAMPLES:-30}"
HERG_SAMPLES="${HERG_SAMPLES:-10}"
SKIP_BUILD="${SKIP_BUILD:-0}"
FAIL_FAST="${FAIL_FAST:-0}"
DRY_RUN="${DRY_RUN:-0}"
OUT_DIR="${OUT_DIR:-ground_truth_results/$(date +%Y%m%d_%H%M%S)}"

# Accept comma- and/or space-separated lists.
split_list() { echo "$1" | tr ',' ' '; }

if [[ -z "${DRUGS:-}" ]]; then
  DRUGS="$(find "$SIM_DIR/chantest_hill" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort | tr '\n' ',')"
fi
read -r -a DRUG_ARR <<< "$(split_list "$DRUGS")"
if [[ ${#DRUG_ARR[@]} -eq 0 ]]; then
  echo "DRUGS resolved to an empty list." >&2
  exit 2
fi

# --- Per-drug concentration lookup (same convention as benchmark.sh) -------
resolve_conc_levels() {
  if [[ "$CONC_LEVELS" == "all" ]]; then
    echo "1x 2x 3x 4x"
  else
    split_list "$CONC_LEVELS"
  fi
}

lookup_drug_conc() {
  local drug="$1"
  local levels
  levels="$(resolve_conc_levels)"
  [[ -f "$CMAX_TABLE" ]] || return 1
  awk -F',' -v drug="$drug" -v levels="$levels" '
    NR == 1 {
      for (i = 1; i <= NF; i++) col[$i] = i
      next
    }
    $1 == drug {
      n = split(levels, lv, " ")
      out = ""
      ok = 1
      for (j = 1; j <= n; j++) {
        idx = col["conc_" lv[j]]
        val = (idx != "") ? $(idx) : ""
        if (val == "") { ok = 0; break }
        out = (out == "") ? val : out "," val
      }
      if (ok) { print out; found = 1 }
    }
    END { exit (found ? 0 : 1) }
  ' "$CMAX_TABLE"
}

mkdir -p "$OUT_DIR/logs" "$OUT_DIR/.results"
CSV_FILE="$OUT_DIR/results.csv"
echo "drug,solver,np,concentrations,hill_file,herg_file,samples,elapsed_seconds,exit_code,log_file" > "$CSV_FILE"

SAMPLES_TMP_DIR="$SIM_DIR/.bench_samples"
mkdir -p "$SAMPLES_TMP_DIR"

if [[ "$DRY_RUN" != "1" ]]; then
  if ! docker compose ps --status running --services 2>/dev/null | grep -qx "drugsim-linux"; then
    echo "drugsim-linux container is not running, starting it (docker compose up -d)..."
    docker compose up -d
  fi

  if [[ -z "$PARALLEL_DRUGS" ]]; then
    container_cores="$(docker compose exec -T drugsim-linux nproc 2>/dev/null | tr -d '\r\n')"
    if [[ -n "$container_cores" ]]; then
      PARALLEL_DRUGS=$((container_cores / NP))
      ((PARALLEL_DRUGS < 1)) && PARALLEL_DRUGS=1
      echo "PARALLEL_DRUGS not set; auto-detected ${container_cores} container cores / NP=${NP} -> running ${PARALLEL_DRUGS} drug(s) at a time."
    else
      PARALLEL_DRUGS=1
      echo "Could not detect container core count; defaulting PARALLEL_DRUGS=1."
    fi
  fi

  if [[ "$SKIP_BUILD" != "1" ]]; then
    echo "Building binary for VARIANT_FLAG=$VARIANT_FLAG once before the batch (RUN_SIM=0)..."
    docker compose exec -T drugsim-linux env \
      "VARIANT_FLAG=${VARIANT_FLAG}" "RUN_SIM=0" \
      bash docker/linux/build_and_run.sh
    build_status=$?
    if [[ $build_status -ne 0 ]]; then
      echo "Build failed (exit $build_status); aborting." >&2
      exit $build_status
    fi
    echo
  fi
else
  [[ -z "$PARALLEL_DRUGS" ]] && PARALLEL_DRUGS=1
fi

echo "Output directory: $OUT_DIR"
echo "Drugs:          ${DRUG_ARR[*]}"
echo "Solver:         ${SOLVER} (time_step_min=${TIME_STEP_MIN})"
echo "Pacing:         $( [[ -n "$PACING" ]] && echo "override=${PACING}/${PACING_WRITE:-$PACING}" || echo "param.txt default (real run)" )"
echo "NP per drug:    $NP"
echo "Parallel drugs: $PARALLEL_DRUGS"
echo "Conc levels:    ${CONC_LEVELS} (per-drug, from ${CMAX_TABLE})"
echo "Samples:        ${SAMPLES} (from HILL_SAMPLES=${HILL_SAMPLES}/HERG_SAMPLES=${HERG_SAMPLES})"
echo

# Runs one drug's full simulation; meant to be launched with `&`. Writes its
# outcome to a per-drug result file since background jobs can't hand values
# back to the parent shell directly.
run_drug() {
  local drug="$1"
  local hill_file herg_file src_hill src_herg avail_hill avail_herg samples_used
  local conc_override cmd log_file result_file timestamp time_file elapsed exit_code

  hill_file="./chantest_hill/${drug}/IC50_samples${HILL_SAMPLES}.csv"
  herg_file="./chantest_herg/${drug}/boot_pars${HERG_SAMPLES}.csv"
  samples_used="all"

  if [[ "$SAMPLES" != "all" && "$SAMPLES" != "0" ]]; then
    src_hill="$SIM_DIR/chantest_hill/${drug}/IC50_samples${HILL_SAMPLES}.csv"
    src_herg="$SIM_DIR/chantest_herg/${drug}/boot_pars${HERG_SAMPLES}.csv"
    avail_hill=$(($(wc -l < "$src_hill") - 1))
    avail_herg=$(($(wc -l < "$src_herg") - 1))
    samples_used=$((avail_hill < avail_herg ? avail_hill : avail_herg))
    if ((SAMPLES < samples_used)); then
      samples_used=$SAMPLES
    elif ((SAMPLES > samples_used)); then
      echo "WARNING: SAMPLES=${SAMPLES} exceeds available rows for ${drug} (hill=${avail_hill}, herg=${avail_herg}); using ${samples_used} instead." >&2
    fi
    hill_file="./.bench_samples/${drug}_hill_n${samples_used}.csv"
    herg_file="./.bench_samples/${drug}_herg_n${samples_used}.csv"
    head -n "$((samples_used + 1))" "$src_hill" > "$SAMPLES_TMP_DIR/${drug}_hill_n${samples_used}.csv"
    head -n "$((samples_used + 1))" "$src_herg" > "$SAMPLES_TMP_DIR/${drug}_herg_n${samples_used}.csv"
  fi

  if ! conc_override="$(lookup_drug_conc "$drug")"; then
    conc_override=""
    echo "WARNING: no concentrations for drug=${drug} (levels=${CONC_LEVELS}) in ${CMAX_TABLE}; falling back to param.txt's own concentrations for this drug." >&2
  fi

  cmd=(docker compose exec -T drugsim-linux env
    "VARIANT_FLAG=${VARIANT_FLAG}"
    "NP=${NP}"
    "SOLVER_OVERRIDE=${SOLVER}"
    "TIME_STEP_MIN_OVERRIDE=${TIME_STEP_MIN}"
    "DRUG_NAME_OVERRIDE=${drug}"
    "HILL_FILE_OVERRIDE=${hill_file}"
    "HERG_FILE_OVERRIDE=${herg_file}")
  [[ -n "$PACING" ]] && cmd+=("PACING_OVERRIDE=${PACING}")
  [[ -n "$PACING_WRITE" ]] && cmd+=("PACING_WRITE_OVERRIDE=${PACING_WRITE}")
  [[ -n "$conc_override" ]] && cmd+=("CONC_OVERRIDE=${conc_override}")
  cmd+=(bash docker/linux/run_only.sh)

  log_file="$OUT_DIR/logs/${drug}.log"
  result_file="$OUT_DIR/.results/${drug}.result"

  if [[ "$DRY_RUN" == "1" ]]; then
    printf '  DRY RUN [%s]: %s\n' "$drug" "${cmd[*]}"
    local dry_conc_csv="${conc_override:-<param.txt default>}"
    dry_conc_csv="${dry_conc_csv//,/;}"
    echo "${drug},${SOLVER},${NP},${dry_conc_csv},${hill_file},${herg_file},${samples_used},0,0,${log_file}" > "$result_file"
    return 0
  fi

  echo "  -> starting ${drug} (NP=${NP})"
  TIMEFORMAT='%R'
  time_file="$(mktemp)"
  { time "${cmd[@]}" >"$log_file" 2>&1; } 2>"$time_file"
  exit_code=$?
  elapsed="$(cat "$time_file")"
  rm -f "$time_file"

  if [[ $exit_code -ne 0 ]]; then
    echo "  <- FAILED ${drug} (exit $exit_code) after ${elapsed}s -- see $log_file"
  else
    echo "  <- done ${drug} in ${elapsed}s"
  fi

  local conc_csv="${conc_override:-<param.txt default>}"
  conc_csv="${conc_csv//,/;}"
  echo "${drug},${SOLVER},${NP},${conc_csv},${hill_file},${herg_file},${samples_used},${elapsed},${exit_code},${log_file}" > "$result_file"
}

batch_start_epoch="$(date +%s)"
total=${#DRUG_ARR[@]}
i=0
overall_failures=0

while ((i < total)); do
  batch=("${DRUG_ARR[@]:i:PARALLEL_DRUGS}")
  echo "Batch: ${batch[*]}"

  pids=()
  for drug in "${batch[@]}"; do
    run_drug "$drug" &
    pids+=($!)
  done
  for pid in "${pids[@]}"; do
    wait "$pid"
  done
  echo

  batch_failures=0
  for drug in "${batch[@]}"; do
    result_file="$OUT_DIR/.results/${drug}.result"
    if [[ -f "$result_file" ]]; then
      exit_code_field="$(cut -d',' -f9 "$result_file")"
      [[ "$exit_code_field" != "0" ]] && batch_failures=$((batch_failures + 1))
    else
      echo "WARNING: no result recorded for ${drug} (job may have crashed the shell)." >&2
      batch_failures=$((batch_failures + 1))
    fi
  done
  overall_failures=$((overall_failures + batch_failures))

  i=$((i + PARALLEL_DRUGS))

  if ((batch_failures > 0 && FAIL_FAST == 1)) && ((i < total)); then
    echo "FAIL_FAST=1 and a run failed in this batch; not launching the remaining $((total - i)) drug(s)." >&2
    break
  fi
done

batch_end_epoch="$(date +%s)"
total_wall_seconds=$((batch_end_epoch - batch_start_epoch))

# Aggregate whatever results were actually produced, in DRUGS order.
sum_elapsed=0
completed=0
for drug in "${DRUG_ARR[@]}"; do
  result_file="$OUT_DIR/.results/${drug}.result"
  [[ -f "$result_file" ]] || continue
  cat "$result_file" >> "$CSV_FILE"
  elapsed_field="$(cut -d',' -f8 "$result_file")"
  sum_elapsed="$(awk -v a="$sum_elapsed" -v b="$elapsed_field" 'BEGIN { print a + b }')"
  completed=$((completed + 1))
done
rm -rf "$OUT_DIR/.results"

if [[ "$DRY_RUN" == "1" ]]; then
  echo "Dry run complete: $total drug(s) would run, $PARALLEL_DRUGS at a time."
  exit 0
fi

echo "Batch complete: $((completed - overall_failures))/$total drug(s) succeeded ($overall_failures failed)."
echo "Wall-clock time: ${total_wall_seconds}s (sum of per-drug times: ${sum_elapsed}s -- speedup from running ${PARALLEL_DRUGS} at a time)"
echo "Raw results:     $CSV_FILE"
echo "Logs:             $OUT_DIR/logs/"
echo "Simulation output (time-series/feature CSVs) is in ${SIM_DIR}/results/<concentration>/, same as any other run."

if ((overall_failures > 0)); then
  exit 1
fi
