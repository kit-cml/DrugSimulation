# DrugSimulation

A high-performance, MPI-parallel in-silico cardiac drug safety simulation framework based on the [CiPA (Comprehensive in vitro Proarrhythmia Assay)](https://cipaproject.org/) initiative. The tool simulates the electrophysiological effects of drugs on human ventricular cardiomyocytes across multiple cell models and drug concentrations, then extracts quantitative biomarkers (CiPA features) for cardiotoxicity assessment.

---

## Table of Contents

- [Overview](#overview)
- [Cell Models](#cell-models)
- [Prerequisites](#prerequisites)
- [Building](#building)
- [Project Structure](#project-structure)
- [Configuration (`param.txt`)](#configuration-paramtxt)
- [Input Files](#input-files)
- [Running a Simulation Locally on Your Machine](#running-a-simulation-locally-on-your-machine)
- [Running in Dockerized Linux (Recommended for macOS)](#running-in-dockerized-linux-recommended-for-macos)
- [Running a Simulation on an HPC Cluster (PBS/Torque)](#running-a-simulation-on-an-hpc-cluster-pbstorque)
- [Post-processing & Reporting](#post-processing--reporting)
- [Output Files](#output-files)
- [Python Scripts](#python-scripts)
- [AI Report Generation](#ai-report-generation)

---

## Overview

DrugSimulation runs a multi-pace action-potential (AP) simulation for each drug sample and concentration. It computes steady-state electrical behaviour over up to 1000 paces, then extracts the following CiPA biomarker features from the final paces:

| Feature | Description |
|---|---|
| `qNet` | Net charge (area under net current) |
| `qInward` | Inward charge |
| `INaLAUC` | Area under the late sodium current |
| `ICaLAUC` | Area under the L-type calcium current |
| `APD90` / `APD50` | Action potential duration at 90% / 50% repolarisation |
| `APDTri` | AP triangulation (APD90 − APD50) |
| `VmMax` / `VmRest` | Peak and resting membrane potential |
| `dVmdtMax` | Maximum upstroke velocity |
| `dVmdtMaxRepol` | Maximum repolarisation rate |
| `CaTD90` / `CaTD50` | Calcium transient duration at 90% / 50% |
| `CaTDTri` | Calcium transient triangulation |
| `CaMax` / `CaRest` | Peak and resting intracellular calcium |

MPI is used to distribute samples across CPU cores. The simulation supports both a full in-silico run and a faster postprocessing-only mode that reads pre-computed steady-state values.

---

## Cell Models

Select the cardiac cell model at compile time using one of the preprocessor macros below:

| Macro | Binary suffix | Cell model |
|---|---|---|
| *(none)* | `ORd-static` | O'Hara–Rudy 2011 (default) |
| `-DCIPAORDV1` | `CiPAORdv1.0` | O'Hara–Rudy CiPA v1.0 (2017) |
| `-DTOR_ORD` | `ToR-ORd` | Tomek–O'Hara–Rudy (ToR-ORd) |
| `-DTOR_ORD_DYNCL` | `ToR-ORd-dynCl` | ToR-ORd with dynamic chloride |
| `-DGRANDI` | `Grandi` | Grandi 2011 atrial |
| `-DORD_STATIC_LAND` | `ORd-static_Land` | O'Hara–Rudy + Land contraction |
| `-DCIPAORDV1_LAND` | `CiPAORdv1.0_Land` | CiPAORdv1.0 + Land contraction |
| `-DTOR_ORD_LAND` | `ToR-ORd_Land` | ToR-ORd + Land contraction |

Add `-DPOSTPROCESSING` to any of the above to build the postprocessing variant.

---

## Prerequisites

| Dependency | Notes |
|---|---|
| C++11 compiler with MPI (`mpicxx` / `mpicc`) | e.g. OpenMPI or MPICH |
| [libCML](../libCML/) | Cardiac Modelling Library (expected at `../libCML/`) |
| [SUNDIALS CVode 5.7.0](https://computing.llnl.gov/projects/sundials) | ODE solver (expected at `../libCML/libs/sundials-5.7.0/`) |
| `libcurl` | Used by the report generator |
| `libjson-c` | Used by the report generator |
| Python ≥ 3.6 | For post-processing scripts |
| Python packages | See [`bin/scripts/requirements.txt`](bin/scripts/requirements.txt) |

Install Python dependencies:

```bash
pip install -r bin/scripts/requirements.txt
```

---

## Building

Native builds are intended for Linux environments that have matching `libCML` and SUNDIALS static archives. On macOS, use the Dockerized Linux workflow below because the bundled archives under `../libCML/libs/sundials-5.7.0/lib64` are Linux target files and will not link natively on macOS.

### Build all variants at once

```bash
bash compile_binaries.sh
```

This builds all eight cell-model variants plus their postprocessing counterparts and places the binaries in `bin/`.

### Build a single variant manually

```bash
# Example: CiPAORdv1.0
make all EXTRA_CXXFLAGS+="-std=c++11 -DCIPAORDV1"

# Example: ToR-ORd postprocessing
make all EXTRA_CXXFLAGS+="-std=c++11 -DTOR_ORD -DPOSTPROCESSING"

# Clean build artefacts for a variant
make EXTRA_CXXFLAGS+="-DCIPAORDV1" clean
```

Compiled binaries are written to `bin/drugsim_<variant>`.

---

## Project Structure

```
DrugSimulation/
├── main.cpp                        # Entry point – MPI init, parameter loading, simulation dispatch
├── Makefile                        # Build system
├── compile_binaries.sh             # Builds all cell-model variants in one step
├── modules/
│   ├── drug_simulation_bench.cpp/hpp   # Top-level simulation orchestrator (MPI work distribution)
│   ├── insilico.cpp/hpp                # Single-sample in-silico AP simulation loop
│   ├── postprocessing.cpp/hpp          # Feature extraction from pre-computed steady states
│   ├── report_drug.cpp/hpp             # Invokes the AI report-generation pipeline
│   └── show_param_logs.cpp/hpp         # Prints parameter summary to stdout
└── bin/
    ├── drugsim_<variant>           # Compiled binaries (git-ignored)
    ├── scripts/                    # Shared utility scripts
    │   ├── combine_and_average_features.py   # Merges per-core feature CSVs; computes per-concentration averages
    │   ├── compute_mean_ci.py                # Bootstrap median + 95% CI across concentrations
    │   ├── plot_features.py                  # Dose–response feature plots
    │   ├── plot_time_series.py               # Time-series plots (Vm, dVm/dt, Cai, currents)
    │   ├── convert_params.sh                 # Converts param.txt keys to camelCase
    │   ├── create_concs_directories.sh       # Creates result sub-directories per concentration
    │   ├── create_escaped_b64_file.sh        # Escapes & Base64-encodes files for API upload
    │   ├── unzip_files.sh / zip_files.sh     # Archive helpers for result files
    │   └── requirements.txt                  # Python package pinned versions
    ├── simulation_cpu/             # Ready-to-run working directory for a standard CPU simulation
    │   ├── param.txt               # Simulation parameters (edit this before running)
    │   ├── exec_bash.sh            # Launch simulation locally via mpiexec
    │   ├── exec_pbs.sh             # Launch simulation on a PBS/Torque HPC cluster
    │   ├── kill_bash.sh            # Kill a running local simulation
    │   ├── generate_report.sh      # Run post-processing and generate figures
    │   ├── exec_perplexity_report_generator.sh  # AI narrative report via Perplexity API
    │   ├── chantest_hill/          # IC50 / Hill coefficient input files per drug
    │   ├── chantest_herg/          # hERG-specific parameter files per drug
    │   ├── population/             # Conductance variability (virtual population) files
    │   ├── initial_states/         # Pre-computed steady-state initial conditions
    │   ├── kitox_error_data/       # Reference IC50 data
    │   └── requests/               # Perplexity API prompt templates
    └── simulation_cpu_postprocessing/  # Working directory for postprocessing-only runs
        └── (same layout as simulation_cpu/)
```

---

## Configuration (`param.txt`)

Edit `bin/simulation_cpu/param.txt` before running. Key parameters:

| Parameter | Example | Description |
|---|---|---|
| `user_name` | `marcell` | Label used in output filenames |
| `number_pacing` | `1000` | Total number of paces to simulate |
| `number_pacing_write` | `10` | Number of final paces written to disk |
| `cycle_length` | `2000` | Pacing cycle length (ms) |
| `cell_model` | `CiPAORdv1.0_endo` | Cell model + cell type (`endo`, `epi`, `myo`) |
| `stimulus_duration` | `0.5` | Stimulus current duration (ms) |
| `stimulus_amplitude_scale` | `1.0` | Scaling factor for stimulus amplitude |
| `solver_type` | `CVode` | ODE solver: `CVode` or `Euler` |
| `time_step_min` | `0.5` | Minimum adaptive time step for CVode (ms) |
| `time_step_max` | `2.0` | Maximum adaptive time step for CVode (ms) |
| `writing_step` | `1.0` | Output writing interval (ms) |
| `drug_name` | `quinidine` | Drug name (used in filenames) |
| `drug_concentrations` | `3237,6474,9711,12948` | Comma-separated concentrations to simulate (µM) |
| `hill_file` | `./chantest_hill/quinidine/IC50_samples30.csv` | IC50 & Hill coefficient samples |
| `herg_file` | `./chantest_herg/quinidine/boot_pars10.csv` | hERG bootstrapped parameters |
| `is_cvar` | `0` | Enable conductance variability (`1` = on) |
| `cvar_file` | `./population/indi_1sample.csv` | Conductance variability population file |
| `number_of_cpu` | `10` | Number of MPI processes for `mpiexec` |

---

## Input Files

### IC50 / Hill coefficients (`hill_file`)

CSV file where each row is one virtual sample. Columns define per-channel IC50 and Hill coefficient pairs:

```
ICaL_IC50, ICaL_h, IK1_IC50, IK1_h, IKs_IC50, IKs_h,
INa_IC50, INa_h, INaL_IC50, INaL_h, Ito_IC50, Ito_h,
hERG_IC50, hERG_h
```

### hERG parameters (`herg_file`)

CSV file with bootstrapped hERG kinetic parameters for the CiPAORdv1.0 model.

### Conductance variability (`cvar_file`)

CSV file providing per-sample scaling factors for ionic conductances (e.g., `GNa_fast`, `GKr`, `PCa`, …). Enables virtual population simulations.

---

## Running a Simulation Locally on Your Machine

This is a complete end-to-end walkthrough for running a simulation on your local machine.

### Step 1 — Build the binaries

From the repository root, build all variants:

```bash
bash compile_binaries.sh
```

Or build only the variant you need (example: CiPAORdv1.0):

```bash
make all EXTRA_CXXFLAGS+="-std=c++11 -DCIPAORDV1"
```

Binaries are placed in `bin/`.

### Step 2 — Install Python dependencies

```bash
pip install -r bin/scripts/requirements.txt
```

### Step 3 — Configure the simulation

```bash
cd bin/simulation_cpu
```

Open `param.txt` and set the parameters for your run. The most important ones:

```ini
user_name          = your_name
cell_model         = CiPAORdv1.0_endo   # choose model + cell type (endo/epi/myo)
drug_name          = quinidine
drug_concentrations = 3237,6474,9711,12948
hill_file          = ./chantest_hill/quinidine/IC50_samples30.csv
herg_file          = ./chantest_herg/quinidine/boot_pars10.csv
number_pacing      = 1000
number_of_cpu      = 4               # set to the number of cores you want to use
```

> **Tip:** Set `number_of_cpu` to the number of physical/logical cores available on your machine. Run `nproc` (Linux) or `sysctl -n hw.logicalcpu` (macOS) to check.

### Step 4 — Launch the simulation

```bash
bash exec_bash.sh
```

`exec_bash.sh` will:
1. Read `param.txt` to determine the cell model, drug, concentrations, and CPU count.
2. Verify that `number_of_cpu` does not exceed the machine's available cores.
3. Create result sub-directories under `./results/<concentration>/`.
4. Run `mpiexec -np <number_of_cpu> ../drugsim_<variant> -input_deck param.txt`, logging all output to `./results/logfile`.
5. Zip time-series, feature, and initial-values CSVs per concentration.
6. Automatically call `generate_report.sh` to produce plots and summary statistics.

Progress is printed to the terminal; detailed logs are in `./results/logfile`.

### Step 5 — Monitor and stop

Tail the log in a separate terminal:

```bash
tail -f ./results/logfile
```

To kill a running simulation cleanly:

```bash
bash kill_bash.sh
```

### Step 6 — View results

After completion, results are in `./results/<concentration>/`:

```
results/
├── 0.00/                          # control (no drug)
├── 3237.00/                       # Cmax1
│   ├── quinidine_3237.00_time_series_smp0.csv
│   ├── quinidine_3237.00_features_core0.csv
│   └── …
├── quinidine-median-ci.csv        # bootstrap median + 95% CI across concentrations
└── logfile
```

Open the generated plots (`.png` / `.pdf`) in the concentration directories, or run `generate_report.sh` again separately if you want to regenerate figures without re-running the full simulation.

---

## Running in Dockerized Linux (Recommended for macOS)

This repository includes a Linux container workflow that preserves the expected original layout where `libCML` is a sibling directory.

The Docker service is pinned to `linux/amd64` so it matches the bundled SUNDIALS archives under `../libCML/libs/sundials-5.7.0/lib64`.

Expected host layout:

```text
/Users/iganarendra/
   DrugSimulation/
   libCML/
```

### 1. Build the container image

From the repository root:

```bash
docker compose up -d --build
```

### 2. Build and run simulation in one command (inside container)

```bash
docker compose exec -T drugsim-linux env VARIANT_FLAG=-DCIPAORDV1 NP=10 RUN_SIM=1 \
   bash docker/linux/build_and_run.sh
```

What this does:
1. Builds `libCML/build/libcml.a` inside Linux.
2. Builds DrugSimulation binary using bundled Linux SUNDIALS libs from `../libCML/libs/sundials-5.7.0/lib64`.
3. Runs `mpiexec` on `bin/drugsim_CiPAORdv1.0` with `param.txt` from `bin/simulation_cpu/`.
4. Writes runtime logs to `bin/simulation_cpu/results/logfile.docker`.

### 3. Build only (skip simulation)

```bash
docker compose exec -T drugsim-linux env RUN_SIM=0 VARIANT_FLAG=-DCIPAORDV1 \
   bash docker/linux/build_and_run.sh
```

### 3a. Run a short pacing test (smoke test) without editing `param.txt`

```bash
docker compose exec -T drugsim-linux env VARIANT_FLAG=-DCIPAORDV1 NP=1 PACING_OVERRIDE=5 PACING_WRITE_OVERRIDE=5 TIME_STEP_MIN_OVERRIDE=0.001 SOLVER_OVERRIDE=Euler \
   RUN_SIM=1 bash docker/linux/build_and_run.sh
```

This creates a temporary parameter file inside the container run, leaves your checked-in `bin/simulation_cpu/param.txt` unchanged, and runs the full simulation flow with only 5 paces using the Euler solver.

### 4. Supported variant flags

- `-DCIPAORDV1`
- `-DTOR_ORD`
- `-DTOR_ORD_DYNCL`
- `-DGRANDI`
- `-DORD_STATIC_LAND`
- `-DCIPAORDV1_LAND`
- `-DTOR_ORD_LAND`
- empty string (`""`) for ORd-static default

### 5. Rollback / cleanup

Container-only cleanup (safe):

```bash
docker compose down --remove-orphans
docker image prune -f
```

Reset build artifacts in both repos:

```bash
cd /Users/iganarendra/DrugSimulation && make clean || true
cd /Users/iganarendra/libCML && make clean || true
```

Revert dockerization file changes in this repo:

```bash
cd /Users/iganarendra/DrugSimulation
git restore --source=HEAD --staged --worktree -- Makefile README.md docker-compose.yml docker/
```

---

## NEW: Run a Hardware Performance Benchmark

Use `docker/linux/benchmark.sh` to run a reproducible timing sweep inside the Linux container.

The benchmark auto-discovers drugs from `bin/simulation_cpu/chantest_hill/`, then sweeps:

- drug
- solver (`CVode` and `Euler` by default)
- MPI process count (`NP_LIST`)
- repeat index (`REPEATS`)

Each run calls `docker/linux/run_only.sh`, times the full simulation wall-clock duration, and appends one row to a benchmark CSV.

### Default benchmark profile

To keep the benchmark fast and comparable, the script uses short simulation overrides unless you change them:

- `PACING=5`
- `PACING_WRITE=5`
- `TIME_STEP_MIN=0.001`
- `SOLVERS="CVode,Euler"`
- `NP_LIST="1"`
- `REPEATS=1`
- `VARIANT_FLAG=-DCIPAORDV1`

It also builds the selected variant once before the sweep unless `SKIP_BUILD=1` is set.

### Basic usage

From the repository root, ready to run, just copy and paste to your terminal:

```bash
bash docker/linux/benchmark.sh
```

**The script will start `drugsim-linux` automatically if needed!**

### Common examples

Single-process baseline:

```bash
bash docker/linux/benchmark.sh
```

Compare 1, 2, 4, and 8 MPI ranks:

```bash
NP_LIST="1,2,4,8" REPEATS=3 bash docker/linux/benchmark.sh
```

Benchmark only selected drugs:

```bash
DRUGS="quinidine,bepridil" NP_LIST="1,4" REPEATS=2 bash docker/linux/benchmark.sh
```

Benchmark only one solver:

```bash
SOLVERS="Euler" NP_LIST="1,4" bash docker/linux/benchmark.sh
```

Reuse an already-built binary:

```bash
SKIP_BUILD=1 NP_LIST="1,4" bash docker/linux/benchmark.sh
```

Preview the sweep without running it:

```bash
DRY_RUN=1 NP_LIST="1,4" bash docker/linux/benchmark.sh
```

### Useful environment variables

| Variable | Default | Meaning |
|---|---|---|
| `DRUGS` | auto-discovered | Comma- or space-separated drug list |
| `SOLVERS` | `CVode,Euler` | Solvers to benchmark |
| `NP_LIST` | `1` | MPI process counts to test |
| `REPEATS` | `1` | Repeats per combination |
| `VARIANT_FLAG` | `-DCIPAORDV1` | Cell-model variant built once for the sweep |
| `PACING` | `5` | Override for `number_pacing` |
| `PACING_WRITE` | `5` | Override for `number_pacing_write` |
| `TIME_STEP_MIN` | `0.001` | Override for `time_step_min` |
| `CONC` | unset | Override `drug_concentrations` for every run |
| `HILL_SAMPLES` | `30` | Uses `chantest_hill/<drug>/IC50_samples<HILL_SAMPLES>.csv` |
| `HERG_SAMPLES` | `10` | Uses `chantest_herg/<drug>/boot_pars<HERG_SAMPLES>.csv` |
| `SKIP_BUILD` | `0` | Set to `1` to skip the one-time build step |
| `FAIL_FAST` | `0` | Set to `1` to stop on the first failed run |
| `OUT_DIR` | `benchmark_results/<timestamp>` | Output directory for this sweep |
| `DRY_RUN` | `0` | Set to `1` to print commands only |

Run `bash docker/linux/benchmark.sh --help` to print the full built-in usage text.

### Output files

Each sweep writes a new directory under `benchmark_results/<timestamp>/` containing:

- `results.csv` — one row per run with elapsed time, exit code, selected parameters, and hardware metadata
- `summary.txt` — mean and standard deviation of elapsed time grouped by `drug`, `solver`, and `np`
- `environment.txt` — host and container details captured once per sweep
- `logs/<drug>_<solver>_np<N>_rep<M>.log` — stdout/stderr for each individual run

The benchmark stamps every CSV row with host OS, CPU model, visible RAM, container architecture, and whether Docker is running under CPU emulation. This is especially important on Apple Silicon hosts because the Docker service is intentionally pinned to `linux/amd64`; emulated timings should not be compared against native Linux runs.

---

## Running a Simulation on an HPC Cluster (PBS/Torque)

1. Navigate to the working directory:
   ```bash
   cd bin/simulation_cpu
   ```

2. Edit `param.txt` (same as above).

3. Submit to the queue:
   ```bash
   qsub exec_pbs.sh
   ```

The PBS script requests 1 node with 10 processors (`ppn=10`) by default. Adjust `#PBS -l nodes=1:ppn=N` and `number_of_cpu` in `param.txt` to match your allocation.

To stop a running local simulation:
```bash
bash kill_bash.sh
```

---

## Post-processing & Reporting

From the `simulation_cpu` working directory, after a simulation completes:

```bash
bash generate_report.sh
```

This script performs four steps:

1. **Time-series plots** — calls `plot_time_series.py` to generate plots of Vm, dVm/dt, Cai, and key ionic currents for each concentration.
2. **Feature plots** — calls `plot_features.py` to produce dose–response plots for all 16 CiPA biomarkers.
3. **Feature aggregation** — calls `combine_and_average_features.py` to merge per-core CSV fragments and compute per-concentration sample averages.
4. **Median & confidence intervals** — calls `compute_mean_ci.py` to compute bootstrap median and 95% CI across concentrations, writing `<drug>-median-ci.csv`.

---

## Output Files

Results are written under `./results/<concentration>/`:

| File pattern | Contents |
|---|---|
| `<drug>_<conc>_time_series_smp<N>.csv` | Full time-series (Vm, Cai, currents) for sample N |
| `<drug>_<conc>_features_core<K>.csv` | CiPA feature values computed by MPI core K |
| `<drug>_<conc>_initial_values_smp<N>.csv` | Steady-state initial conditions for sample N |
| `<drug>-<conc>-features-all.csv` | Combined features for all samples at a concentration |
| `sample-averages-across-concentrations.csv` | Mean feature value per sample across concentrations |
| `<drug>-median-ci.csv` | Bootstrap median + 95% CI for each feature |
| `<drug>-depol-failure-count.csv` | Count of depolarisation failures per concentration |

---

## Python Scripts

All Python utilities live in [`bin/scripts/`](bin/scripts/).

| Script | Purpose |
|---|---|
| `plot_time_series.py` | Plot Vm, dVm/dt, Cai, INaL, ICaL, IKr, etc. vs. time |
| `plot_features.py` | Plot 16 CiPA features as dose–response curves |
| `combine_and_average_features.py` | Merge per-core feature CSVs; compute sample averages |
| `compute_mean_ci.py` | Bootstrap median and 95% CI from averaged samples |
| `convert_params.sh` | Convert `param.txt` keys to camelCase for API use |
| `create_concs_directories.sh` | Create output directory tree per concentration |
| `create_escaped_b64_file.sh` | Escape special characters and Base64-encode a file |
| `zip_files.sh` / `unzip_files.sh` | Zip/unzip simulation output files |

---

## AI Report Generation

DrugSimulation includes an optional pipeline to automatically generate a scientific manuscript (LaTeX/PDF) using the [Perplexity](https://www.perplexity.ai/) Sonar Pro API.

After running `generate_report.sh`, execute:

```bash
bash exec_perplexity_report_generator.sh
```

This pipeline:
1. Converts and Base64-encodes the parameter file, median-CI CSV, time-series CSVs, and a LaTeX template.
2. Sends structured prompts to the Perplexity API for each manuscript section (Methods, Results parts 1–7, Discussion) in both English and Korean.
3. Assembles the API responses into a complete LaTeX document.

Prompt templates are stored in [`bin/simulation_cpu/requests/`](bin/simulation_cpu/requests/). Responses are saved to `bin/simulation_cpu_postprocessing/responses/`.

> **Note:** A Perplexity API key must be configured in the environment before running the report generator.
