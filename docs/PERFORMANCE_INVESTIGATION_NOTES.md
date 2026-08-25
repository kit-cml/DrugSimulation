# Why single-core / single-sample runs are slow

Investigation date: 2026-07-22
Branch: `feat/add-docker`

## Question

A single simulation "task" (one hill/herg sample, one drug concentration, `np=1`,
no MPI parallelism) should in principle be a quick single-core job. In practice
some runs take from tens of seconds to several hours. This note documents the
causes found, ranked by size/actionability.

## Findings

### 1. Compiled with `-O0` (no optimization)

[Makefile:40](../Makefile#L40):

```make
BASE_CXXFLAGS := -std=c++11 -g -O0
```

The binary ships with zero compiler optimization and debug symbols. The
simulation core is a tight scalar ODE loop executed potentially billions of
times per run (see #4) — exactly the code shape that benefits most from
`-O2`/`-O3` (typically 10-30x for this kind of workload). This is the
single biggest, easiest, most fixable lever found.

### 2. Container runs x86_64 under emulation on Apple Silicon

[docker-compose.yml:3](../docker-compose.yml#L3):

```yaml
platform: linux/amd64
```

Confirmed via `benchmark_results/20260722_113116/environment.txt`:

```
host_os=Darwin
host_arch=arm64
host_cpu_model=Apple M3 Pro
container_arch=x86_64
emulated=yes
```

Docker Desktop is running the x86_64 container through QEMU translation on an
arm64 host. Scalar floating-point-heavy code under emulation carries a
further real overhead on top of #1.

### 3. Euler solver needs a tiny fixed step over the full real protocol

From `modules/insilico.cpp`:

```
tmax = number_pacing * cycle_length
```

Default `param.txt`: `number_pacing = 1000`, `cycle_length = 2000` (ms) →
`tmax = 2,000,000 ms` of simulated time **per concentration, per sample**.

With `solver_type = Euler`, the timestep is fixed at `time_step_min` for the
whole run (no adaptivity). `docker/linux/generate_ground_truth.sh` notes
Euler "must stay small for stability," e.g. `time_step_min = 0.001`–`0.005`.
At `dt = 0.001`, that's **~2 billion sequential Euler steps** for one
concentration — inherently serial; parallelism across samples/concentrations
doesn't help a single one.

CVode (the actual default `solver_type` in `param.txt`) uses adaptive
stepping and can take much larger steps when the system isn't stiff, which is
presumably why it's the shipped default rather than Euler.

### 4. Non-linear blowup beyond step-count math (bepridil @ 132 µM)

Same drug (bepridil), same solver (Euler), same real pacing (1000), same
`np=1` — only the fixed step size differs:

| Run | dt | concentrations | elapsed | per-concentration |
|---|---|---|---|---|
| `ground_truth_results/20260722_114811/logs/bepridil.log` | 0.005 | 4 (33,66,99,132) | 22.99 s | ~5.75 s |
| `benchmark_results/20260722_113116/logs/bepridil_Euler_np1_rep1.log` | 0.001 (5x smaller) | 2 (33,132) | 16,011.45 s | ~8,006 s |

A 5x smaller step should cost ~5x more time (more steps, same cost/step).
Observed: **~1,400x** more per concentration. Something other than step count
is dominating.

Correlated evidence:
- The slow run logs `WARNING: ... IS NOT DEPOLARIZED PROPERLY` at conc=132.
- The time series the investigation started from,
  `bin/simulation_cpu/results/132.00/bepridil_132.00_time_series_smp0.csv`,
  never fires an action potential — `Vm` sits flat at `-88.05 mV` for the
  entire 2000 ms trace.

**Working hypothesis (unconfirmed):** in this failing-to-fire / near-resting
regime, some gating variables may spend long stretches at subnormal/denormal
magnitudes. Scalar FP arithmetic on denormals can be 10-100x slower on
typical x86 FPUs (without flush-to-zero), which combined with emulation (#2)
could plausibly account for the gap. Not yet verified — would need
profiling or a flush-to-zero flag test to confirm.

### 5. A 3-hour-old run was caught in progress during this investigation

At the time of investigation, container `drugsim-linux` showed:

```
PID 492  99.9% CPU  177:52 TIME  drugsim_CiPAORdv1.0 -input_deck .docker-test-param.452.txt
```

Started 09:04 UTC, still running ~3 hours later. That input deck
(`.docker-test-param.452.txt`, a leftover from `run_only.sh`'s per-PID temp
file) specifies: `solver_type = Euler`, `time_step_min = 0.001`,
`number_pacing = 1000` (full real protocol), drug `chlorpromazine` at two
concentrations (38, 152 µM), one hill/herg sample, `np=1` — i.e. the same
expensive shape as #3/#4.

## Recommended next steps

1. Rebuild with `-O2` (or `-O3`) and re-run the same benchmark drugs to
   quantify the win before touching anything else.
2. For quick single-sample smoke checks, don't run the full
   `number_pacing=1000` Euler-at-0.001 protocol — use `solver_type=CVode`
   (adaptive) or override `PACING`/`PACING_WRITE` down, as
   `docker/linux/generate_ground_truth.sh` and `benchmark.sh` already support.
3. Consider a native `linux/arm64` build/container to drop the emulation tax
   on Apple Silicon hosts, if libCML and other deps support it.
4. If pursuing #4's denormal hypothesis: profile the bepridil/132µM Euler
   run (e.g. `perf stat`, or compile with `-ffast-math`/FTZ enabled) to see
   if denormal handling is actually the hot spot.
