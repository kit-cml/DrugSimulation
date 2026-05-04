# Solution 2 (macOS Native) Self-Execution and Rollback Guide

## Goal
Build and run DrugSimulation on macOS using one consistent clang-based toolchain, while keeping safe rollback points after every risky step.

This guide assumes:
- DrugSimulation repo path: `/Users/iganarendra/DrugSimulation`
- libCML repo path: `/Users/iganarendra/libCML`
- macOS with Xcode command line tools installed

---

## High-Level Decision Tree

1. If both repos are clean, continue.
2. Else, snapshot/stash/backup first, then continue.
3. Build SUNDIALS 5.x from source inside libCML.
4. If SUNDIALS build fails, stop and fix dependency/toolchain issue before touching simulation code.
5. Patch `tm` enum collisions in libCML (clang compatibility).
6. Build libCML with clang.
7. If libCML build fails, rollback patch set and retry from a clean point.
8. Build DrugSimulation with matching clang path.
9. If link fails with ABI/symbol mismatch, verify both sides were truly rebuilt with clang and SUNDIALS 5.x.
10. Run simulation (skip Python/report stage).
11. If runtime fails, inspect logs and rollback only the minimal layer needed.

---

## Step 0: Preflight Checks

Run:

```bash
xcode-select -p
clang++ --version
cmake --version
make --version
mpicxx --version
```

If `xcode-select -p` fails:

```bash
xcode-select --install
```

If `cmake` is missing:

```bash
brew install cmake
```

Rollback impact: none (tool installation only).

---

## Step 1: Safety Snapshot (Mandatory)

### 1.1 Check git status in both repos

```bash
cd /Users/iganarendra/DrugSimulation
git status --short

cd /Users/iganarendra/libCML
git status --short
```

### 1.2 Branch by condition

- If status is clean in both repos:
  - Create safety tags/branches anyway.
- Else:
  - Stash or commit before proceeding.

Recommended:

```bash
cd /Users/iganarendra/DrugSimulation
git checkout -b backup/pre-solution2-$(date +%Y%m%d-%H%M%S)

git add -A
git commit -m "backup: pre solution2 attempt" || true

git stash push -u -m "pre-solution2-drugsim" || true

cd /Users/iganarendra/libCML
git checkout -b backup/pre-solution2-$(date +%Y%m%d-%H%M%S)

git add -A
git commit -m "backup: pre solution2 attempt" || true

git stash push -u -m "pre-solution2-libcml" || true
```

### 1.3 Optional file-system backup (extra safety)

```bash
cd /Users/iganarendra
mkdir -p backup_solution2

tar -czf backup_solution2/DrugSimulation-pre-solution2.tar.gz DrugSimulation
tar -czf backup_solution2/libCML-pre-solution2.tar.gz libCML
```

Rollback impact: complete restore possible even if git history gets messy.

---

## Step 2: Build Native SUNDIALS 5.x for macOS (inside libCML)

Reason: bundled prebuilt archives are Linux ELF, not usable on macOS.

### 2.1 Configure build dir

```bash
cd /Users/iganarendra/libCML/libs/sundials-5.7.0
rm -rf build-macos-clang
mkdir -p build-macos-clang
cd build-macos-clang
```

### 2.2 CMake configure

```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=clang \
  -DCMAKE_CXX_COMPILER=clang++ \
  -DBUILD_SHARED_LIBS=OFF \
  -DEXAMPLES_ENABLE=OFF \
  -DSUNDIALS_BUILD_WITH_MONITORING=OFF
```

If configure fails:
- If message indicates missing CMake modules/toolchain:
  - confirm Xcode tools and CMake install, then rerun.
- Else:
  - stop here and fix before moving on.

### 2.3 Build

```bash
cmake --build . -j$(sysctl -n hw.ncpu)
```

### 2.4 Verify generated static libs are Mach-O

```bash
find . -name "libsundials*.a" | head -n 10
file src/cvode/libsundials_cvode.a
file src/nvector/serial/libsundials_nvecserial.a
```

Expected: output contains `current ar archive random library` on macOS (not ELF object references from Linux builds).

Rollback for Step 2:

```bash
cd /Users/iganarendra/libCML/libs/sundials-5.7.0
rm -rf build-macos-clang
```

---

## Step 3: Patch `tm` Enum Collisions in libCML (clang fix)

Reason: some model headers use `tm` as enum identifier; clang conflicts with `struct tm` contexts more strictly.

### 3.1 Locate collisions

```bash
cd /Users/iganarendra/libCML
rg -n "\btm\b" cellmodels | head -n 200
```

### 3.2 Apply minimal rename pattern

Use one consistent replacement for enum token, for example `tm_state`.

Important rules:
- Rename only enum member identifiers and direct references.
- Do not rename unrelated variable names unless compilation requires it.
- Keep changes scoped to files that fail or are clearly part of the collision set.

### 3.3 Verify no accidental broad replacement

```bash
git diff --stat
git diff -- cellmodels
```

If diff looks too broad:

```bash
git restore --source=HEAD --staged --worktree -- cellmodels
```

Then redo with targeted edits.

Rollback for Step 3:

```bash
cd /Users/iganarendra/libCML
git restore --source=HEAD --staged --worktree -- .
# or selective restore for touched files only
```

---

## Step 4: Build libCML with clang + native SUNDIALS includes/libs

### 4.1 Clean old build artifacts

```bash
cd /Users/iganarendra/libCML
rm -rf build
mkdir -p build
```

### 4.2 Build command (template)

Adjust include/lib paths if your SUNDIALS build output differs.

```bash
make \
  CXX=clang++ \
  CC=clang \
  CXXFLAGS="-fPIC -fpermissive -std=c++11 -I./libs/sundials-5.7.0/include -I./libs/sundials-5.7.0/build-macos-clang/include" \
  LDFLAGS="-L./libs/sundials-5.7.0/build-macos-clang/src/cvode -L./libs/sundials-5.7.0/build-macos-clang/src/nvector/serial" \
  build/libcml.a
```

If build fails:
- If failure is still `tm`-related:
  - go back to Step 3 and finish renames in remaining files.
- If failure is missing SUNDIALS headers:
  - verify include path and generated headers in build folder.
- If failure is unrelated model source issue:
  - fix one model source at a time, rebuild.

### 4.3 Verify archive exists and is current

```bash
ls -lh build/libcml.a
```

Rollback for Step 4:

```bash
cd /Users/iganarendra/libCML
rm -rf build
# optionally restore source edits if needed
# git restore --source=HEAD --staged --worktree -- .
```

---

## Step 5: Wire DrugSimulation to Correct CML/SUNDIALS Artifacts

Goal: ensure DrugSimulation links against the newly rebuilt clang-compatible `libcml.a` and matching native SUNDIALS static libs.

### 5.1 Build clean

```bash
cd /Users/iganarendra/DrugSimulation
make clean || true
```

### 5.2 Build with explicit linker overrides (preferred, less invasive)

Use your Makefile variable names as applicable:

```bash
make CXX=mpicxx \
  LDFLAGS="-g ../libCML/build/libcml.a \
  ../libCML/libs/sundials-5.7.0/build-macos-clang/src/cvode/libsundials_cvode.a \
  ../libCML/libs/sundials-5.7.0/build-macos-clang/src/nvector/serial/libsundials_nvecserial.a \
  -lcurl -ljson-c"
```

If your build needs additional SUNDIALS components, append the exact generated `.a` paths from your local build.

If link fails with `std::__1` vs `std::__cxx11` symbols:
- libCML was not fully rebuilt with clang, or stale objects were reused.
- Do full clean of libCML build dir and rebuild Step 4.

If link fails with SUNDIALS symbols missing:
- wrong SUNDIALS lib set or wrong order.
- inspect undefined symbol names and add corresponding SUNDIALS modules from your build output.

Rollback for Step 5:

```bash
cd /Users/iganarendra/DrugSimulation
make clean || true
git restore --source=HEAD --staged --worktree -- Makefile
```

(If you only used command-line overrides and did not edit Makefile, rollback is just `make clean`.)

---

## Step 6: Run the Simulation (Skip Python)

From simulation folder:

```bash
cd /Users/iganarendra/DrugSimulation/bin/simulation_cpu
```

Run executable directly or with MPI depending on your project expectation.

Examples:

```bash
../drugsim_CiPAORdv1.0
```

or

```bash
mpirun -np 10 ../drugsim_CiPAORdv1.0
```

Use `param.txt` already present in this folder.

If run fails immediately:
- Check binary exists and has execute permission.
- Confirm current working directory is `bin/simulation_cpu`.
- Validate input files referenced by `param.txt` exist.

Rollback for Step 6: none required (runtime only).

---

## Step 7: Validation Checklist

Mark complete only when all are true:

- `libCML/build/libcml.a` exists and was rebuilt in current attempt.
- DrugSimulation binary (for your target model) is produced.
- No `std::__1` vs `std::__cxx11` unresolved symbol errors.
- Simulation starts and progresses through pacing/concentrations.
- Output files appear in expected output directories.

---

## Fast Rollback Recipes

## A) Undo local edits only (keep history)

```bash
cd /Users/iganarendra/libCML
git restore --source=HEAD --staged --worktree -- .
rm -rf build libs/sundials-5.7.0/build-macos-clang

cd /Users/iganarendra/DrugSimulation
git restore --source=HEAD --staged --worktree -- .
make clean || true
```

## B) Return to stash snapshot

```bash
cd /Users/iganarendra/libCML
git stash list
git stash apply stash@{0}

cd /Users/iganarendra/DrugSimulation
git stash list
git stash apply stash@{0}
```

## C) Nuclear restore from tar backups

```bash
cd /Users/iganarendra
mv DrugSimulation DrugSimulation.failed.$(date +%Y%m%d-%H%M%S)
mv libCML libCML.failed.$(date +%Y%m%d-%H%M%S)

tar -xzf backup_solution2/DrugSimulation-pre-solution2.tar.gz
tar -xzf backup_solution2/libCML-pre-solution2.tar.gz
```

Use this only if git state is beyond easy repair.

---

## Practical Execution Strategy (Recommended)

1. Perform Steps 0-2 fully.
2. Commit checkpoint in libCML:

```bash
cd /Users/iganarendra/libCML
git add -A
git commit -m "checkpoint: native sundials built"
```

3. Perform Step 3 patching and commit:

```bash
git add -A
git commit -m "fix: rename tm enum identifiers for clang"
```

4. Perform Step 4 and commit if build succeeds.
5. Perform Step 5 in DrugSimulation, preferably with CLI overrides first (avoid permanent Makefile edits).
6. Run Step 6.

This gives clean, reversible checkpoints at each major stage.

---

## When to Stop and Ask for Help

Stop immediately if:
- You see massive unintended diffs in libCML model files.
- You get a new class of compiler errors unrelated to `tm` or SUNDIALS.
- You accidentally mix Homebrew SUNDIALS 7.x headers with SUNDIALS 5.x libs.

At that point, rollback to last checkpoint and re-run from the previous stable step.
