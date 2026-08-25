# ARM64 Clang Compilation & Segmentation Fault Fix Report

## Overview
This document outlines the root cause analysis, fix plan, and implementation details for resolving Clang compilation errors on ARM64 and postprocessing segmentation faults in `libCML` and `DrugSimulation`.

---

## 1. Clang Namespace Collision (`tm` vs `struct tm`)

### Detailed Root Cause
- **Standard Library Header Pollution**: Apple Clang on macOS/ARM64 uses libc++ standard library headers. When standard headers like `<vector>`, `<string>`, or `<iostream>` are included, they internally cascade-include `<locale>` and `<ctime>`.
- **Global Identifier Clash**: `<ctime>` declares `struct tm` representing calendar time (containing fields like `tm_sec`, `tm_min`, `tm_hour`, etc.). 
- **Enum Member Collision**: Code-generated cell models in `libCML` (such as `Ohara_Rudy_2011`, `Shorten_2007`, etc.) defined enumerations with an explicit member identifier named `tm`. Under C++ rules, non-type enum constants pollute the enclosing scope or directly conflict with global type names when imported via `using` declarations or namespace resolutions in libc++ locale facets (`std::locale`).
- **Compiler Diagnostic**: Clang strict parsing threw fatal errors:
  ```text
  error: must use 'struct' tag to refer to type 'tm' in this scope
  note: struct 'tm' is hidden by a non-type declaration of 'tm' here
  ```

### Detailed Resolution
- **Global Refactoring**: Performed systematic renaming of the enum identifier `tm` to `cell_tm` across all generated cell model header files (`enum_*.hpp`) and corresponding `.cpp` calculation implementations in `libCML/src/types/cellmodels/`.
- **Scope Isolation**: Ensured enum constants do not collide with POSIX/standard C library structs by prefixing domain-specific abbreviations (`cell_`).

---

## 2. Initial Values Segmentation Fault on ARM64

### Detailed Root Cause
- **Vector Out-of-Bounds Access**: During long pacing simulations (`insilico`), the model extracts state vectors to populate `p_features.initial_values`. 
- **Pacing Count Mismatch**: Under certain ARM64 compilation flags and optimizations, `end_of_cycle_funct` recorded reference features at pacing offsets that did not guarantee `temp_features.initial_values` was fully populated before being copied to `p_features`.
- **Dereferencing Invalid Pointers**: When writing out states via `fprintf(fp_initial_values, ...)` in `modules/insilico.cpp`, if `initial_values` size was mismatched or empty while `vm_peak` condition evaluated to true, it resulted in a null-pointer or out-of-bounds vector subscript segmentation fault (`SIGSEGV`).

### Detailed Resolution
- **Hoisted State Capture**: Shifted the `initial_values` population logic inside `end_of_cycle_funct` (`modules/insilico.cpp`) to occur deterministically when pacing threshold is met:
  ```cpp
  if (*pace_count >= p_param->number_pacing - 252) {
    temp_features.initial_values.insert(
      temp_features.initial_values.begin(), 
      p_cell->STATES, 
      p_cell->STATES + p_cell->states_size
    );
  }
  ```
- **Strict Size Validation Guards**: Added explicit preconditions before writing states to files in both `modules/insilico.cpp` and `modules/postprocessing.cpp`:
  ```cpp
  if (p_features.vm_peak >= 0 && p_features.initial_values.size() == p_cell->states_size) {
    for (short idx = 0; idx < p_cell->states_size; idx++) {
      fprintf(fp_initial_values, "%.32lf\n", p_features.initial_values[idx]);
    }
  }
  ```
  This prevents undefined behavior and segmentation faults even if simulation convergence or pacing bounds vary.

---

## Verification
- Verified compilation and test runs.
- Confirmed stable execution on ARM64 platforms.
