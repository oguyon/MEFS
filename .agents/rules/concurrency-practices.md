---
description: Thread safety and concurrency patterns for MEFS.
---

# Concurrency Practices

MEFS supports multi-threading via OpenMP. Follow these patterns to avoid race conditions.

## OpenMP Parallelization
- Ensure loops parallelized via `#pragma omp parallel for` do not have race conditions on shared state.
- Use atomic updates for shared statistics or counters:
  ```c
  #ifdef _OPENMP
  #pragma omp atomic
  #endif
  state->framedist_calls++;
  ```
- Keep loop index variables private to each thread (declaring them inside the `for` statement achieves this automatically in C99/C11/C17).

## Volatile Keyword
- Use `volatile sig_atomic_t` for signal handling flags (e.g., `stop_requested`):
  ```c
  extern volatile sig_atomic_t stop_requested;
  ```
- Do not use `volatile` on data arrays or pointers, as it disables compiler optimization and SIMD vectorization entirely.
