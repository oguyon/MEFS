---
trigger: always_on
---

# Common Agent Mistakes

Consolidated checklist of pitfalls that AI agents frequently hit when generating code for MEFS. Check this list before finalizing any generated code.

## CMake & Build System
1. **Forgetting to add new source files to CMake.**
   If you add a new source file, it must be added to the appropriate target in `CMakeLists.txt`.
2. **Not linking required dependencies.**
   When compiling new executables, link against standard libraries and required dependencies if they are used (e.g., `${CFITSIO_LIBRARIES}`, `${FFTW3_LIBRARIES}`, `m`).

## Memory & Resource Management
3. **Allocation inside compute loops.**
   Never call `malloc()`, `calloc()`, or `realloc()` inside core computation loops. Allocate buffers during setup/initialization and reuse them.
4. **Leaking file descriptors or handles.**
   Always close FITS handles, FFTW plans, and file pointers on error paths. Use the `goto cleanup` pattern to free resources in reverse allocation order.

## Compiler & Math Correctness
6. **Implicit double promotions.**
   Ensure that float math uses float variants (`sqrtf`, `powf`) and float literals (`0.5f`) when working with single-precision floats, or standard double variants when using double-precision floats.
7. **Type mismatches in loop indices.**
   Always match the loop index type to the bound variable type (e.g., `for (uint64_t ii = 0; ii < xysize; ii++)` when `xysize` is `uint64_t`) to prevent breaking compiler SIMD auto-vectorization.
8. **Implicit header includes.**
   Every `.c` file must include exactly the headers it uses. Do not rely on header side-effects.
9. **Lines > 100 characters.**
   Limit line length in C source code, scripts, and documentation files to 100 characters.
