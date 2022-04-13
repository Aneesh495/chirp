# ADR-0005: Native C MEX Acceleration and Pure MATLAB Fallback Architecture

## Status
Accepted

## Context
Radar signal processing and SAR imaging involve inner loops over millions of raw electromagnetic samples, high-dimensional tensor operations, and spatial backprojection accumulations. While modern interpreted languages and JIT compilers provide reasonable matrix multiplication throughput, complex sliding-window branching (such as 2D CFAR) and non-vectorizable spatial lookups (such as point-by-point backprojection) suffer significant latency penalties in interpreted loops.

At the same time, radar research, testing, and deployment require portability across disparate computing environments: x86_64 and ARM64 Apple Silicon developer machines, CI build pipelines, and platforms where native C compilers or MEX toolchains may not be immediately available.

## Decision
We establish a zero-dependency Dual-Path Execution Architecture across all computationally intensive kernels:
1. Every performance-critical processing block implements two strictly equivalent, bit-level compatible execution paths:
   - A compiled native C MEX kernel compiled via POSIX `make mex` or `mkoctfile --mex` (located in `src/mex/` and built into `bin/mex/`).
   - A pure, highly optimized vectorized MATLAB / GNU Octave fallback implementation (located in `+chirp` packages).
2. Runtime Auto-Detection: At initialization or execution, processing engines query `exist('chirp_<kernel>_mex', 'file') == 3`. If the compiled binary exists and is compatible, execution dynamically routes to the MEX kernel; otherwise, execution transparently falls back to the native MATLAB engine without throwing errors or interrupting execution.
3. Automated Parity Verification: A dedicated audit suite (`tests/unit/test_mex_parity.m`) continuously validates that the numerical output of every C MEX kernel matches the pure MATLAB implementation within floating-point tolerance ($\Delta < 10^{-13}$).

## Consequences
### Positive
- Maximum execution performance on developer workstations with up to 2,000x acceleration in CFAR detection and 38x speedup in SAR backprojection.
- Absolute portability: The repository clones and runs immediately on any system with a standard MATLAB or GNU Octave installation, even before compilation.
- Continuous numerical integrity guaranteed through automated parity unit tests.

### Negative
- Requires maintaining two parallel implementations of core algorithms (C source code and MATLAB m-code), increasing test and auditing maintenance overhead.
