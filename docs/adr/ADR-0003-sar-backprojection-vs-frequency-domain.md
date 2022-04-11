# ADR-0003: Time-Domain Backprojection versus Frequency-Domain SAR Algorithms

## Status
Accepted

## Context
Synthetic Aperture Radar (SAR) synthesizes an effective antenna aperture orders of magnitude larger than physical antenna dimensions by coherently accumulating pulse returns along a moving flight path. Image formation requires focusing phase histories collected across slow-time and fast-time into a 2D spatial reflectivity map $I(x, y)$.

Two major architectural paradigms exist:
1. Frequency-Domain Algorithms (Range-Doppler Algorithm, Omega-K, Polar Format Algorithm): Utilize Fast Fourier Transforms to achieve $O(N^2 \log N)$ computational efficiency. However, they rely on planar flight path assumptions, straight linear trajectories, constant velocity, and narrow synthetic beam angles. Non-linear flight paths, atmospheric turbulence, and high squint angles require complex Range Cell Migration Correction (RCMC) and Stolt interpolation approximations.
2. Time-Domain Backprojection (BPA): Reconstructs ground pixels individually by computing exact geometric slant ranges from flight trajectory coordinates to each ground coordinate, phase-rotating and accumulating contributions coherently.

## Decision
We select Time-Domain Backprojection (BPA) as the primary, high-fidelity SAR image formation engine, complemented by frequency-domain Polar Format and Range-Doppler engines for specialized real-time workflows:
1. Time-Domain Backprojection is implemented with exact Pythagorean 3D slant range calculation for arbitrary, non-linear, non-planar platform trajectories (Stripmap, Spotlight, and Circular SAR).
2. The backprojection reconstruction loop is parallelized and accelerated via a native C MEX kernel (`chirp_backprojection_mex`), which evaluates sub-pixel range interpolation and complex phase rotation in continuous floating-point memory.
3. Frequency-domain Polar Format Algorithm (PFA) is provided for spotlight SAR scenes where narrow angular excursions permit fast 2D polar-to-Cartesian resampling.

## Consequences
### Positive
- Perfect geometric focus without approximations, valid for wide aperture angles ($> 30^\circ$) and circular SAR configurations.
- Seamlessly accommodates arbitrary flight trajectories, platform drift, altitude variations, and non-uniform velocity without requiring separate motion compensation filtering.
- Pixel-independent computation maps directly to SIMD vectorization and parallel multicore processing, achieving over 10 Million pixels per second in C MEX.

### Negative
- Computational complexity scales as $O(N_X \cdot N_Y \cdot N_{\text{pulses}})$, which is computationally heavier than frequency-domain FFT techniques for very large grids.
- Requires accurate GPS/INS trajectory coordinates for each transmitted pulse.
