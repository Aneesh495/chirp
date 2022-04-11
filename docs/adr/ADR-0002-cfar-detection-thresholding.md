# ADR-0002: 2D CFAR Detection and Adaptive Thresholding Strategy

## Status
Accepted

## Context
Raw Range-Doppler maps generated from FMCW radar cubes contain non-stationary clutter, thermal noise floors, and multi-target reflections spanning 100 dB of dynamic range. Standard fixed-threshold detection fails because thermal noise power fluctuates across temperature, interference raises local noise floors, and ground clutter varies by orders of magnitude across physical environments.

Constant False Alarm Rate (CFAR) processors adaptively compute detection thresholds on a cell-by-cell basis by estimating background noise and clutter power in sliding reference windows surrounding each cell under test (CUT). However, conventional 2D sliding windows suffer from two primary operational failure modes:
1. Masking: Closely spaced adjacent targets contaminate each other's training reference cells, raising the threshold and suppressing legitimate detections.
2. Clutter edge false alarms: Moving across sharp physical boundaries (such as road edges or tunnel entrances) causes classical cell-averaging detectors to trigger false detections.

## Decision
We implement a unified, multi-algorithm 2D CFAR engine with selectable statistical estimators:
1. Cell-Averaging (CA-CFAR): Optimal maximum-likelihood estimator for homogeneous, independent and identically distributed Rayleigh clutter and thermal noise.
2. Greatest-Of (GO-CFAR): Evaluates leading and lagging half-windows independently and selects the maximum, preventing false alarm inflation across clutter step edges.
3. Smallest-Of (SO-CFAR): Evaluates leading and lagging half-windows independently and selects the minimum, resolving closely spaced adjacent targets without mutual masking.
4. Ordered-Statistic (OS-CFAR): Sorts reference cells in ascending order and selects the $k$-th ranked sample ($k \approx 0.75 N$), completely immune to up to $N - k$ interfering target reflections.
5. Integral Image Prefix-Sum Acceleration: In CA-CFAR, we construct a 2D integral image (summed-area table) in $O(N_R \cdot N_V)$ preprocessing time, allowing arbitrary rectangular sliding-window noise sums to be computed in exact $O(1)$ constant time per cell.

## Consequences
### Positive
- Strict statistical false alarm rate control ($P_{\text{fa}} = 10^{-4}$ to $10^{-6}$) maintained across homogeneous, multi-target, and boundary scenes.
- 2D Integral image optimization delivers $25,000+$ frames per second in compiled C MEX kernels.
- Flexible algorithmic abstraction allows runtime switching between low SNR loss (CA) and multi-target robustness (OS/SO).

### Negative
- OS-CFAR requires sorting reference cells, which has higher computational complexity ($O(N \log N)$) than linear accumulation.
- 2D sliding windows require guard cells around the CUT to prevent energy spillage from broad mainlobes from contaminating the noise estimate.
