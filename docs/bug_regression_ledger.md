# Bug Regression and Mathematical Defect Ledger

## Overview

This engineering ledger records all defects, mathematical edge cases, and runtime compatibility issues encountered during the development of CHIRP. In accordance with zero-defect verification principles, each entry details the failure mode, mathematical root cause, exact resolution applied, and permanent automated test suite validating the fix against regressions.

---

### Defect 01: Method Chaining Indexing Conflict in GNU Octave
- **Component**: `+chirp/+dsp/RangeDopplerMap.m`
- **Location**: Line 48 (`findPeaks` method)
- **Symptom**: Runtime syntax crash `error: nonCoherentIntegration: function called with too many inputs`.
- **Root Cause**: In GNU Octave, invoking `obj.nonCoherentIntegration(:)` is parsed by the language grammar as an argument pass `nonCoherentIntegration(obj, ':')` rather than indexing the returned matrix.
- **Resolution**: Assigned the intermediate matrix result to a local variable `pwr = obj.nonCoherentIntegration();` before computing `max(pwr(:))`.
- **Regression Verification**: `tests/unit/test_dsp.m` (Assertion 6).

---

### Defect 02: Non-Strict Peak Detector Flooding Plateau With False Targets
- **Component**: `+chirp/+dsp/RangeDopplerMap.m`
- **Location**: Line 78 (`findPeaks` dilation loop)
- **Symptom**: RangeDopplerMap returned 15,030 peaks on a synthetic flat noise floor.
- **Root Cause**: The 2D neighborhood maximum comparison used `>=` without masking the center cell under test. On flat plateaus or identical adjacent noise values, every cell equaled the neighborhood maximum, declaring every bin a spectral peak.
- **Resolution**: Explicitly masked the center cell `neighborhood(d+1, d+1) = -Inf;` and enforced strict inequality `val > max(neighborhood(:))` combined with descending power sorting.
- **Regression Verification**: `tests/unit/test_dsp.m` (Assertion 6: `findPeaks must return exactly 1 peak on plateau`).

---

### Defect 03: Root-MUSIC Conjugate Phase Reversal in Steering Vector Angle Extraction
- **Component**: `+chirp/+aoa/RootMusicEstimator.m`
- **Location**: Line 64
- **Symptom**: Targets simulated at $+15.0^\circ$ and $-12.0^\circ$ reconstructed at $-15.0^\circ$ and $+12.0^\circ$.
- **Root Cause**: The forward spatial steering vector uses $a_k = \exp(-j 2\pi d k \sin(\theta) / \lambda)$. Polynomial root angle extraction via $\text{angle}(z_i)$ yielded $+2\pi d \sin(\theta) / \lambda$, which requires negation to map back to forward physical angles.
- **Resolution**: Inverted the extracted phase angle: `sinTheta = -(phases * lambda) / (2.0 * pi * d);`.
- **Regression Verification**: `tests/unit/test_aoa.m` (Assertion 4: Root-MUSIC theta1 and theta2 accuracy).

---

### Defect 04: Octave Class Parsing Bug on Empty Object Initialization
- **Component**: `+chirp/+tracking/GnnTracker.m`, `+chirp/+tracking/JpdaTracker.m`
- **Location**: Property initialization
- **Symptom**: `error: no such method or property 'empty' in class chirp.tracking.RadarTrack`.
- **Root Cause**: GNU Octave does not automatically synthesize the static `.empty(1, 0)` method for nested package classdefs under all operating environments.
- **Resolution**: Initialized property as `Tracks = [];` and guarded array operations with `isempty(obj.Tracks)`.
- **Regression Verification**: `tests/unit/test_tracking.m` (Assertion 1).

---

### Defect 05: 2D Matrix Index Linear Subscript Out-of-Bounds in Micro-Doppler Feature Extractor
- **Component**: `+chirp/+microdoppler/MicroDopplerFeatures.m`
- **Location**: Line 52
- **Symptom**: `error: freqAxis(477): out of bound 256`.
- **Root Cause**: `find(spectrogramPower == maxVal)` returned a 2D linear index into the $256 \times 64$ matrix. The linear index was passed directly into the 1D frequency axis vector `freqAxis` (length 256), throwing an out-of-bounds exception.
- **Resolution**: Decomposed 2D linear index into row and column coordinates using `[pkRow, ~] = find(spectrogramPower == maxVal, 1);` then indexed `freqAxis(pkRow)`.
- **Regression Verification**: `tests/unit/test_microdoppler.m` (Assertion 2).

---

### Defect 06: Zero-Cadence Penalty Misclassification of Constant Velocity Pedestrians
- **Component**: `+chirp/+microdoppler/KinematicClassifier.m`
- **Location**: Lines 72-84
- **Symptom**: Simulated pedestrian walking signature misclassified as drone or cyclist.
- **Root Cause**: When walking cadence estimation had low confidence or was near DC, an unconstrained penalty heavily depressed pedestrian likelihood while unconstrained rotor bounds elevated drone scores.
- **Resolution**: Conditioned cadence penalties strictly on high confidence detections, constrained drone rotor thresholds to speeds exceeding 25 m/s, and prioritized torso Doppler extent matching walking bounds (0.5 to 2.5 m/s).
- **Regression Verification**: `tests/unit/test_microdoppler.m` (Assertion 4).

---

### Defect 07: Speed of Light Hardcoded Approximation Mismatch in Range Resolution
- **Component**: `tests/unit/test_waveform.m`
- **Location**: Line 22
- **Symptom**: Assertion failure expecting range resolution to equal literal `0.300000000 m`.
- **Root Cause**: Exact physical constant $c = 299,792,458\text{ m/s}$ divided by $2 \times 500\text{ MHz}$ produces $0.299792458\text{ m}$. Comparing against rounded $0.3\text{ m}$ exceeded tolerance $10^{-6}$.
- **Resolution**: Computed theoretical resolution using exact fundamental speed of light `c / (2.0 * B)`.
- **Regression Verification**: `tests/unit/test_waveform.m` (Assertion 1).

---

### Defect 08: Monopulse Single-Snapshot Angle Extraction Null on Subarray Pair
- **Component**: `+chirp/+aoa/MonopulseEstimator.m`
- **Location**: Lines 35-50
- **Symptom**: Monopulse estimation returned $0.00^\circ$ for a target located at $2.50^\circ$ off boresight.
- **Root Cause**: The monopulse ratio was taking only the real part $\text{Re}(\Delta \Sigma^*) / |\Sigma|^2$. For two spatially separated subarrays, the phase center displacement creates phase-comparison monopulse, where the angle information resides entirely in the imaginary component $\text{Im}(\Delta / \Sigma)$.
- **Resolution**: Evaluated both real and imaginary components, dynamically selecting the dominant projection to support amplitude-comparison and phase-comparison monopulse seamlessly.
- **Regression Verification**: `tests/unit/test_aoa.m` (Assertion 6).

---

### Defect 09: Lexicographical Reordering in BFS Queue for Doppler-Aware DBSCAN
- **Component**: `+chirp/+pointcloud/VelocityAwareDbscan.m`
- **Location**: Lines 78-86
- **Symptom**: DBSCAN separated two clusters into 4 fragmented clusters; identical points were assigned different cluster IDs.
- **Root Cause**: During neighbor expansion, the queue update performed `seedQueue = unique([seedQueue, newPoints]);`. The `unique` function sorts elements lexicographically. This reordered elements behind and ahead of `qIdx`, causing unvisited points to jump past the pointer and form new clusters prematurely.
- **Resolution**: Replaced `unique` with strict FIFO queue order and boolean membership array `inQueue(1:N)`, preserving sequence and guaranteeing deterministic $O(1)$ membership checks.
- **Regression Verification**: `tests/unit/test_pointcloud.m` (Assertion 4).

---

### Defect 10: Nonconformant Vector-Matrix Subtraction in Hermitian Covariance Invariant
- **Component**: `tests/audit/test_mathematical_invariants.m`
- **Location**: Line 20
- **Symptom**: `error: operator -: nonconformant arguments (op1 is 64x1, op2 is 8x8)`.
- **Root Cause**: Flattened column vector `R(:)` was subtracted from 2D reshaped matrix `conj(reshape(R.', size(R)))`.
- **Resolution**: Replaced with standard matrix conjugate transpose subtraction: `max(max(abs(R - R')))`.
- **Regression Verification**: `tests/audit/test_mathematical_invariants.m` (Invariant 2).

---

### Defect 11: Stepped-Frequency SFCW Range Aliasing Beyond Unambiguous Range
- **Component**: `tests/unit/test_waveform.m`
- **Location**: Line 40
- **Symptom**: Target at 12.0 meters reconstructed at 0.0 meters.
- **Root Cause**: Stepped frequency continuous wave has maximum unambiguous range $R_{\text{max}} = c / (2 \Delta f)$. With $\Delta f = 25\text{ MHz}$, $R_{\text{max}} = 6.0\text{ meters}$. A test target at $12.0\text{ m}$ was exactly $2 \times R_{\text{max}}$, aliasing modulo $R_{\text{max}}$ to $0.0\text{ m}$.
- **Resolution**: Set frequency step to $\Delta f = 10\text{ MHz}$ ($R_{\text{max}} \approx 15.0\text{ m}$), placing the $12.0\text{ m}$ target inside unambiguous range.
- **Regression Verification**: `tests/unit/test_waveform.m` (Assertion 6).

---

### Defect 12: Missing 'contains' Builtin and Character Encoding Warnings in GNU Octave
- **Component**: `tests/audit/test_no_em_dashes.m`
- **Location**: Lines 15, 24
- **Symptom**: `Error: 'contains' undefined` and range conversion warnings.
- **Root Cause**: In GNU Octave, string function `contains()` is not implemented in core language; `char(hex2dec('2014'))` overflows 8-bit character bytes.
- **Resolution**: Used `~isempty(strfind(str, pattern))` and scanned UTF-8 bytes via `char([226, 128, 148])` with `fileread`.
- **Regression Verification**: `tests/audit/test_no_em_dashes.m`.
