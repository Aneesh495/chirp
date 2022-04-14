# FMCW Radar Signal Processing Mathematical Formulation

## 1. Frequency-Modulated Continuous-Wave (FMCW) Waveforms

An FMCW radar periodically transmits a linear frequency ramp (chirp) whose instantaneous carrier frequency increases linearly with slope $S$ over chirp duration $T_c$:

$$f_{\text{TX}}(t) = f_c + S t, \quad 0 \le t \le T_c$$

where $f_c$ is the start carrier frequency (e.g. 77.0 GHz) and the sweep bandwidth is $B = S \cdot T_c$ (e.g. 1.0 GHz). The instantaneous transmit phase is the integral of angular frequency:

$$\phi_{\text{TX}}(t) = 2\pi \int_0^t f_{\text{TX}}(\tau) d\tau = 2\pi \left( f_c t + \frac{1}{2} S t^2 \right)$$

For a target located at range $R(t) = R_0 + v t$, the reflected echo arrives with round-trip time delay:

$$\tau(t) = \frac{2 R(t)}{c} = \frac{2 (R_0 + v t)}{c}$$

The received signal is an attenuated, delayed replica of the transmit waveform:

$$s_{\text{RX}}(t) = A_{\text{RX}} \exp\left( j 2\pi \left( f_c (t - \tau) + \frac{1}{2} S (t - \tau)^2 \right) \right)$$

## 2. De-Chirp Mixing and Intermediate Frequency (IF) Signal

The received RF signal is mixed with a replica of the transmit chirp in a quadrature mixer and low-pass filtered, producing the complex baseband Intermediate Frequency (IF) or beat signal:

$$s_{\text{IF}}(t) = s_{\text{RX}}(t) \cdot s_{\text{TX}}^*(t) = A \exp\left( -j 2\pi \left( f_c \tau + S t \tau - \frac{1}{2} S \tau^2 \right) \right)$$

Neglecting the residual video phase $\frac{1}{2} S \tau^2 \approx 0$ for short ranges, the beat frequency $f_b$ and phase $\psi$ are:

$$f_b = S \tau = \frac{2 S R_0}{c} + \frac{2 f_c v}{c}$$

$$\psi_0 = \frac{4\pi R_0}{\lambda}$$

For typical short chirps ($T_c < 50\ \mu\text{s}$), the velocity contribution to fast-time beat frequency is negligible compared to the range term, yielding:

$$f_{\text{beat}} \approx \frac{2 S R}{c} \implies R = \frac{c \cdot f_{\text{beat}}}{2 S}$$

The physical range resolution is governed strictly by the sweep bandwidth $B$:

$$\Delta R = \frac{c}{2 B}$$

## 3. Two-Dimensional Range-Doppler Processing

Radar data collected over $N_{\text{adc}}$ fast-time samples and $N_{\text{chirp}}$ slow-time chirps forms a 2D matrix $S[n, m]$.

### Step 1: Fast-Time Range FFT
A windowed discrete Fourier transform along the fast-time dimension $n$ extracts the beat frequency spectrum:

$$X_R[k, m] = \sum_{n=0}^{N_{\text{adc}}-1} w_R[n] \cdot S[n, m] \exp\left( -j \frac{2\pi k n}{N_{\text{FFT,R}}} \right)$$

Window functions $w_R[n]$ suppress sidelobes to prevent weak distant targets from being masked by strong close-in reflections:
- **Hann**: First sidelobe at -31.5 dB, 1.50 bin mainlobe width.
- **Blackman-Harris (4-term)**: First sidelobe at -92.0 dB, 2.00 bin mainlobe width.
- **Chebyshev**: Equiripple sidelobes with user-specified attenuation (e.g. -80.0 dB).

### Step 2: Slow-Time Doppler FFT
Between successive chirps separated by pulse repetition interval $T_{\text{PRI}}$, a target moving with radial velocity $v$ moves distance $\Delta R = v \cdot T_{\text{PRI}}$. This produces a slow-time phase progression:

$$\Delta \phi = \frac{4\pi \Delta R}{\lambda} = \frac{4\pi v T_{\text{PRI}}}{\lambda}$$

Taking the windowed FFT across chirp indices $m$ resolves the radial velocity spectrum:

$$X_{RD}[k, l] = \sum_{m=0}^{N_{\text{chirp}}-1} w_D[m] \cdot X_R[k, m] \exp\left( -j \frac{2\pi l m}{N_{\text{FFT,D}}} \right)$$

The Doppler frequency bin $l$ maps directly to radial velocity:

$$f_D = \frac{2 v}{\lambda} \implies v = \frac{\lambda \cdot f_D}{2}$$

The maximum unambiguous radial velocity is:

$$V_{\text{max}} = \frac{\lambda}{4 T_{\text{PRI}}}$$

And the velocity resolution is:

$$\Delta v = \frac{\lambda}{2 N_{\text{chirp}} T_{\text{PRI}}} = \frac{\lambda}{2 T_{\text{frame}}}$$

## 4. Doppler Ambiguity Resolution via Chinese Remainder Theorem

When a target velocity exceeds $V_{\text{max}}$, the measured Doppler folds modulo $2 V_{\text{max}}$:

$$v_{\text{meas}} = \text{mod}(v_{\text{true}} + V_{\text{max}}, 2 V_{\text{max}}) - V_{\text{max}}$$

To resolve velocities up to multiple folds of $V_{\text{max}}$, CHIRP transmits alternating frames with staggered chirp periods $T_1$ and $T_2$ having unambiguous velocities $V_{\text{max},1}$ and $V_{\text{max},2}$. The expanded unambiguous velocity is:

$$V_{\text{expanded}} = \frac{V_{\text{max},1} \cdot V_{\text{max},2}}{|V_{\text{max},1} - V_{\text{max},2}|}$$

The true velocity is recovered by solving the integer ambiguity hypothesis pair:

$$\hat{v} = v_1 + k_1 (2 V_{\text{max},1}) = v_2 + k_2 (2 V_{\text{max},2})$$

## 5. Moving Target Indication (MTI) Clutter Rejection

Stationary ground clutter and zero-Doppler reflections from guardrails and road surfaces are filtered using slow-time finite impulse response high-pass filters:

- **Single Delay Canceler**: $y[m] = x[m] - x[m-1]$, frequency response $H(z) = 1 - z^{-1}$.
- **Double Delay Canceler**: $y[m] = x[m] - 2x[m-1] + x[m-2]$, frequency response $H(z) = (1 - z^{-1})^2$.
- **Polynomial Fitting**: Estimates low-frequency clutter trends via least-squares polynomial subtraction across slow-time.
