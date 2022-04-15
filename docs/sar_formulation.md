# Synthetic Aperture Radar (SAR) Mathematical Formulation

## 1. Principles of Synthetic Aperture Imaging

In conventional real-aperture radar, the angular resolution along cross-range is limited by the diffraction limit of the physical antenna aperture $D_a$:

$$\theta_{3\text{dB}} \approx \frac{\lambda}{D_a} \implies \Delta x = R \cdot \theta_{3\text{dB}} \approx \frac{\lambda R}{D_a}$$

At long stand-off distances ($R = 10\text{ km}$), achieving a 1-meter cross-range resolution at $X$-band ($\lambda = 0.03\text{ m}$) would require a physical antenna aperture of $D_a = 300\text{ meters}$, which is impossible to fly on an aircraft or satellite.

Synthetic Aperture Radar (SAR) synthesizes an equivalent virtual antenna array of length $L_{\text{syn}}$ by transmitting coherent pulses as the radar moves along a trajectory. The synthesized beamwidth is:

$$\theta_{\text{syn}} = \frac{\lambda}{2 L_{\text{syn}}}$$

where the factor of 2 accounts for the two-way round-trip phase progression. In stripmap mode, the maximum integration baseline is limited by the physical footprint: $L_{\text{syn}} = R \cdot \frac{\lambda}{D_a}$. Substituting into cross-range resolution:

$$\Delta x = R \cdot \theta_{\text{syn}} = R \cdot \frac{\lambda}{2 L_{\text{syn}}} = R \cdot \frac{\lambda}{2 \left( R \frac{\lambda}{D_a} \right)} = \frac{D_a}{2}$$

The physical cross-range resolution of a focused SAR system is completely independent of platform range and carrier wavelength, depending solely on half the physical antenna dimension $D_a / 2$.

## 2. Time-Domain Backprojection Algorithm (BPA)

Time-Domain Backprojection is an exact, geometry-preserving tomographic image reconstruction algorithm. Let the platform position at slow-time pulse index $p$ be $\mathbf{r}_p = [x_p, y_p, z_p]^T$ and let an arbitrary ground grid pixel coordinate be $\mathbf{r}_{\text{pix}} = [x_g, y_g, z_g]^T$.

The exact instantaneous Euclidean slant range is:

$$R(p, \mathbf{r}_{\text{pix}}) = \|\mathbf{r}_p - \mathbf{r}_{\text{pix}}\|_2 = \sqrt{(x_p - x_g)^2 + (y_p - y_g)^2 + (z_p - z_g)^2}$$

Let $s_{\text{rc}}(r, p)$ denote the range-compressed baseband profile at pulse $p$. The focused complex reflectivity image $I(\mathbf{r}_{\text{pix}})$ is reconstructed by evaluating the coherent matched-filter backprojection integral:

$$I(\mathbf{r}_{\text{pix}}) = \sum_{p=1}^{N_{\text{pulses}}} s_{\text{rc}}\left( R(p, \mathbf{r}_{\text{pix}}), p \right) \cdot \exp\left( +j \frac{4\pi}{\lambda} R(p, \mathbf{r}_{\text{pix}}) \right) \cdot w(p)$$

where:
- $s_{\text{rc}}\left( R(p, \mathbf{r}_{\text{pix}}), p \right)$ is evaluated at continuous slant range $R$ via sub-pixel linear or cubic interpolation.
- $\exp\left( +j \frac{4\pi}{\lambda} R \right)$ removes the geometric phase delay, aligning all pulse contributions coherently in the complex plane.
- $w(p)$ is a slow-time apodization window (e.g. Taylor or Kaiser) controlling azimuth sidelobes.

### Computational Complexity
For a ground reconstruction grid of $N_X \times N_Y$ pixels and $N_P$ slow-time pulses, the total operational count is $O(N_X \cdot N_Y \cdot N_P)$. While computationally intensive, BPA requires zero approximations, accommodates arbitrary non-linear flight paths, and parallelizes across image pixels.

## 3. Range-Doppler Algorithm (RDA) with RCMC

The Range-Doppler Algorithm operates in the hybrid Range-Doppler $(r, f_\eta)$ domain to correct Range Cell Migration (RCM) via 2D separable filtering.

### Slant Range Taylor Expansion
For a straight-and-level flight along the $X$-axis with constant velocity $v$:

$$R(\eta) = \sqrt{R_0^2 + (v \eta)^2} \approx R_0 + \frac{v^2 \eta^2}{2 R_0}$$

The range migration across slow-time aperture time $\eta$ is:

$$\Delta R_{\text{rcm}}(\eta) = R(\eta) - R_0 \approx \frac{v^2 \eta^2}{2 R_0}$$

### Range Cell Migration Correction (RCMC)
Transforming slow-time $\eta$ to Doppler frequency $f_\eta$ via azimuth FFT:

$$\Delta R_{\text{rcm}}(f_\eta) = \frac{\lambda^2 R_0 f_\eta^2}{8 v^2}$$

In the Range-Doppler domain $(r, f_\eta)$, each range line is shifted along fast-time range by $-\Delta R_{\text{rcm}}(f_\eta)$ via sinc interpolation, straightening curved hyperbolic trajectories into vertical columns.

### Azimuth Matched Filtering
After RCMC, azimuth compression is applied in the Doppler frequency domain by multiplying with the matched filter transfer function:

$$H_{\text{az}}(f_\eta) = \exp\left( -j \pi \frac{\lambda R_0}{2 v^2} f_\eta^2 \right)$$

An inverse azimuth FFT yields the focused 2D image $I(x, y)$.

## 4. Polar Format Algorithm (PFA)

In Spotlight SAR, the antenna continuously steers toward a fixed ground scene center $\mathbf{r}_c = [0, R_0, 0]^T$. The received signal samples trace a polar grid in spatial frequency space $(k_X, k_Y)$:

$$k_X = \frac{4\pi}{\lambda} \sin(\theta), \quad k_Y = \frac{4\pi}{\lambda} \cos(\theta)$$

where $\theta(p)$ is the instantaneous platform squint angle. PFA resamples these polar spatial frequencies onto a uniform 2D Cartesian wavenumber grid $(k_X, k_Y)$ using 1D separable chirp z-transforms or sinc interpolators, followed by a standard 2D inverse FFT to focus the spotlight scene.

## 5. Phase Gradient Autofocus (PGA)

Uncompensated platform flight turbulence, wind gusts, and GPS/INS drift introduce unknown phase errors $\phi_e(p)$ across azimuth slow-time. Phase Gradient Autofocus (PGA) estimates and cancels arbitrary higher-order phase errors without parametric polynomial assumptions:

1. **Circular Shift**: Shifts dominant isolated point-like scatterers in each range bin to the zero-Doppler axis.
2. **Windowing**: Applies an adaptive rectangular window to eliminate surrounding clutter and noise.
3. **Phase Gradient Estimation**: Computes the maximum-likelihood derivative of phase error from conjugate products:

$$\dot{\phi}_e(p) = \text{arg}\left( \sum_{m} y[m, p] \cdot y^*[m, p-1] \right)$$

4. **Integration & Iteration**: Integrates the gradient to recover $\phi_e(p)$, de-rotates the phase history, and iterates until the residual phase error variance converges below a threshold.
