# ADR-0001: MIMO Virtual Array Geometry and Co-Array Synthesis

## Status
Accepted

## Context
Automotive radar sensors operating in the 77 to 81 GHz millimeter-wave band face severe physical aperture constraints. A physical phased array with half-wavelength spacing ($\lambda/2 \approx 1.95\text{ mm}$) requires dozens of physical receiver chains to achieve the narrow beamwidths ($< 2^\circ$) required for long-range autonomous driving on multi-lane highways. Physical receiver channels scale cost, die area, power consumption, and thermal dissipation linearly.

Multiple-Input Multiple-Output (MIMO) radar technology overcomes this hardware limitation by transmitting orthogonal waveforms from multiple spatially separated transmit antennas ($N_{\text{TX}}$) and receiving reflected echoes across a sparse receive array ($N_{\text{RX}}$). The effective aperture is synthesized through spatial convolution of the transmit and receive phase centers, forming an equivalent virtual single-input multiple-output (SIMO) array with $N_{\text{virt}} = N_{\text{TX}} \times N_{\text{RX}}$ phase centers.

## Decision
We adopt a Time-Division Multiplexed (TDM) and Phase-Coded MIMO architecture with a non-redundant spatial geometry:
1. Physical receive antennas are arranged as a uniform linear array (ULA) with element spacing $d_{\text{RX}} = \lambda / 2$.
2. Physical transmit antennas are spaced at large multiples of the receive array aperture, specifically $d_{\text{TX}} = N_{\text{RX}} \cdot d_{\text{RX}} = 4 \cdot (\lambda / 2) = 2\lambda$.
3. Virtual channel spatial convolution yields a continuous, filled, non-redundant 12-element virtual ULA spanning total aperture length $D_{\text{virt}} = (N_{\text{TX}} \cdot N_{\text{RX}} - 1) \cdot (\lambda / 2) = 11 \cdot (\lambda / 2) \approx 21.4\text{ mm}$.
4. Orthogonal multiplexing is achieved primarily via TDM chirp interleaving, supplemented by binary phase-modulated Hadamard codes for Doppler division multiplexing (DDMA).

## Consequences
### Positive
- Triple angular resolution without increasing receive RF front-end complexity or ADC channels (12 virtual channels realized using only 4 physical ADCs).
- Linear steering vector formulation remains uniform, enabling closed-form polynomial root-finding via Root-MUSIC and rotational shift invariance via ESPRIT.
- Zero spatial grating lobes within the primary field of view ($\pm 60^\circ$) due to continuous half-wavelength virtual spacing.

### Negative
- TDM interleaving reduces the maximum unambiguous Doppler velocity by factor $N_{\text{TX}}$ ($V_{\text{max,TDM}} = V_{\text{max}} / N_{\text{TX}}$).
- Motion-induced phase errors between successive TX chirps introduce virtual element phase perturbations, requiring Doppler compensation before spatial beamforming.
- DDMA requires slow-time phase decoding matrices, slightly increasing DSP memory footprint.
