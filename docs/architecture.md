# System Architecture and Design Specification

## Overview

CHIRP is a modular, high-performance Radar Signal Processing and Synthetic Aperture Radar (SAR) perception engine written in MATLAB with native C MEX acceleration kernels. The architecture is engineered to provide an end-to-end processing pipeline transforming raw analog-to-digital converter (ADC) electromagnetic reflections into track-level 4D situational awareness and focused sub-meter radar imagery.

## System Data Flow Pipeline

The end-to-end data processing pipeline flows through eight sequential stages:

```mermaid
flowchart TD
    subgraph S1["1. Raw Data Ingestion & Synthesis"]
        ADC["Raw Baseband ADC Data / DCA1000 UDP"] --> Dechirp["De-Chirping & IQ Demodulation"]
        Sim["Waveform Synthesis & Multipath Sim"] --> Dechirp
    end

    subgraph S2["2. Core Radar DSP"]
        Dechirp --> Cube["4D Radar Cube [N_adc x N_chirp x N_virt]"]
        Cube --> RFFT["Range-FFT + Windowing (Hann/Blackman)"]
        RFFT --> DFFT["Doppler-FFT (Velocity Disambiguation)"]
        DFFT --> NCI["Non-Coherent Channel Integration"]
    end

    subgraph S3["3. Adaptive Detection (CFAR)"]
        NCI --> RDM["Range-Doppler Map (RDM)"]
        RDM --> CFAR["2D CFAR Detector (CA / GO / SO / OS)"]
        CFAR --> Peaks["Detections Mask & Sub-Cell Interp"]
    end

    subgraph S4["4. Direction of Arrival (AoA)"]
        Peaks --> AoA["Spatial Covariance R = E[x x^H]"]
        AoA --> SuperRes["Super-Resolution (MUSIC / Root-MUSIC / ESPRIT)"]
        SuperRes --> Spherical["3D Spherical Coordinates (r, az, el, v)"]
    end

    subgraph S5["5. 4D Point Cloud Perception"]
        Spherical --> Cloud["PointCloud4D (X, Y, Z, V, SNR)"]
        Cloud --> DBSCAN["Velocity-Aware DBSCAN Clustering"]
        DBSCAN --> BBox["Oriented Bounding Box Estimation"]
    end

    subgraph S6["6. State Estimation & Multi-Target Tracking"]
        BBox --> GNN["Global Nearest Neighbor / JPDA Matcher"]
        GNN --> IMM["IMM Filter (CV + CA + CT Models)"]
        IMM --> Tracks["Track Lifecycle Manager (Tentative/Confirmed/Coasted)"]
    end

    subgraph S7["7. Imaging & Micro-Doppler Auxiliary"]
        DFFT --> STFT["STFT / Synchrosqueezing Spectrogram"]
        STFT --> Classifier["Kinematic Target Classifier"]
        Cube --> SAR["Time-Domain Backprojection SAR Engine"]
        SAR --> FocusedSAR["Focused 2D Ground Reflectivity Map"]
    end

    subgraph S8["8. Visualization & Export"]
        Tracks --> GUI["ANSI Terminal & Vector SVG Renderers"]
        FocusedSAR --> GUI
        Classifier --> GUI
        Cloud --> IO["PLY / LAS / GeoJSON / RCUBE Writers"]
    end
```

## Package Architecture

The codebase is organized into modular packages under the `+chirp` namespace:

| Package | Description | Core Classes |
| :--- | :--- | :--- |
| `+chirp/+utils` | Physical constants, fast bitwise math, matrix utilities, CRC32 | `PhysicalConstants`, `FastMath`, `MatrixUtils`, `Crc32`, `StatUtils` |
| `+chirp/+waveform` | Waveform definitions, ambiguity functions, phase modulation | `FmcwChirp`, `AmbiguityFunction`, `PhaseCodeModulator`, `WaveformGenerator` |
| `+chirp/+sim` | Multi-channel physical radar and RF channel simulation | `RadarConfig`, `TargetModel`, `SwerlingModel`, `MimoArray`, `RadarSimulator` |
| `+chirp/+dsp` | 2D Range-Doppler processing, MTI filtering, phase calibration | `RangeProcessor`, `DopplerProcessor`, `RangeDopplerMap`, `MtiFilter`, `DopplerAmbiguityResolver` |
| `+chirp/+aoa` | Direction of Arrival beamforming and super-resolution | `BartlettBeamformer`, `CaponBeamformer`, `MusicEstimator`, `RootMusicEstimator`, `EspritEstimator` |
| `+chirp/+cfar` | Adaptive thresholding detectors and sub-cell interpolation | `CaCfar`, `GoCfar`, `SoCfar`, `OsCfar`, `TmCfar`, `ThresholdCalculator` |
| `+chirp/+pointcloud`| 4D Point cloud representations, spatial indexing, clustering | `PointCloud4D`, `SpatialKdTree`, `VelocityAwareDbscan`, `BoundingBoxEstimator` |
| `+chirp/+tracking` | Multi-target state estimation, Kalman filters, data association| `ExtendedKalmanFilter`, `UnscentedKalmanFilter`, `ImmFilter`, `GnnTracker`, `JpdaTracker` |
| `+chirp/+microdoppler`| Time-frequency transforms, gait cadence, target classification | `StftAnalyzer`, `SynchrosqueezingTransform`, `MicroDopplerFeatures`, `KinematicClassifier` |
| `+chirp/+sar` | Synthetic aperture radar image reconstruction algorithms | `SarGeometry`, `BackprojectionEngine`, `RangeDopplerAlgorithm`, `PolarFormatAlgorithm` |
| `+chirp/+io` | Hardware capture readers, standard format exporters | `Dca1000Reader`, `RadarCubeCodec`, `PlyWriter`, `LasWriter`, `GeoJsonWriter` |
| `+chirp/+viz` | Terminal ANSI scopes, vector SVG visualizers | `TerminalHeatmap`, `TerminalPointCloud`, `SvgPointCloudRenderer`, `PolarPlotRenderer` |

## Memory Layout and Tensor Representations

Radar cubes are organized in column-major continuous memory buffers:
- **Fast-Time (ADC Samples)**: Dimension 1 ($N_{\text{adc}}$ samples along the frequency chirp sweep).
- **Slow-Time (Chirp Index)**: Dimension 2 ($N_{\text{chirps}}$ pulses per frame).
- **Spatial (Virtual Channel)**: Dimension 3 ($N_{\text{virt}} = N_{\text{TX}} \times N_{\text{RX}}$ virtual co-array phase centers).

This memory layout ensures that the first processing step (1D Range-FFT) operates on contiguous memory columns, maximizing CPU L1/L2 cache hit rates and enabling SIMD auto-vectorization.
