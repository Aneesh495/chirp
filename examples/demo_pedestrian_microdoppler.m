function demo_pedestrian_microdoppler()
    % DEMO_PEDESTRIAN_MICRODOPPLER Biomechanical human gait micro-Doppler analysis.
    %
    % Simulates walking human kinematics with swinging limbs and torso reflections,
    % computing high-resolution STFT spectrogram, cadence estimation, and kinematic
    % feature classification, exporting vector SVG time-frequency visualization.
    
    fprintf('=================================================================\n');
    fprintf('  CHIRP DEMO: PEDESTRIAN MICRO-DOPPLER GAIT CLASSIFICATION      \n');
    fprintf('=================================================================\n\n');
    
    thisDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(thisDir);
    addpath(fullfile(projectRoot, 'src'));
    addpath(fullfile(projectRoot, 'bin', 'mex'));
    
    % 1. Radar Kinematics and Simulation Configuration
    fc = 77.0e9;
    lambda = chirp.utils.PhysicalConstants.freq2wavelength(fc);
    prf = 2000.0; % 2 kHz PRF for slow-time sampling
    duration = 1.5; % 1.5 seconds of gait observation
    t = 0:(1.0 / prf):duration;
    
    vTorso = 1.35; % Walking speed (m/s)
    strideCadence = 1.85; % Gait cadence (Hz, ~110 steps/min)
    
    fprintf('Kinematic Walking Parameters:\n');
    fprintf('  Carrier Frequency : %.2f GHz (Wavelength = %.2f mm)\n', fc / 1e9, lambda * 1e3);
    fprintf('  Walking Speed     : %.2f m/s (Torso Doppler = %.1f Hz)\n', vTorso, 2.0 * vTorso / lambda);
    fprintf('  Gait Cadence      : %.2f Hz (Stride period = %.2f s)\n\n', strideCadence, 1.0 / strideCadence);
    
    % 2. Synthesize Multi-Component Biomechanical Echo
    % Torso return: constant baseline Doppler with slight oscillation
    fTorso = (2.0 * vTorso / lambda) * (1.0 + 0.05 * sin(2.0 * pi * 2.0 * strideCadence * t));
    phaseTorso = 2.0 * pi * cumsum(fTorso) / prf;
    sigTorso = 1.0 * exp(1j * phaseTorso);
    
    % Left and Right Arms: anti-phase swings with velocity peak at 2.2 m/s
    fArms = (2.0 * 0.8 / lambda) * sin(2.0 * pi * strideCadence * t);
    sigArms = 0.3 * exp(1j * (phaseTorso + 2.0 * pi * cumsum(fArms) / prf)) + ...
              0.3 * exp(1j * (phaseTorso - 2.0 * pi * cumsum(fArms) / prf));
              
    % Left and Right Legs: intermittent forward thrusts up to 2.8 m/s
    fLegs = (2.0 * 1.4 / lambda) * max(0.0, sin(2.0 * pi * strideCadence * t)).^2;
    sigLegs = 0.4 * exp(1j * (phaseTorso + 2.0 * pi * cumsum(fLegs) / prf));
    
    % Combined echo with thermal background noise
    sigTotal = sigTorso + sigArms + sigLegs + 0.05 * (randn(size(t)) + 1j * randn(size(t)));
    
    % 3. Compute High-Resolution STFT Spectrogram
    fprintf('Computing Short-Time Fourier Transform (STFT) Spectrogram...\n');
    stft = chirp.microdoppler.StftAnalyzer('WindowLength', 128, 'HopSize', 16, 'Nfft', 256);
    [pwrDb, tAx, fAx] = stft.compute(sigTotal, prf);
    
    pwrLin = 10.0.^(pwrDb / 10.0);
    
    % 4. Extract Biomechanical Features
    fprintf('Extracting Micro-Doppler Spectral Features...\n');
    feats = chirp.microdoppler.MicroDopplerFeatures.extract(pwrLin, tAx, fAx);
    
    fprintf('  Torso Doppler Freq : %6.1f Hz (Estimated Speed = %.2f m/s)\n', ...
        feats.torsoVelocityHz, feats.torsoVelocityHz * lambda / 2.0);
    fprintf('  Micro-Doppler BW   : %6.1f Hz (Velocity Span = %.2f m/s)\n', ...
        feats.meanBandwidthHz, feats.meanBandwidthHz * lambda / 2.0);
    fprintf('  Dominant Cadence   : %6.2f Hz (True: %.2f Hz)\n', ...
        feats.cadenceHz, strideCadence);
    fprintf('  Spectral Entropy   : %6.3f\n', feats.spectralEntropy);
    
    % 5. Kinematic Target Classification
    [className, confidence, classScores] = chirp.microdoppler.KinematicClassifier.classify(feats, lambda);
    
    fprintf('\n--- Kinematic Classification Result ---\n');
    fprintf('  Predicted Class    : %s\n', upper(className));
    fprintf('  Confidence Metric  : %5.1f %%\n', confidence * 100.0);
    fprintf('  Score Breakdown    : Pedestrian = %.2f | Cyclist = %.2f | Vehicle = %.2f | Drone = %.2f\n', ...
        classScores.pedestrian, classScores.cyclist, classScores.vehicle, classScores.drone);
        
    % 6. Render Vector SVG Spectrogram Visualization
    outDir = fullfile(projectRoot, 'output');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    
    svgPath = fullfile(outDir, 'pedestrian_microdoppler.svg');
    chirp.viz.MicroDopplerRenderer.renderToFile(svgPath, pwrDb, tAx, fAx, feats);
    fprintf('\nExported vector SVG spectrogram visualization to: %s\n', svgPath);
    
    fprintf('\nPedestrian Micro-Doppler demo completed successfully.\n');
end
