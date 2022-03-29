function test_microdoppler()
    % TEST_MICRODOPPLER Unit tests for micro-Doppler time-frequency analysis.
    
    fprintf('Running test_microdoppler...\n');
    
    fs = 1000.0;
    t = 0:(1/fs):0.5;
    % Two-component signal: 100 Hz carrier + 200 Hz tone
    sig = exp(1j * 2*pi * 100 * t) + 0.5 * exp(1j * 2*pi * 200 * t);
    
    % 1. STFT Spectrogram
    stft = chirp.microdoppler.StftAnalyzer('WindowLength', 64, 'HopSize', 8, 'Nfft', 128);
    [pwrDb, tAx, fAx] = stft.compute(sig, fs);
    assert(size(pwrDb, 1) == 128, 'STFT frequency bins mismatch');
    assert(length(tAx) == size(pwrDb, 2), 'STFT time bins mismatch');
    
    % 2. Micro-Doppler Feature extraction
    feats = chirp.microdoppler.MicroDopplerFeatures.extract(10.^(pwrDb / 10.0), tAx, fAx);
    assert(feats.meanBandwidthHz > 50.0, 'Extracted bandwidth must capture tone separation');
    assert(feats.spectralEntropy > 0.0 && feats.spectralEntropy < 1.0, 'Entropy out of range [0, 1]');
    
    % 3. Cadence extraction from synthetic periodic envelope
    tEnv = 0:0.01:2.0; % 100 Hz envelope sampling
    envPeriodic = 2.0 * sin(2*pi * 2.0 * tEnv); % 2 Hz cadence
    [cadence, conf] = chirp.microdoppler.CadenceExtractor.extractFromEnvelope(envPeriodic, 100.0);
    fprintf('  Extracted cadence: %.2f Hz (True: 2.00 Hz), Confidence: %.2f\n', cadence, conf);
    assert(abs(cadence - 2.0) < 0.25, 'Cadence extraction error');
    assert(conf > 0.5, 'Cadence confidence too low');
    
    % 4. Kinematic Classifier
    lambda = 0.00389;
    pedFeats = struct('meanBandwidthHz', 700.0, 'maxDopplerHz', 750.0, 'cadenceHz', 1.8, ...
                      'torsoVelocityHz', 720.0, 'spectralEntropy', 0.65);
    [cName, conf, scores] = chirp.microdoppler.KinematicClassifier.classify(pedFeats, lambda);
    assert(strcmp(cName, 'pedestrian'), 'Failed to classify pedestrian signature');
    assert(conf > 0.6, 'Pedestrian classification confidence too low');
    
    fprintf('  PASSED: test_microdoppler\n');
end
