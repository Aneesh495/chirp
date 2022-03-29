function test_dsp()
    % TEST_DSP Unit tests for core radar signal processing algorithms.
    
    fprintf('Running test_dsp...\n');
    
    % 1. Window metrics
    w = chirp.dsp.WindowFunction.create('blackman-harris', 128);
    metrics = chirp.dsp.WindowFunction.evaluateMetrics(w);
    assert(metrics.enbw > 1.8 && metrics.enbw < 2.2, 'Blackman-Harris ENBW out of nominal range');
    assert(metrics.scallopingDb > 0.5, 'Scalloping loss must be positive');
    
    % 2. 1D Range-FFT
    chirpObj = chirp.waveform.WaveformGenerator.automotiveLRR();
    chirpObj.NumAdcSamples = 128;
    [beatSig, fb] = chirpObj.dechirpTarget(24.0, 0.0, 1.0);
    cube = reshape(beatSig, [128, 1, 1]);
    
    rp = chirp.dsp.RangeProcessor('Nfft', 256);
    [prof, rAxis] = rp.process(cube, chirpObj);
    [~, maxBin] = max(abs(prof(:, 1, 1)));
    fprintf('  Range-FFT peak range: %.2f m (True: 24.0 m)\n', rAxis(maxBin));
    assert(abs(rAxis(maxBin) - 24.0) < 0.6, 'Range-FFT peak location error');
    
    % 3. MTI delay-line cancelers
    % Stationary DC signal across chirps (0 Doppler)
    dcData = ones(64, 32, 2);
    mti2 = chirp.dsp.MtiFilter.twoPulseCanceler(dcData);
    % After first sample, all differences must be exactly 0
    diffSlice = mti2(:, 2:end, :);
    assert(max(abs(diffSlice(:))) < 1e-12, 'Two-pulse canceler failed to reject stationary DC clutter');
    
    % 4. IQ Calibration
    rawI = 2.0 * cos(linspace(0, 10*pi, 1000)) + 0.5; % DC offset = 0.5, Gain = 2.0
    rawQ = 1.0 * sin(linspace(0, 10*pi, 1000) + 0.1) - 0.3; % DC = -0.3, Phase err = 0.1 rad
    rawIq = rawI + 1j * rawQ;
    
    calib = chirp.dsp.IqCalibrator();
    calib = calib.estimateBlind(rawIq);
    assert(abs(calib.DcOffsetI - 0.5) < 0.1, 'Estimated DC I mismatch');
    assert(abs(calib.DcOffsetQ - (-0.3)) < 0.1, 'Estimated DC Q mismatch');
    assert(abs(calib.GainMismatch - 2.0) < 0.2, 'Estimated Gain mismatch');
    
    % 5. Doppler Ambiguity Resolution
    lambda = 0.00389;
    resolver = chirp.dsp.DopplerAmbiguityResolver(lambda, 20000, 25000);
    % True velocity = 28.0 m/s (exceeds Vmax1 = 19.45 m/s)
    trueV = 28.0;
    v1_aliased = mod(trueV + resolver.Vmax1, 2 * resolver.Vmax1) - resolver.Vmax1;
    v2_aliased = mod(trueV + resolver.Vmax2, 2 * resolver.Vmax2) - resolver.Vmax2;
    resolvedV = resolver.resolveChineseRemainder(v1_aliased, v2_aliased, 4);
    assert(abs(resolvedV - trueV) < 0.1, 'Chinese Remainder Doppler de-aliasing failed');
    
    % 6. Regression test for RangeDopplerMap findPeaks strict peak test
    pwrFlat = 1e-4 * ones(32, 32);
    pwrFlat(15, 15) = 1.0; % Single isolated peak: 0 dB vs -80 dB
    rdm = chirp.dsp.RangeDopplerMap(pwrFlat, 1:32, 1:32);
    [pks, rEst, vEst] = rdm.findPeaks([], 2);
    assert(length(pks) == 1, 'Bug regression: findPeaks must return exactly 1 peak on plateau');
    assert(rEst(1) == 15 && vEst(1) == 15, 'Peak coordinate mismatch');
    
    fprintf('  PASSED: test_dsp\n');
end
