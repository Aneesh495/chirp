function test_waveform()
    % TEST_WAVEFORM Unit tests for radar waveform modeling and synthesis.
    
    fprintf('Running test_waveform...\n');
    
    % 1. FMCW Chirp parameters
    chirpObj = chirp.waveform.WaveformGenerator.automotiveLRR();
    assert(abs(chirpObj.Slope - 1.0e13) < 1e6, 'Slope mismatch');
    c = chirp.utils.PhysicalConstants.C;
    expectedRes = c / (2.0 * chirpObj.SweepBandwidth);
    assert(abs(chirpObj.RangeResolution - expectedRes) < 1e-12, 'Range resolution mismatch');
    assert(chirpObj.MaxUnambiguousRange > 100.0, 'Max unambiguous range too small');
    
    % 2. De-chirped beat signal
    rTarget = 35.0;
    vTarget = 10.0;
    [beatSig, fb] = chirpObj.dechirpTarget(rTarget, vTarget, 1.0);
    assert(length(beatSig) == chirpObj.NumAdcSamples, 'Beat signal length mismatch');
    
    % Analytical beat frequency: fb = 2*S*R/c + 2*fc*v/c
    expectedFb = (2.0 * chirpObj.Slope * rTarget) / c + (2.0 * chirpObj.StartFrequency * vTarget) / c;
    assert(abs(fb - expectedFb) < 1e-3, 'Analytical beat frequency calculation mismatch');
    
    % 3. Barker code autocorrelation and Peak Sidelobe Level (PSL)
    b13 = chirp.waveform.PhaseCodeModulator.barker(13);
    psl13 = chirp.waveform.PhaseCodeModulator.peak_sidelobe_ratio_db(b13);
    % Barker-13 theoretical PSL: 20*log10(1/13) = -22.28 dB
    assert(abs(psl13 - (-22.28)) < 0.1, 'Barker-13 PSL mismatch');
    
    % 4. Frank code length
    frank4 = chirp.waveform.PhaseCodeModulator.frank(4);
    assert(length(frank4) == 16, 'Frank code length mismatch');
    
    % 5. Zadoff-Chu sequence CAZAC property (constant modulus)
    zc = chirp.waveform.PhaseCodeModulator.zadoff_chu(1, 31);
    assert(length(zc) == 31, 'Zadoff-Chu length mismatch');
    modulusD = abs(abs(zc) - 1.0);
    assert(max(modulusD) < 1e-12, 'Zadoff-Chu not constant amplitude');
    
    % 6. Stepped-Frequency Continuous-Wave (SFCW)
    sfcw = chirp.waveform.FrequencySteppedChirp('NumSteps', 64, 'FrequencyStep', 10.0e6);
    rx = sfcw.simulateReturns(12.0, 1.0);
    [prof, rAx] = sfcw.reconstructProfile(rx);
    [~, maxIdx] = max(abs(prof));
    assert(abs(rAx(maxIdx) - 12.0) < 0.5, 'SFCW target reconstruction error');
    
    fprintf('  PASSED: test_waveform\n');
end
