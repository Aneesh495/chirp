function test_end_to_end_perception()
    % TEST_END_TO_END_PERCEPTION Integration test verifying full radar perception stack.
    %
    % Executes complete pipeline:
    %   Electromagnetic Waveform -> 4D Radar Cube -> Range-Doppler DSP ->
    %   2D CFAR -> Super-Resolution AoA -> 4D Point Cloud -> DBSCAN -> GNN Tracking
    
    fprintf('Running test_end_to_end_perception...\n');
    
    % Step 1: System and waveform configuration
    cfg = chirp.sim.RadarConfig('NumChirps', 32, 'TxPowerDbm', 15.0);
    chirpObj = chirp.waveform.WaveformGenerator.automotiveLRR();
    chirpObj.NumAdcSamples = 256;
    mimo = chirp.sim.MimoArray.automotive3Tx4Rx(chirpObj.Wavelength);
    
    sim = chirp.sim.RadarSimulator(cfg, chirpObj, mimo);
    sim.EnableNoise = true;
    
    % Step 2: Define active physical obstacles
    % Target 1: Vehicle cruising ahead
    tgt1 = chirp.sim.TargetModel(...
        'Position', [1.5; 35.0; 0.0], ...
        'Velocity', [0.0; 12.0; 0.0], ...
        'Rcs', 10.0, ...
        'TargetId', 1 ...
    );
    
    % Target 2: Crossing obstacle
    tgt2 = chirp.sim.TargetModel(...
        'Position', [-3.0; 20.0; 0.0], ...
        'Velocity', [0.0; -4.0; 0.0], ...
        'Rcs', 3.0, ...
        'TargetId', 2 ...
    );
    
    dt = 0.05; % 20 Hz frame rate
    tracker = chirp.tracking.GnnTracker();
    rp = chirp.dsp.RangeProcessor('WindowType', 'blackman-harris');
    dp = chirp.dsp.DopplerProcessor('WindowType', 'chebyshev');
    ca = chirp.cfar.CaCfar(1e-4, 8, 2, 4, 1);
    music = chirp.aoa.MusicEstimator(mimo, chirpObj.Wavelength, true);
    
    % Execute 6 consecutive tracking frames
    for frame = 1:6
        time = (frame - 1) * dt;
        
        % Synthesize radar cube for current frame
        cube = sim.simulateFrame([tgt1, tgt2]);
        
        % 1D Range-FFT
        [rangeProf, rAxis] = rp.process(cube, chirpObj);
        
        % 2D Doppler-FFT
        [rdmCube, vAxis] = dp.process(rangeProf, chirpObj);
        
        % Non-coherent integration Range-Doppler map
        rdm = chirp.dsp.RangeDopplerMap(rdmCube, rAxis, vAxis);
        pwrDb = rdm.nonCoherentIntegration();
        pwrLin = 10.^(pwrDb / 10.0);
        
        % 2D CFAR detection
        [mask, ~, ~] = ca.detect(pwrLin);
        [rIdx, vIdx] = find(mask);
        
        % Estimate AoA for each detection
        nDet = length(rIdx);
        detectionsSpherical = zeros(4, nDet);
        
        for k = 1:nDet
            r = rAxis(rIdx(k));
            v = vAxis(vIdx(k));
            snapshot = squeeze(rdmCube(rIdx(k), vIdx(k), :));
            
            [specM, azGrid] = music.estimateAzimuth(snapshot, 1, -30:1:30, 0);
            [~, pkAz] = max(specM);
            azDeg = azGrid(pkAz);
            azRad = azDeg * (pi / 180.0);
            
            detectionsSpherical(:, k) = [r; azRad; 0.0; v];
        end
        
        % Multi-target tracking step
        activeTracks = tracker.step(detectionsSpherical, time, dt);
        
        % Propagate targets for next frame
        tgt1 = tgt1.propagate(dt);
        tgt2 = tgt2.propagate(dt);
    end
    
    fprintf('  End-to-End Tracking Results (%d active tracks):\n', length(activeTracks));
    assert(~isempty(activeTracks), 'Zero tracks recovered from end-to-end radar cube');
    
    for t = 1:length(activeTracks)
        trk = activeTracks(t);
        fprintf('    Track #%d | Status: %s | Hits: %d | Pos: [%.2f, %.2f, %.2f] m | Vel: [%.2f, %.2f, %.2f] m/s\n', ...
            trk.TrackId, trk.Status, trk.Hits, ...
            trk.State(1), trk.State(2), trk.State(3), ...
            trk.State(4), trk.State(5), trk.State(6));
    end
    
    fprintf('  PASSED: test_end_to_end_perception\n');
end
