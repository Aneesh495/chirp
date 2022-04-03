function demo_4d_automotive_radar()
    % DEMO_4D_AUTOMOTIVE_RADAR End-to-end 4D automotive radar perception demonstration.
    %
    % Simulates a 77 GHz FMCW MIMO radar on a vehicle navigating a highway with
    % multiple dynamic targets, executing the complete signal chain:
    %   1. FMCW chirp waveform synthesis and MIMO channel propagation
    %   2. 2D Range-Doppler processing (Windowed Range FFT + Doppler FFT)
    %   3. 2D CFAR detection (Cell-Averaging with guard bands)
    %   4. Super-resolution Azimuth/Elevation angle estimation (MUSIC)
    %   5. Doppler-penalized spatial clustering (Velocity-Aware DBSCAN)
    %   6. Track filtering and maintenance (GNN Tracker)
    %   7. ASCII terminal dashboard and vector SVG point cloud rendering
    
    fprintf('=================================================================\n');
    fprintf('  CHIRP DEMO: 4D AUTOMOTIVE RADAR PERCEPTION PIPELINE            \n');
    fprintf('=================================================================\n\n');
    
    thisDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(thisDir);
    addpath(fullfile(projectRoot, 'src'));
    addpath(fullfile(projectRoot, 'bin', 'mex'));
    
    % 1. System and Waveform Configuration
    cfg = chirp.sim.RadarConfig('NumChirps', 32, 'TxPowerDbm', 15.0);
    chirpObj = chirp.waveform.WaveformGenerator.automotiveLRR();
    chirpObj.NumAdcSamples = 256;
    mimo = chirp.sim.MimoArray.automotive3Tx4Rx(chirpObj.Wavelength);
    
    fprintf('Radar Configuration:\n');
    fprintf('  Carrier Frequency : %.2f GHz (Wavelength = %.2f mm)\n', chirpObj.StartFrequency / 1e9, chirpObj.Wavelength * 1e3);
    fprintf('  Sweep Bandwidth   : %.1f MHz (Range Res = %.2f m)\n', chirpObj.SweepBandwidth / 1e6, chirpObj.RangeResolution);
    apertureMm = (max(mimo.VirtualPositions(1, :)) - min(mimo.VirtualPositions(1, :))) * 1e3;
    fprintf('  Virtual Channels  : %d elements (Aperture = %.1f mm)\n', mimo.NumVirtual, apertureMm);
    vMax = chirpObj.Wavelength / (4.0 * chirpObj.TotalChirpPeriod);
    fprintf('  Max Unambig Range : %.1f m | Max Radial Velocity : %.1f m/s\n\n', chirpObj.MaxUnambiguousRange, vMax);
    
    % 2. Define Active Highway Scenario Targets
    % Target 1: Vehicle cruising ahead in right lane
    tgt1 = chirp.sim.TargetModel(...
        'Position', [1.5; 35.0; 0.0], ...
        'Velocity', [0.0; 12.0; 0.0], ...
        'Rcs', 10.0, ...
        'TargetId', 1 ...
    );
    
    % Target 2: Passing vehicle in left lane
    tgt2 = chirp.sim.TargetModel(...
        'Position', [-3.0; 20.0; 0.0], ...
        'Velocity', [0.0; -4.0; 0.0], ...
        'Rcs', 8.0, ...
        'TargetId', 2 ...
    );
    
    % Target 3: Crossing pedestrian near shoulder
    tgt3 = chirp.sim.TargetModel(...
        'Position', [-4.0; 12.0; 0.0], ...
        'Velocity', [1.2; 0.0; 0.0], ...
        'Rcs', 1.0, ...
        'TargetId', 3 ...
    );
    
    sim = chirp.sim.RadarSimulator(cfg, chirpObj, mimo);
    sim.EnableNoise = true;
    
    dt = 0.05; % 20 Hz frame rate
    tracker = chirp.tracking.GnnTracker();
    rp = chirp.dsp.RangeProcessor('WindowType', 'blackman-harris');
    dp = chirp.dsp.DopplerProcessor('WindowType', 'chebyshev');
    ca = chirp.cfar.CaCfar(1e-4, 8, 2, 4, 1);
    music = chirp.aoa.MusicEstimator(mimo, chirpObj.Wavelength, true);
    dbscan = chirp.pointcloud.VelocityAwareDbscan(2.0, 2.5, 2);
    
    outDir = fullfile(projectRoot, 'output');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    
    fprintf('Simulating 6 consecutive radar frames (0.30 s elapsed)...\n');
    
    for frame = 1:6
        time = (frame - 1) * dt;
        
        % Synthesize radar cube for current frame
        cube = sim.simulateFrame([tgt1, tgt2, tgt3]);
        
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
        pts3D = zeros(3, nDet);
        vels = zeros(1, nDet);
        snrs = zeros(1, nDet);
        
        for k = 1:nDet
            r = rAxis(rIdx(k));
            v = vAxis(vIdx(k));
            snapshot = squeeze(rdmCube(rIdx(k), vIdx(k), :));
            
            [specM, azGrid] = music.estimateAzimuth(snapshot, 1, -30:1:30, 0);
            [~, pkAz] = max(specM);
            azDeg = azGrid(pkAz);
            azRad = azDeg * (pi / 180.0);
            
            detectionsSpherical(:, k) = [r; azRad; 0.0; v];
            pts3D(:, k) = [r * sin(azRad); r * cos(azRad); 0.0];
            vels(k) = v;
            snrs(k) = pwrLin(rIdx(k), vIdx(k));
        end
        
        % 4D Point Cloud Clustering
        cloud = chirp.pointcloud.PointCloud4D(pts3D, vels, snrs, ones(1, nDet), time);
        [labels, nClusters] = dbscan.cluster(cloud);
        
        % Multi-target tracking step
        activeTracks = tracker.step(detectionsSpherical, time, dt);
        
        % Propagate targets for next frame
        tgt1 = tgt1.propagate(dt);
        tgt2 = tgt2.propagate(dt);
        tgt3 = tgt3.propagate(dt);
        
        if frame == 6
            fprintf('\n--- Frame 6 Final Perception Output ---\n');
            fprintf('  Raw CFAR Detections   : %d points\n', nDet);
            fprintf('  Identified Clusters   : %d dynamic clusters\n', nClusters);
            fprintf('  Active Tracked Objects: %d tracks\n\n', length(activeTracks));
            
            for trkIdx = 1:length(activeTracks)
                trk = activeTracks(trkIdx);
                fprintf('  Track #%d [%9s] | Pos: [%5.2f, %5.2f, %5.2f] m | Vel: [%5.2f, %5.2f, %5.2f] m/s | Hits: %d\n', ...
                    trk.TrackId, trk.Status, trk.State(1), trk.State(2), trk.State(3), trk.State(4), trk.State(5), trk.State(6), trk.Hits);
            end
            
            % Render Terminal Heatmap
            fprintf('\nRange-Doppler Power Spectrum Heatmap:\n');
            termHeat = chirp.viz.TerminalHeatmap(50, 18);
            termHeat.render(pwrDb, rAxis, vAxis);
            
            % Export vector SVG point cloud
            svgPath = fullfile(outDir, 'perception_pointcloud.svg');
            chirp.viz.SvgPointCloudRenderer.renderToFile(svgPath, cloud);
            fprintf('\nExported vector SVG visualization to: %s\n', svgPath);
        end
    end
    
    fprintf('\n4D Automotive Radar demo completed successfully.\n');
end
