function chirp_cli(varargin)
    % CHIRP Unified Command-Line Interface for FMCW Radar & SAR Engine.
    %
    % Usage:
    %   chirp version
    %   chirp sim [options]
    %   chirp process [options]
    %   chirp cfar [options]
    %   chirp track [options]
    %   chirp sar [options]
    %   chirp microdoppler [options]
    %   chirp bench
    %   chirp verify
    
    % Ensure src directory is in path
    thisDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(thisDir);
    addpath(fullfile(projectRoot, 'src'));
    addpath(fullfile(projectRoot, 'bin', 'mex'));
    
    if nargin == 0 || strcmp(varargin{1}, '--help') || strcmp(varargin{1}, '-h')
        printHelp();
        return;
    end
    
    subcmd = lower(varargin{1});
    args = varargin(2:end);
    
    switch subcmd
        case 'version'
            v = chirp.version();
            fprintf('chirp version %d.%d.%d (%s)\n', v.major, v.minor, v.patch, v.release);
            fprintf('%s\n', v.description);
            hasMex = (exist('chirp_cfar_2d_mex', 'file') == 3);
            if hasMex
                fprintf('Acceleration: Native C MEX acceleration kernels compiled and active.\n');
            else
                fprintf('Acceleration: Pure vectorized MATLAB / Octave execution mode.\n');
            end
            
        case 'sim'
            fprintf('[chirp sim] Synthesizing 4D FMCW radar cube with 2 active targets...\n');
            cfg = chirp.sim.RadarConfig('NumChirps', 64);
            chirpObj = chirp.waveform.WaveformGenerator.automotiveLRR();
            chirpObj.NumAdcSamples = 256;
            mimo = chirp.sim.MimoArray.automotive3Tx4Rx(chirpObj.Wavelength);
            
            tgt1 = chirp.sim.TargetModel('Position', [2.0; 25.0; 0.0], 'Velocity', [0.0; 12.0; 0.0], 'Rcs', 10.0);
            tgt2 = chirp.sim.TargetModel('Position', [-4.0; 45.0; 0.0], 'Velocity', [0.0; -8.0; 0.0], 'Rcs', 5.0);
            
            sim = chirp.sim.RadarSimulator(cfg, chirpObj, mimo);
            sim.EnableNoise = true;
            [cube, tGrid] = sim.simulateFrame([tgt1, tgt2]);
            
            outPath = 'output/simulated_radar_cube.rcube';
            [outDir, ~, ~] = fileparts(outPath);
            if ~exist(outDir, 'dir'), mkdir(outDir); end
            chirp.io.RadarCubeCodec.save(outPath, cube);
            fprintf('Simulation complete. Radar cube [%d x %d x %d] written to %s\n', ...
                size(cube, 1), size(cube, 2), size(cube, 3), outPath);
                
        case 'process'
            fprintf('[chirp process] Executing 1D Range-FFT and 2D Doppler-FFT...\n');
            cubePath = 'output/simulated_radar_cube.rcube';
            if ~exist(cubePath, 'file')
                fprintf('No simulated cube found. Running sim first...\n');
                chirp_cli('sim');
            end
            [cube, dims] = chirp.io.RadarCubeCodec.load(cubePath);
            chirpObj = chirp.waveform.WaveformGenerator.automotiveLRR();
            chirpObj.NumAdcSamples = dims(1);
            
            rp = chirp.dsp.RangeProcessor();
            [rangeProf, rAxis] = rp.process(cube, chirpObj);
            
            dp = chirp.dsp.DopplerProcessor();
            [rdmCube, vAxis] = dp.process(rangeProf, chirpObj);
            
            rdm = chirp.dsp.RangeDopplerMap(rdmCube, rAxis, vAxis);
            pwrDb = rdm.nonCoherentIntegration();
            
            chirp.viz.TerminalHeatmap.render(pwrDb, rAxis, vAxis, 'Range-Doppler Power Spectrum', 60, 18);
            
        case 'cfar'
            fprintf('[chirp cfar] Running 2D CA-CFAR detection and point cloud extraction...\n');
            cubePath = 'output/simulated_radar_cube.rcube';
            if ~exist(cubePath, 'file'), chirp_cli('sim'); end
            [cube, dims] = chirp.io.RadarCubeCodec.load(cubePath);
            chirpObj = chirp.waveform.WaveformGenerator.automotiveLRR();
            chirpObj.NumAdcSamples = dims(1);
            
            rp = chirp.dsp.RangeProcessor();
            [rangeProf, rAxis] = rp.process(cube, chirpObj);
            dp = chirp.dsp.DopplerProcessor();
            [rdmCube, vAxis] = dp.process(rangeProf, chirpObj);
            rdm = chirp.dsp.RangeDopplerMap(rdmCube, rAxis, vAxis);
            pwrDb = rdm.nonCoherentIntegration();
            pwrLin = 10.^(pwrDb / 10.0);
            
            ca = chirp.cfar.CaCfar(1e-4, 8, 2, 4, 1);
            [mask, ~, ~] = ca.detect(pwrLin);
            [rIdx, vIdx] = find(mask);
            
            fprintf('Detected %d CFAR candidate points.\n', length(rIdx));
            
            % Angle of arrival estimation per detected point
            lambda = chirpObj.Wavelength;
            mimo = chirp.sim.MimoArray.automotive3Tx4Rx(lambda);
            music = chirp.aoa.MusicEstimator(mimo, lambda, true);
            
            ptsPos = zeros(3, length(rIdx));
            ptsVel = zeros(1, length(rIdx));
            ptsSnr = zeros(1, length(rIdx));
            ptsRcs = ones(1, length(rIdx));
            
            for k = 1:length(rIdx)
                r = rAxis(rIdx(k));
                v = vAxis(vIdx(k));
                spatialSnapshot = squeeze(rdmCube(rIdx(k), vIdx(k), :));
                
                [specM, azGrid] = music.estimateAzimuth(spatialSnapshot, 1, -45:1:45, 0);
                [~, pkAz] = max(specM);
                azDeg = azGrid(pkAz);
                
                [x, y, z] = chirp.utils.PhysicalConstants.sph2cartDeg(r, azDeg, 0);
                ptsPos(:, k) = [x; y; z];
                ptsVel(k) = v;
                ptsSnr(k) = pwrDb(rIdx(k), vIdx(k));
            end
            
            cloud = chirp.pointcloud.PointCloud4D(ptsPos, ptsVel, ptsSnr, ptsRcs, 0.0);
            chirp.viz.TerminalPointCloud.render(cloud, [-15, 15], [0, 60], 51, 15);
            
        case 'track'
            fprintf('[chirp track] Running multi-target kinematic tracker demonstration...\n');
            dt = 0.05;
            tracker = chirp.tracking.GnnTracker();
            
            % 10-frame trajectory
            p1 = [2.0; 20.0; 0.0]; v1 = [0.0; 12.0; 0.0];
            p2 = [-5.0; 35.0; 0.0]; v2 = [0.0; -6.0; 0.0];
            
            for f = 1:10
                p1 = p1 + v1 * dt;
                p2 = p2 + v2 * dt;
                
                r1 = norm(p1); az1 = atan2(p1(1), p1(2)); vr1 = sum(p1 .* v1)/r1;
                r2 = norm(p2); az2 = atan2(p2(1), p2(2)); vr2 = sum(p2 .* v2)/r2;
                
                z1 = [r1 + 0.05*randn(); az1; 0; vr1 + 0.05*randn()];
                z2 = [r2 + 0.05*randn(); az2; 0; vr2 + 0.05*randn()];
                
                activeTracks = tracker.step([z1, z2], (f - 1) * dt, dt);
            end
            
            fprintf('Active Tracker Results (%d confirmed tracks):\n', length(activeTracks));
            for t = 1:length(activeTracks)
                trk = activeTracks(t);
                fprintf('  Track #%d | Status: %-10s | Pos: [%6.2f, %6.2f, %6.2f] m | Vel: [%6.2f, %6.2f, %6.2f] m/s\n', ...
                    trk.TrackId, trk.Status, trk.State(1), trk.State(2), trk.State(3), ...
                    trk.State(4), trk.State(5), trk.State(6));
            end
            
        case 'sar'
            fprintf('[chirp sar] Running Synthetic Aperture Radar Backprojection...\n');
            lambda = 0.03;
            geo = chirp.sar.SarGeometry('FlightAltitude', 300.0, 'PlatformVelocity', 50.0, 'ApertureLength', 40.0);
            
            nPulses = 32;
            uCoords = linspace(-20.0, 20.0, nPulses);
            platPos = zeros(3, nPulses);
            for k = 1:nPulses, platPos(:, k) = geo.platformPosition(uCoords(k)); end
            
            tgtGround = [0.0; 500.0; 0.0];
            nRange = 64;
            rAxis = linspace(400.0, 700.0, nRange);
            dr = rAxis(2) - rAxis(1);
            rProfiles = zeros(nRange, nPulses);
            
            for k = 1:nPulses
                R = norm(platPos(:, k) - tgtGround);
                b = round((R - rAxis(1)) / dr) + 1;
                if b >= 1 && b <= nRange
                    rProfiles(b, k) = exp(-1j * (4.0 * pi / lambda) * R);
                end
            end
            
            bp = chirp.sar.BackprojectionEngine(lambda, true);
            gX = linspace(-15.0, 15.0, 31);
            gY = linspace(485.0, 515.0, 31);
            sarImg = bp.reconstruct2D(rProfiles, rAxis, platPos, gX, gY, 0.0);
            
            imgDb = 20.0 * log10(max(abs(sarImg), eps));
            chirp.viz.TerminalHeatmap.render(imgDb, gY, gX, 'SAR Focused Image (Backprojection)', 50, 16);
            
        case 'microdoppler'
            fprintf('[chirp microdoppler] Analyzing pedestrian micro-Doppler signature...\n');
            lambda = chirp.utils.PhysicalConstants.freq2wavelength(77e9);
            fs = 2000.0;
            t = 0:(1/fs):1.0;
            f_torso = 2.0 * 1.4 / lambda;
            f_limb = 2.0 * 1.2 / lambda;
            sig = exp(1j * 2*pi * f_torso * t) + 0.4 * exp(1j * 2*pi * (f_torso * t + (f_limb / (2*pi*1.8)) * sin(2*pi*1.8*t)));
            
            stft = chirp.microdoppler.StftAnalyzer('WindowLength', 64, 'HopSize', 16, 'Nfft', 128);
            [pwrDb, tAx, fAx] = stft.compute(sig, fs);
            
            feats = chirp.microdoppler.MicroDopplerFeatures.extract(10.^(pwrDb/10), tAx, fAx);
            [cName, conf, scores] = chirp.microdoppler.KinematicClassifier.classify(feats, lambda);
            
            chirp.viz.TerminalHeatmap.render(pwrDb, fAx, tAx, 'Micro-Doppler Spectrogram', 50, 14);
            fprintf('Classification: %s (Confidence: %.1f%%)\n', upper(cName), conf * 100);
            
        case 'bench'
            addpath(fullfile(projectRoot, 'benchmarks'));
            run_benchmarks;
            
        case 'verify'
            addpath(fullfile(projectRoot, 'tests'));
            run_all_tests;
            
        otherwise
            fprintf('Unknown subcommand: %s\n', subcmd);
            printHelp();
    end
end

function printHelp()
    fprintf('chirp - High-Resolution FMCW 4D Imaging Radar & SAR Perception Engine\n\n');
    fprintf('Usage: chirp <command> [options]\n\n');
    fprintf('Commands:\n');
    fprintf('  version       Display engine version and acceleration compilation status\n');
    fprintf('  sim           Synthesize multi-target 4D radar cube with thermal noise\n');
    fprintf('  process       Execute 1D Range-FFT and 2D Doppler-FFT spectrum\n');
    fprintf('  cfar          Run 2D CA-CFAR detection and extract 4D point cloud\n');
    fprintf('  track         Execute EKF/GNN kinematic multi-target tracking\n');
    fprintf('  sar           Run time-domain SAR backprojection image reconstruction\n');
    fprintf('  microdoppler  Generate time-frequency spectrogram and classify gait\n');
    fprintf('  bench         Run computational throughput and latency benchmarks\n');
    fprintf('  verify        Run full verification suite and negative controls\n');
end
