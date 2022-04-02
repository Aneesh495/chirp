function results = bench_sar_backprojection()
    % BENCH_SAR_BACKPROJECTION Benchmark Synthetic Aperture Radar Backprojection reconstruction.
    %
    % Evaluates exact time-domain Backprojection focusing latency and pixel throughput
    % comparing vectorized MATLAB/Octave engine vs compiled C MEX kernel.
    
    fprintf('Running SAR Backprojection Benchmark...\n');
    
    thisDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(thisDir);
    addpath(fullfile(projectRoot, 'src'));
    addpath(fullfile(projectRoot, 'bin', 'mex'));
    
    gridConfigs = [ ...
        32,  32, 32; ...
        64,  64, 64; ...
        128, 128, 64 ...
    ];
    
    results = struct('Config', {}, 'Pixels', {}, 'TimeNativeMs', {}, 'TimeMexMs', {}, ...
                     'PixelsPerSecNative', {}, 'PixelsPerSecMex', {}, 'Speedup', {});
    
    hasMex = (exist('chirp_backprojection_mex', 'file') == 3);
    
    c = chirp.utils.PhysicalConstants.C;
    fc = 10.0e9;
    lambda = c / fc;
    fs = 100.0e6;
    
    for i = 1:size(gridConfigs, 1)
        nX = gridConfigs(i, 1);
        nY = gridConfigs(i, 2);
        nPulses = gridConfigs(i, 3);
        totalPixels = nX * nY;
        
        cfgStr = sprintf('%dx%d (P=%d)', nX, nY, nPulses);
        
        % Generate synthetic platform trajectory and range profiles
        trajX = linspace(-20.0, 20.0, nPulses);
        trajY = zeros(1, nPulses);
        trajZ = 200.0 * ones(1, nPulses);
        trajPos = [trajX; trajY; trajZ];
        
        nR = 256;
        rangeProfiles = randn(nR, nPulses) + 1j * randn(nR, nPulses);
        rangeAxis = linspace(180.0, 250.0, nR);
        
        gridX = linspace(-15.0, 15.0, nX);
        gridY = linspace(180.0, 220.0, nY);
        gridZ = 0.0;
        
        engine = chirp.sar.BackprojectionEngine(lambda, false);
        
        % 1. Native Backprojection benchmark
        tStart = tic();
        imgNative = engine.reconstruct2D(rangeProfiles, rangeAxis, trajPos, gridX, gridY, gridZ);
        elapsedNative = toc(tStart);
        timeNativeMs = elapsedNative * 1000.0;
        pxSecNative = totalPixels / elapsedNative;
        
        % 2. C MEX Backprojection benchmark
        if hasMex
            % Warmup
            dummy = chirp_backprojection_mex(rangeProfiles, rangeAxis, trajPos, gridX, gridY, gridZ, lambda);
            
            tStart = tic();
            imgMex = chirp_backprojection_mex(rangeProfiles, rangeAxis, trajPos, gridX, gridY, gridZ, lambda);
            elapsedMex = toc(tStart);
            timeMexMs = elapsedMex * 1000.0;
            pxSecMex = totalPixels / elapsedMex;
            speedup = elapsedNative / elapsedMex;
        else
            timeMexMs = NaN;
            pxSecMex = NaN;
            speedup = 1.0;
        end
        
        results(i).Config = cfgStr;
        results(i).Pixels = totalPixels;
        results(i).TimeNativeMs = timeNativeMs;
        results(i).TimeMexMs = timeMexMs;
        results(i).PixelsPerSecNative = pxSecNative;
        results(i).PixelsPerSecMex = pxSecMex;
        results(i).Speedup = speedup;
        
        fprintf('  Config: %16s | Native: %6.1f ms (%7.0f px/s) | MEX: %6.1f ms (%7.0f px/s) | Speedup: %4.2fx\n', ...
            cfgStr, timeNativeMs, pxSecNative, timeMexMs, pxSecMex, speedup);
    end
end
