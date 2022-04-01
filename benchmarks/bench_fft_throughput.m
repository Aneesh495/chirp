function results = bench_fft_throughput()
    % BENCH_FFT_THROUGHPUT Benchmark Range and Doppler FFT throughput across radar cubes.
    %
    % Evaluates pure matrix FFT vs compiled C MEX range FFT kernel across small,
    % medium, and large automotive 4D radar cube configurations.
    
    fprintf('Running FFT Throughput Benchmark...\n');
    
    thisDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(thisDir);
    addpath(fullfile(projectRoot, 'src'));
    addpath(fullfile(projectRoot, 'bin', 'mex'));
    
    cubeConfigs = [ ...
        256,  64, 12; ...
        512, 128, 12; ...
        1024, 256, 16  ...
    ];
    
    results = struct('Config', {}, 'Samples', {}, 'TimeNativeMs', {}, 'TimeMexMs', {}, ...
                     'ThroughputNativeMSps', {}, 'ThroughputMexMSps', {}, 'Speedup', {});
    
    hasMex = (exist('chirp_range_fft_mex', 'file') == 3);
    
    for i = 1:size(cubeConfigs, 1)
        nR = cubeConfigs(i, 1);
        nV = cubeConfigs(i, 2);
        nC = cubeConfigs(i, 3);
        totalSamples = nR * nV * nC;
        
        cfgStr = sprintf('%dx%dx%d', nR, nV, nC);
        cube = randn(nR, nV, nC) + 1j * randn(nR, nV, nC);
        win = 0.5 - 0.5 * cos(2 * pi * (0:nR-1).' / (nR - 1));
        
        % 1. Native FFT Range processing benchmark
        nWarmup = 2;
        for w = 1:nWarmup
            dummy = fft(cube .* win, nR, 1);
        end
        
        nIters = 10;
        tStart = tic();
        for it = 1:nIters
            resNative = fft(cube .* win, nR, 1);
        end
        elapsedNative = toc(tStart) / nIters;
        
        timeNativeMs = elapsedNative * 1000.0;
        throughputNative = (totalSamples / elapsedNative) / 1e6;
        
        % 2. C MEX Range processing benchmark
        if hasMex
            for w = 1:nWarmup
                dummy = chirp_range_fft_mex(cube, win, int32(nR));
            end
            
            tStart = tic();
            for it = 1:nIters
                resMex = chirp_range_fft_mex(cube, win, int32(nR));
            end
            elapsedMex = toc(tStart) / nIters;
            timeMexMs = elapsedMex * 1000.0;
            throughputMex = (totalSamples / elapsedMex) / 1e6;
            speedup = elapsedNative / elapsedMex;
        else
            timeMexMs = NaN;
            throughputMex = NaN;
            speedup = 1.0;
        end
        
        results(i).Config = cfgStr;
        results(i).Samples = totalSamples;
        results(i).TimeNativeMs = timeNativeMs;
        results(i).TimeMexMs = timeMexMs;
        results(i).ThroughputNativeMSps = throughputNative;
        results(i).ThroughputMexMSps = throughputMex;
        results(i).Speedup = speedup;
        
        fprintf('  Config: %12s | Native: %6.2f ms (%6.1f MS/s) | MEX: %6.2f ms (%6.1f MS/s) | Speedup: %4.2fx\n', ...
            cfgStr, timeNativeMs, throughputNative, timeMexMs, throughputMex, speedup);
    end
end
