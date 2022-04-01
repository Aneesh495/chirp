function results = bench_cfar_performance()
    % BENCH_CFAR_PERFORMANCE Benchmark 2D CFAR detector latency across algorithms and kernels.
    %
    % Evaluates Cell-Averaging, Greatest-Of, Smallest-Of, Ordered-Statistic, Trimmed-Mean,
    % and compiled C MEX sliding-window implementations on 2D Range-Doppler maps.
    
    fprintf('Running CFAR Performance Benchmark...\n');
    
    thisDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(thisDir);
    addpath(fullfile(projectRoot, 'src'));
    addpath(fullfile(projectRoot, 'bin', 'mex'));
    
    nR = 128;
    nV = 32;
    rdmPwr = abs(randn(nR, nV) + 1j * randn(nR, nV)).^2;
    rdmPwr(30, 10) = 100.0;
    rdmPwr(70, 20) = 150.0;
    
    gR = 2; gV = 1;
    tR = 4; tV = 2;
    pfa = 1e-4;
    
    algorithms = {'CA-CFAR', 'GO-CFAR', 'SO-CFAR', 'OS-CFAR'};
    hasMex = (exist('chirp_cfar_2d_mex', 'file') == 3);
    if hasMex
        algorithms{end+1} = 'CA-CFAR (C MEX)';
    end
    
    results = struct('Algorithm', {}, 'LatencyMs', {}, 'Fps', {}, 'Detections', {});
    nIters = 2;
    
    for i = 1:length(algorithms)
        algo = algorithms{i};
        
        switch algo
            case 'CA-CFAR'
                det = chirp.cfar.CaCfar(pfa, tR, gR, tV, gV);
                fn = @() det.detect(rdmPwr);
            case 'GO-CFAR'
                det = chirp.cfar.GoCfar(pfa, tR, gR, tV, gV);
                fn = @() det.detect(rdmPwr);
            case 'SO-CFAR'
                det = chirp.cfar.SoCfar(pfa, tR, gR, tV, gV);
                fn = @() det.detect(rdmPwr);
            case 'OS-CFAR'
                det = chirp.cfar.OsCfar(pfa, tR, gR, tV, gV, 24);
                fn = @() det.detect(rdmPwr);
            case 'CA-CFAR (C MEX)'
                fn = @() chirp_cfar_2d_mex(rdmPwr, pfa, int32(tR), int32(gR), int32(tV), int32(gV));
        end
        
        % Warmup
        [mask, ~] = fn();
        
        tStart = tic();
        for it = 1:nIters
            [mask, ~] = fn();
        end
        elapsed = toc(tStart) / nIters;
        latencyMs = elapsed * 1000.0;
        fps = 1.0 / elapsed;
        numDets = sum(mask(:));
        
        results(i).Algorithm = algo;
        results(i).LatencyMs = latencyMs;
        results(i).Fps = fps;
        results(i).Detections = numDets;
        
        fprintf('  Algo: %18s | Latency: %6.2f ms | FPS: %6.1f | Detections: %d\n', ...
            algo, latencyMs, fps, numDets);
    end
end
