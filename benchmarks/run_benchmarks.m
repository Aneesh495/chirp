function run_benchmarks()
    % RUN_BENCHMARKS Master benchmark execution harness for chirp radar engine.
    %
    % Executes throughput and latency benchmarks across DSP transforms, CFAR detectors,
    % and SAR backprojection image focusers, reporting summary metrics.
    
    fprintf('=================================================================\n');
    fprintf('  CHIRP RADAR & SAR ENGINE - PERFORMANCE BENCHMARK HARNESS       \n');
    fprintf('=================================================================\n\n');
    
    thisDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(thisDir);
    addpath(fullfile(projectRoot, 'src'));
    addpath(fullfile(projectRoot, 'bin', 'mex'));
    addpath(thisDir);
    
    tTotal = tic();
    
    % 1. FFT Throughput
    fprintf('\n--- 1. Range & Doppler FFT Processing Throughput ---\n');
    bench_fft_throughput();
    
    % 2. CFAR Detectors Latency
    fprintf('\n--- 2. 2D CFAR Detection Algorithms Latency ---\n');
    bench_cfar_performance();
    
    % 3. SAR Backprojection
    fprintf('\n--- 3. SAR Time-Domain Backprojection Imaging Throughput ---\n');
    bench_sar_backprojection();
    
    elapsedTotal = toc(tTotal);
    fprintf('\n=================================================================\n');
    fprintf('  Benchmark harness completed successfully in %.2f seconds.\n', elapsedTotal);
    fprintf('=================================================================\n\n');
end
