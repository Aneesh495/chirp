function run_all_tests()
    % RUN_ALL_TESTS Master test harness running all unit, integration, and audit tests.
    
    startTime = tic();
    fprintf('=================================================================\n');
    fprintf('  CHIRP RADAR & SAR ENGINE - AUTOMATED VERIFICATION HARNESS      \n');
    fprintf('=================================================================\n\n');
    
    testDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(testDir);
    addpath(fullfile(projectRoot, 'src'));
    addpath(fullfile(projectRoot, 'bin', 'mex'));
    addpath(fullfile(testDir, 'unit'));
    addpath(fullfile(testDir, 'integration'));
    addpath(fullfile(testDir, 'audit'));
    
    testSuites = {
        'test_waveform',
        'test_radar_sim',
        'test_dsp',
        'test_aoa',
        'test_cfar',
        'test_pointcloud',
        'test_tracking',
        'test_microdoppler',
        'test_sar',
        'test_io',
        'test_mex_parity',
        'test_end_to_end_perception',
        'test_negative_controls',
        'test_mathematical_invariants',
        'test_no_em_dashes'
    };
    
    numPassed = 0;
    numFailed = 0;
    failedTests = {};
    
    for k = 1:length(testSuites)
        testName = testSuites{k};
        t0 = tic();
        try
            feval(testName);
            elapsed = toc(t0);
            fprintf('>>> %-32s [ PASS ] (%.3f s)\n\n', testName, elapsed);
            numPassed = numPassed + 1;
        catch err
            elapsed = toc(t0);
            fprintf('>>> %-32s [ FAIL ] (%.3f s)\n', testName, elapsed);
            fprintf('    Error: %s\n\n', err.message);
            numFailed = numFailed + 1;
            failedTests{end + 1} = testName;
        end
    end
    
    totalTime = toc(startTime);
    fprintf('=================================================================\n');
    fprintf('  VERIFICATION CAMPAIGN SUMMARY                                  \n');
    fprintf('=================================================================\n');
    fprintf('  Total Suites Executed : %d\n', length(testSuites));
    fprintf('  Passed                : %d\n', numPassed);
    fprintf('  Failed                : %d\n', numFailed);
    fprintf('  Total Duration        : %.3f seconds\n', totalTime);
    fprintf('=================================================================\n');
    
    if numFailed > 0
        error('run_all_tests:Failed', '%d test suites failed verification!', numFailed);
    else
        fprintf('All verification gates passed with zero defects.\n\n');
    end
end
