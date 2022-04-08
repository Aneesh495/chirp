function demo_mimo_virtual_array()
    % DEMO_MIMO_VIRTUAL_ARRAY MIMO virtual array synthesis and super-resolution AoA.
    %
    % Simulates physical 3-Tx and 4-Rx antennas, constructs the 12-element virtual
    % uniform linear array, and evaluates 5 direction-of-arrival estimators:
    %   1. Classical Bartlett Fourier Beamformer
    %   2. Capon Minimum Variance Distortionless Response (MVDR)
    %   3. Multiple Signal Classification (MUSIC)
    %   4. Polynomial Root-MUSIC
    %   5. Estimation of Signal Parameters via Rotational Invariance (ESPRIT)
    % Demonstrates resolution of two coherent targets spaced closer than Rayleigh limit.
    
    fprintf('=================================================================\n');
    fprintf('  CHIRP DEMO: MIMO VIRTUAL ARRAY & SUPER-RESOLUTION AOA          \n');
    fprintf('=================================================================\n\n');
    
    thisDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(thisDir);
    addpath(fullfile(projectRoot, 'src'));
    addpath(fullfile(projectRoot, 'bin', 'mex'));
    
    % 1. Synthesize Automotive MIMO Array
    fc = 77.0e9;
    lambda = chirp.utils.PhysicalConstants.freq2wavelength(fc);
    d = lambda / 2.0;
    
    mimo = chirp.sim.MimoArray.automotive3Tx4Rx(lambda);
    
    fprintf('MIMO Array Geometry:\n');
    fprintf('  Physical TX Antennas: %d elements (Spacing = 4 * d = 2.0 * lambda)\n', mimo.NumTx);
    fprintf('  Physical RX Antennas: %d elements (Spacing = 1 * d = 0.5 * lambda)\n', mimo.NumRx);
    fprintf('  Virtual Array Total : %d co-array elements\n', mimo.NumVirtual);
    
    virtX_d = mimo.VirtualPositions(1, :) / d;
    fprintf('  Virtual Positions (in lambda/2 units): [ ');
    fprintf('%.1f ', virtX_d);
    fprintf(']\n');
    
    % Rayleigh Angular Resolution Limit: theta_res = lambda / D = 1 / (N - 1)
    rayleighLimitDeg = (lambda / (mimo.NumVirtual * d)) * (180.0 / pi);
    fprintf('  Theoretical Rayleigh 3-dB Diffraction Limit: %.2f degrees\n\n', rayleighLimitDeg);
    
    % 2. Synthesize Two Targets Below Rayleigh Diffraction Limit
    % Rayleigh limit is ~9.5 deg. We space targets by only 5.0 deg (-2.5 deg and +2.5 deg)!
    theta1 = -2.5;
    theta2 =  2.5;
    fprintf('Synthesizing two targets at theta1 = %+.1f deg and theta2 = %+.1f deg (Delta = %.1f deg < Rayleigh)\n', ...
        theta1, theta2, theta2 - theta1);
        
    sv1 = mimo.steeringVector(theta1, 0.0, lambda);
    sv2 = mimo.steeringVector(theta2, 0.0, lambda);
    
    nSnapshots = 100;
    s1 = exp(1j * 2 * pi * rand(1, nSnapshots));
    s2 = exp(1j * 2 * pi * rand(1, nSnapshots));
    noise = 0.05 * (randn(mimo.NumVirtual, nSnapshots) + 1j * randn(mimo.NumVirtual, nSnapshots));
    
    X = sv1 * s1 + sv2 * s2 + noise;
    R = (X * X') / nSnapshots;
    
    % 3. Evaluate Direction-of-Arrival Estimators
    azSearchGrid = -30.0:0.25:30.0;
    
    % A. Bartlett Fourier Beamformer
    bb = chirp.aoa.BartlettBeamformer(mimo, lambda);
    specBartlett = bb.estimateAzimuth(X, azSearchGrid, 0.0);
    
    % B. Capon MVDR Beamformer
    cb = chirp.aoa.CaponBeamformer(mimo, lambda, 1e-4);
    specCapon = cb.estimateAzimuth(X, azSearchGrid, 0.0);
    
    % C. MUSIC Super-Resolution
    music = chirp.aoa.MusicEstimator(mimo, lambda, true);
    specMusic = music.estimateAzimuth(X, 2, azSearchGrid, 0.0);
    
    % D. Root-MUSIC
    rmusic = chirp.aoa.RootMusicEstimator(d, lambda);
    angsRoot = rmusic.estimate(X, 2);
    
    % E. ESPRIT TLS
    esprit = chirp.aoa.EspritEstimator(d, lambda);
    angsEsprit = esprit.estimate(X, 2);
    
    fprintf('\n--- Estimation Results ---\n');
    fprintf('  True Target Angles  : %+.2f deg, %+.2f deg\n', theta1, theta2);
    fprintf('  Root-MUSIC Estimates: %+.2f deg, %+.2f deg (Error: %.2f deg)\n', ...
        angsRoot(1), angsRoot(2), max(abs(angsRoot - [theta1, theta2])));
    fprintf('  ESPRIT Estimates    : %+.2f deg, %+.2f deg (Error: %.2f deg)\n', ...
        angsEsprit(1), angsEsprit(2), max(abs(angsEsprit - [theta1, theta2])));
        
    % 4. Render Terminal PPI Polar Plot Display
    rRange = linspace(5.0, 50.0, 32);
    rangeAzMatrix = zeros(length(rRange), length(azSearchGrid));
    % Place the two targets at range = 25.0 m
    tgtRangeBin = 16;
    rangeAzMatrix(tgtRangeBin, :) = specMusic;
    
    fprintf('\nRange-Azimuth PPI Radar Scope Display:\n');
    chirp.viz.PolarPlotRenderer.renderTerminal(rangeAzMatrix, rRange, azSearchGrid, 55, 18);
    
    fprintf('\nMIMO Virtual Array and AoA demo completed successfully.\n');
end
