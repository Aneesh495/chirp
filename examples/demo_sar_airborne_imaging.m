function demo_sar_airborne_imaging()
    % DEMO_SAR_AIRBORNE_IMAGING Synthetic Aperture Radar (SAR) 2D imaging demonstration.
    %
    % Simulates an airborne SAR platform illuminating ground scatterers along a
    % synthetic baseline, evaluating exact time-domain Backprojection vs frequency-domain
    % Range-Doppler reconstruction algorithms.
    
    fprintf('=================================================================\n');
    fprintf('  CHIRP DEMO: AIRBORNE SYNTHETIC APERTURE RADAR (SAR) IMAGING    \n');
    fprintf('=================================================================\n\n');
    
    thisDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(thisDir);
    addpath(fullfile(projectRoot, 'src'));
    addpath(fullfile(projectRoot, 'bin', 'mex'));
    
    % 1. Platform Kinematics and Geometry
    c = chirp.utils.PhysicalConstants.C;
    fc = 10.0e9; % X-band 10 GHz
    lambda = c / fc;
    
    altitude = 400.0;     % Altitude (m)
    velocity = 75.0;      % Flight speed (m/s)
    apertureLen = 40.0;   % Synthetic aperture baseline (m)
    nPulses = 64;         % Slow-time pulse count
    
    geo = chirp.sar.SarGeometry('FlightAltitude', altitude, ...
                                'PlatformVelocity', velocity, ...
                                'ApertureLength', apertureLen);
                                
    fprintf('SAR Mission Geometry:\n');
    fprintf('  Carrier Frequency : %.2f GHz (Wavelength = %.2f cm)\n', fc / 1e9, lambda * 100.0);
    fprintf('  Platform Altitude : %.1f m | Velocity: %.1f m/s\n', altitude, velocity);
    fprintf('  Synthetic Baseline: %.1f m across %d azimuth pulses\n\n', apertureLen, nPulses);
    
    % 2. Synthesize Synthetic Aperture Trajectory
    uCoords = linspace(-apertureLen / 2.0, apertureLen / 2.0, nPulses);
    platPos = zeros(3, nPulses);
    for k = 1:nPulses
        platPos(:, k) = geo.platformPosition(uCoords(k));
    end
    
    % 3. Ground Target Scatterers (X = cross-range, Y = ground-range)
    groundTargets = [ ...
        -10.0, 500.0, 0.0, 1.0; ...  % Target 1: Left scatterer
          0.0, 500.0, 0.0, 2.0; ...  % Target 2: Center scatterer
         12.0, 510.0, 0.0, 1.5       % Target 3: Right/deeper scatterer
    ];
    
    nRange = 128;
    rAxis = linspace(450.0, 750.0, nRange);
    dr = rAxis(2) - rAxis(1);
    
    % Simulate range profiles across all pulses
    rProfiles = zeros(nRange, nPulses);
    for p = 1:nPulses
        pPos = platPos(:, p);
        for t = 1:size(groundTargets, 1)
            tgtPos = groundTargets(t, 1:3).';
            rcs = groundTargets(t, 4);
            dist = norm(pPos - tgtPos);
            binIdx = round((dist - rAxis(1)) / dr) + 1;
            if binIdx >= 1 && binIdx <= nRange
                phase = exp(-1j * (4.0 * pi / lambda) * dist);
                rProfiles(binIdx, p) = rProfiles(binIdx, p) + rcs * phase;
            end
        end
    end
    
    % Add low thermal noise floor
    rProfiles = rProfiles + 0.02 * (randn(size(rProfiles)) + 1j * randn(size(rProfiles)));
    
    % 4. Time-Domain Backprojection Image Reconstruction
    fprintf('Executing Exact Time-Domain Backprojection Image Formation...\n');
    gX = linspace(-25.0, 25.0, 64);
    gY = linspace(480.0, 530.0, 64);
    
    bp = chirp.sar.BackprojectionEngine(lambda, false);
    tBp = tic();
    sarImage = bp.reconstruct2D(rProfiles, rAxis, platPos, gX, gY, 0.0);
    elapsedBp = toc(tBp);
    
    pwrBp = abs(sarImage).^2;
    pwrDb = 10.0 * log10(max(pwrBp / max(pwrBp(:)), 1e-6));
    
    fprintf('  Backprojection completed in %6.2f ms (64x64 grid = 4096 pixels)\n', elapsedBp * 1000.0);
    
    % 5. C MEX Acceleration Comparison
    if exist('chirp_backprojection_mex', 'file') == 3
        tMex = tic();
        imgMex = chirp_backprojection_mex(rProfiles, rAxis, platPos, gX, gY, 0.0, lambda);
        elapsedMex = toc(tMex);
        speedup = elapsedBp / elapsedMex;
        fprintf('  C MEX Backprojection completed in %6.2f ms (Speedup: %4.1fx)\n', elapsedMex * 1000.0, speedup);
    end
    
    % 6. Render Terminal ASCII Map of Reconstructed Scene
    fprintf('\nReconstructed 2D SAR Image Map (Cross-Range X vs Ground-Range Y):\n');
    term = chirp.viz.TerminalHeatmap(50, 16);
    term.render(pwrDb, gX, gY);
    
    fprintf('\nAirborne SAR Imaging demo completed successfully.\n');
end
