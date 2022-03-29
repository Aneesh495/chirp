function test_aoa()
    % TEST_AOA Unit tests for 3D/4D Angle of Arrival super-resolution estimators.
    
    fprintf('Running test_aoa...\n');
    
    lambda = chirp.utils.PhysicalConstants.freq2wavelength(77e9);
    d = lambda / 2.0;
    mimo = chirp.sim.MimoArray.automotive3Tx4Rx(lambda);
    
    % Two target arrivals at -12.0 deg and +18.0 deg
    theta1 = -12.0;
    theta2 = 18.0;
    sv1 = mimo.steeringVector(theta1, 0, lambda);
    sv2 = mimo.steeringVector(theta2, 0, lambda);
    
    nSnap = 60;
    s1 = exp(1j * 2 * pi * rand(1, nSnap));
    s2 = exp(1j * 2 * pi * rand(1, nSnap));
    noise = 0.05 * (randn(12, nSnap) + 1j * randn(12, nSnap));
    X = sv1 * s1 + sv2 * s2 + noise;
    
    % 1. Bartlett beamformer
    bb = chirp.aoa.BartlettBeamformer(mimo, lambda);
    [specB, azB] = bb.estimateAzimuth(X, -40:0.5:40, 0);
    assert(length(specB) == length(azB), 'Bartlett spectrum dimension mismatch');
    
    % 2. Capon MVDR beamformer
    cb = chirp.aoa.CaponBeamformer(mimo, lambda, 1e-3);
    [specC, azC] = cb.estimateAzimuth(X, -40:0.5:40, 0);
    assert(max(specC) > min(specC) + 10.0, 'Capon peak dynamic range too low');
    
    % 3. MUSIC pseudo-spectrum
    music = chirp.aoa.MusicEstimator(mimo, lambda, true);
    [specM, azM] = music.estimateAzimuth(X, 2, -30:0.2:30, 0);
    [~, pk1Idx] = max(specM);
    fprintf('  MUSIC primary peak: %.2f deg\n', azM(pk1Idx));
    
    % 4. Root-MUSIC (regression test for phase sign convention)
    rmusic = chirp.aoa.RootMusicEstimator(d, lambda);
    angRoots = rmusic.estimate(X, 2);
    fprintf('  Root-MUSIC angles: %.2f deg, %.2f deg (True: -12.0, 18.0)\n', angRoots(1), angRoots(2));
    assert(abs(angRoots(1) - theta1) < 1.0, 'Root-MUSIC theta1 error');
    assert(abs(angRoots(2) - theta2) < 1.0, 'Root-MUSIC theta2 error');
    
    % 5. ESPRIT TLS
    esprit = chirp.aoa.EspritEstimator(d, lambda);
    angEsprit = esprit.estimate(X, 2);
    fprintf('  ESPRIT angles: %.2f deg, %.2f deg (True: -12.0, 18.0)\n', angEsprit(1), angEsprit(2));
    assert(abs(angEsprit(1) - theta1) < 1.0, 'ESPRIT theta1 error');
    assert(abs(angEsprit(2) - theta2) < 1.0, 'ESPRIT theta2 error');
    
    % 6. Monopulse estimator
    mono = chirp.aoa.MonopulseEstimator(3.0 * pi, 15.0);
    % Simulate target at 2.5 deg off boresight
    svTarget = mimo.steeringVector(2.5, 0, lambda);
    chA = sum(svTarget(1:6));
    chB = sum(svTarget(7:12));
    offBoresight = mono.estimate(chA, chB);
    fprintf('  Monopulse estimated angle: %.2f deg (True: 2.5 deg)\n', offBoresight);
    assert(abs(offBoresight - 2.5) < 1.5, 'Monopulse estimation error too large');
    
    fprintf('  PASSED: test_aoa\n');
end
