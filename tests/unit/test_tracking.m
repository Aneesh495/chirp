function test_tracking()
    % TEST_TRACKING Unit tests for kinematic filtering and multi-target tracking.
    
    fprintf('Running test_tracking...\n');
    
    % 1. Motion model state transition
    dt = 0.1;
    [F_cv, Q_cv] = chirp.tracking.MotionModel.constantVelocity(dt, 0.5);
    assert(size(F_cv, 1) == 6 && size(F_cv, 2) == 6, 'CV state dimension mismatch');
    assert(all(eig(Q_cv) >= -1e-12), 'Process covariance Q must be positive semi-definite');
    
    % 2. Extended Kalman Filter (EKF) step
    ekf = chirp.tracking.ExtendedKalmanFilter();
    initState = [0.0; 20.0; 0.0; 1.0; 8.0; 0.0];
    initCov = eye(6);
    
    [pState, pCov] = ekf.predict(initState, initCov, dt);
    assert(abs(pState(1) - 0.1) < 1e-6, 'Predicted X position mismatch');
    assert(abs(pState(2) - 20.8) < 1e-6, 'Predicted Y position mismatch');
    
    % Analytical observation: target at [0.1; 20.8; 0] moving at [1; 8; 0]
    rTrue = norm(pState(1:3));
    azTrue = atan2(pState(1), pState(2));
    vrTrue = sum(pState(1:3) .* pState(4:6)) / rTrue;
    z = [rTrue; azTrue; 0.0; vrTrue];
    
    [uState, uCov, innov, S] = ekf.update(pState, pCov, z);
    assert(abs(innov(1)) < 1e-6, 'Innovation for exact measurement must be near zero');
    assert(trace(uCov) < trace(pCov), 'Measurement update must reduce covariance uncertainty');
    
    % 3. Hungarian Matcher assignment
    costMatrix = [1.0, 5.0, 10.0; ...
                  4.0, 2.0, 8.0; ...
                  9.0, 7.0, 3.0];
    [matches, unR, unC] = chirp.tracking.HungarianMatcher.match(costMatrix, 15.0);
    assert(size(matches, 1) == 3, 'All 3 rows should be assigned');
    % Optimal pairing is (1,1), (2,2), (3,3)
    assert(matches(1, 1) == 1 && matches(1, 2) == 1, 'Match (1, 1) expected');
    assert(matches(2, 1) == 2 && matches(2, 2) == 2, 'Match (2, 2) expected');
    assert(matches(3, 1) == 3 && matches(3, 2) == 3, 'Match (3, 3) expected');
    
    % 4. TrackLifecycleManager duplicate merging
    t1 = chirp.tracking.RadarTrack(1, [1.0; 20.0; 0; 0; 5; 0], 0.1 * eye(6), 0);
    t2 = chirp.tracking.RadarTrack(2, [1.05; 20.05; 0; 0; 5; 0], 0.1 * eye(6), 0); % Duplicate
    t3 = chirp.tracking.RadarTrack(3, [10.0; 40.0; 0; 0; 5; 0], 0.1 * eye(6), 0); % Distinct
    
    tlm = chirp.tracking.TrackLifecycleManager();
    cleaned = tlm.pruneAndMerge([t1, t2, t3]);
    assert(length(cleaned) == 2, 'Duplicate track was not merged');
    
    fprintf('  PASSED: test_tracking\n');
end
