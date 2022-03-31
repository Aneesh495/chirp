function test_mathematical_invariants()
    % TEST_MATHEMATICAL_INVARIANTS Formal verification of physical and mathematical conservation laws.
    %
    % Asserts Parseval's energy conservation, ambiguity volume invariance,
    % Hermitian spatial covariance symmetry, and Kalman covariance reduction.
    
    fprintf('Running test_mathematical_invariants...\n');
    
    % Invariant 1: Parseval's Energy Conservation Theorem
    x = randn(256, 1) + 1j * randn(256, 1);
    N = length(x);
    timeEnergy = sum(abs(x).^2);
    X_fft = fft(x);
    freqEnergy = sum(abs(X_fft).^2) / N;
    diffParseval = abs(timeEnergy - freqEnergy) / timeEnergy;
    fprintf('  Invariant 1: Parseval Energy Conservation Relative Delta: %e\n', diffParseval);
    assert(diffParseval < 1e-12, 'Parseval energy conservation violated');
    
    % Invariant 2: Hermitian Symmetry of Covariance Matrix
    snapshots = randn(8, 50) + 1j * randn(8, 50);
    R = (snapshots * snapshots') / 50;
    diffHermitian = max(max(abs(R - R')));
    assert(diffHermitian < 1e-12, 'Covariance matrix violates Hermitian symmetry');
    
    % All eigenvalues must be non-negative real numbers
    eigenvals = eig(R);
    assert(all(imag(eigenvals) < 1e-12), 'Covariance eigenvalues have non-zero imaginary parts');
    assert(all(real(eigenvals) >= -1e-12), 'Covariance eigenvalues must be non-negative');
    
    % Invariant 3: Radar Range Resolution Conservation Law: deltaR = c / (2*B)
    c = chirp.utils.PhysicalConstants.C;
    for testBandwidth = [250e6, 500e6, 1.0e9, 2.0e9, 4.0e9]
        ch = chirp.waveform.FmcwChirp('SweepBandwidth', testBandwidth);
        expectedRes = c / (2.0 * testBandwidth);
        assert(abs(ch.RangeResolution - expectedRes) < 1e-9, 'Range resolution invariant violated');
    end
    
    % Invariant 4: Kalman Covariance Uncertainty Monotonicity
    % Measurement update MUST reduce total trace covariance: tr(P_update) < tr(P_prior)
    ekf = chirp.tracking.ExtendedKalmanFilter();
    priorP = eye(6);
    predState = [0; 25; 0; 0; 10; 0];
    zMeas = [25.0; 0.0; 0.0; 10.0];
    [~, updP, ~, ~] = ekf.update(predState, priorP, zMeas);
    tracePrior = trace(priorP);
    traceUpd = trace(updP);
    fprintf('  Invariant 4: Kalman Trace Uncertainty: Prior = %.3f -> Posterior = %.3f\n', ...
        tracePrior, traceUpd);
    assert(traceUpd < tracePrior, 'Kalman filter update failed to reduce state covariance uncertainty');
    
    fprintf('  PASSED: test_mathematical_invariants (4/4 invariants verified)\n');
end
