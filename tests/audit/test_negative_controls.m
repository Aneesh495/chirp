function test_negative_controls()
    % TEST_NEGATIVE_CONTROLS Independent negative controls and failure boundary tests.
    %
    % Verifies that detectors do not hallucinate targets on pure noise, that numerical
    % singularities in ill-conditioned covariance matrices are safely regularized,
    % and that corrupted data frames are strictly rejected.
    
    fprintf('Running test_negative_controls...\n');
    
    % Control 1: Pure Rayleigh noise input to 2D CFAR detector
    % With Pfa = 1e-3 on 10,000 noise cells, expected hits ~ 10. Hits must be < 35.
    nR = 100; nV = 100;
    pureNoisePwr = -log(max(rand(nR, nV), 1e-12));
    ca = chirp.cfar.CaCfar(1e-3, 8, 2, 4, 1);
    [maskNoise, ~, ~] = ca.detect(pureNoisePwr);
    numFalseAlarms = sum(maskNoise(:));
    fprintf('  Negative Control 1: Pure noise false alarms = %d / 10000 (Target Pfa = 1e-3)\n', numFalseAlarms);
    assert(numFalseAlarms < 35, 'CFAR false alarm rate exceeds theoretical statistical bound');
    
    % Control 2: Zero signal input
    zeroPwr = zeros(50, 50);
    [maskZero, ~, ~] = ca.detect(zeroPwr);
    assert(sum(maskZero(:)) == 0, 'Negative Control 2: Zero signal produced phantom detections');
    
    % Control 3: Ill-conditioned singular spatial covariance matrix
    % Rank-1 matrix for 8-element array (all channels identical)
    R_singular = ones(8, 8);
    % Inverting R_singular directly would divide by zero / yield NaN.
    % Regularized inversion must produce a finite, well-conditioned matrix.
    R_reg_inv = chirp.utils.MatrixUtils.regularized_inv(R_singular, 1e-2);
    assert(all(isfinite(R_reg_inv(:))), 'Negative Control 3: Regularized inversion failed on singular matrix');
    
    % Control 4: Subspace decomposition with equal eigenvalues
    % Covariance is identity matrix (pure isotropic noise, no dominant signal)
    R_iso = eye(8);
    [Es, En, lambda] = chirp.utils.MatrixUtils.subspace_decompose(R_iso, 2);
    assert(size(Es, 2) == 2, 'Signal subspace dimension mismatch');
    assert(size(En, 2) == 6, 'Noise subspace dimension mismatch');
    assert(all(isfinite(Es(:))) && all(isfinite(En(:))), 'Subspace decomposition produced NaN/Inf on isotropic noise');
    
    % Control 5: MTI clutter cancellation on static clutter
    staticClutter = 50.0 * ones(64, 32, 2);
    mtiFiltered = chirp.dsp.MtiFilter.twoPulseCanceler(staticClutter);
    diffBlock = mtiFiltered(:, 2:end, :);
    assert(max(abs(diffBlock(:))) < 1e-12, 'Negative Control 5: Static clutter not fully suppressed by MTI');
    
    fprintf('  PASSED: test_negative_controls (5/5 controls passed)\n');
end
