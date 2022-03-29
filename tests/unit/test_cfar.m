function test_cfar()
    % TEST_CFAR Unit tests for adaptive Constant False Alarm Rate detectors.
    
    fprintf('Running test_cfar...\n');
    
    % 1. Threshold scaling factor alpha derivations
    alphaCA = chirp.cfar.ThresholdCalculator.caCfarAlpha(1e-4, 32);
    % Analytical: 32 * (10000^(1/32) - 1) ~ 10.655
    assert(abs(alphaCA - 10.655) < 0.1, 'CA-CFAR alpha derivation error');
    
    alphaOS = chirp.cfar.ThresholdCalculator.osCfarAlpha(1e-4, 32, 24);
    assert(alphaOS > 1.0 && alphaOS < 100.0, 'OS-CFAR alpha out of expected range');
    
    % 2. 1D CFAR detector
    noiseVec = -log(max(rand(1, 200), 1e-12));
    noiseVec(50) = 40.0; % Injected strong peak
    
    cfar1d = chirp.cfar.CfarDetector1D('Mode', 'CA', 'Pfa', 1e-4, 'NumTrain', 16, 'NumGuard', 2);
    [mask1d, ~, ~] = cfar1d.detect(noiseVec);
    assert(mask1d(50) == true, '1D CFAR target not detected');
    
    % 3. 2D CA-CFAR with integral image acceleration
    nR = 64; nV = 32;
    rdmPwr = -log(max(rand(nR, nV), 1e-12));
    rdmPwr(25, 12) = 50.0;
    
    ca2d = chirp.cfar.CaCfar(1e-4, 8, 2, 4, 1);
    [mask2d, thresh2d, nf2d] = ca2d.detect(rdmPwr);
    assert(mask2d(25, 12) == true, '2D CA-CFAR target not detected');
    
    % 4. 2D OS-CFAR detector
    os2d = chirp.cfar.OsCfar(1e-4, 6, 2, 3, 1);
    [maskOS, ~, ~] = os2d.detect(rdmPwr);
    assert(maskOS(25, 12) == true, '2D OS-CFAR target not detected');
    
    % 5. Subcell Peak Interpolator
    rAxis = linspace(0, 50, nR);
    vAxis = linspace(-15, 15, nV);
    [rFine, vFine, ~] = chirp.cfar.SubcellPeakInterpolator.interpolate2D(rdmPwr, 25, 12, rAxis, vAxis);
    assert(abs(rFine - rAxis(25)) < (rAxis(2) - rAxis(1)), 'Subcell range offset too large');
    assert(abs(vFine - vAxis(12)) < (vAxis(2) - vAxis(1)), 'Subcell velocity offset too large');
    
    fprintf('  PASSED: test_cfar\n');
end
