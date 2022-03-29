function test_sar()
    % TEST_SAR Unit tests for Synthetic Aperture Radar reconstruction algorithms.
    
    fprintf('Running test_sar...\n');
    
    % 1. SAR geometry kinematics
    geo = chirp.sar.SarGeometry('FlightAltitude', 400.0, 'PlatformVelocity', 60.0, 'ApertureLength', 30.0);
    pos = geo.platformPosition(0.0);
    assert(abs(pos(1) - 0.0) < 1e-6 && abs(pos(3) - 400.0) < 1e-6, 'Platform position mismatch');
    
    rSlant = geo.slantRange(pos, [0; 300; 0]);
    % Distance sqrt(300^2 + 400^2) = 500.0 m
    assert(abs(rSlant - 500.0) < 1e-6, 'Pythagorean 3-4-5 slant range mismatch');
    
    % 2. Backprojection reconstruction focal accuracy
    lambda = 0.03; % 10 GHz
    nPulses = 32;
    uCoords = linspace(-15.0, 15.0, nPulses);
    platPos = zeros(3, nPulses);
    for k = 1:nPulses, platPos(:, k) = geo.platformPosition(uCoords(k)); end
    
    tgtGround = [0.0; 600.0; 0.0];
    nRange = 64;
    rAxis = linspace(500.0, 900.0, nRange);
    dr = rAxis(2) - rAxis(1);
    rProfiles = zeros(nRange, nPulses);
    
    for k = 1:nPulses
        R = norm(platPos(:, k) - tgtGround);
        b = round((R - rAxis(1)) / dr) + 1;
        if b >= 1 && b <= nRange
            rProfiles(b, k) = exp(-1j * (4.0 * pi / lambda) * R);
        end
    end
    
    bp = chirp.sar.BackprojectionEngine(lambda, false);
    gX = linspace(-10.0, 10.0, 21);
    gY = linspace(590.0, 610.0, 21);
    sarImg = bp.reconstruct2D(rProfiles, rAxis, platPos, gX, gY, 0.0);
    
    [~, maxIdx] = max(abs(sarImg(:)));
    [pkY, pkX] = ind2sub(size(sarImg), maxIdx);
    fprintf('  SAR Backprojection focal peak: x = %.2f m, y = %.2f m (True: [0.0, 600.0])\n', ...
        gX(pkX), gY(pkY));
    assert(abs(gX(pkX) - 0.0) <= 2.0, 'SAR focal error in cross-range');
    assert(abs(gY(pkY) - 600.0) <= 2.0, 'SAR focal error in ground-range');
    
    fprintf('  PASSED: test_sar\n');
end
