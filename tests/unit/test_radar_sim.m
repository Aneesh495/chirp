function test_radar_sim()
    % TEST_RADAR_SIM Unit tests for radar simulation, physics, and MIMO arrays.
    
    fprintf('Running test_radar_sim...\n');
    
    % 1. Link budget and noise floor
    cfg = chirp.sim.RadarConfig('TxPowerDbm', 12.0, 'NoiseFigureDb', 10.0);
    assert(abs(cfg.TxPowerWatts - 0.0158489) < 1e-6, 'Tx power in watts mismatch');
    assert(cfg.NoiseFloorWatts > 0, 'Noise floor must be positive');
    
    % 2. Radar range equation
    lambda = 0.00389; % 77 GHz
    prx = cfg.receivedPower(20.0, 5.0, lambda);
    assert(prx > 0, 'Received power must be positive');
    % Inverse fourth power test: doubling range reduces power by 16x (12 dB)
    prxDouble = cfg.receivedPower(40.0, 5.0, lambda);
    ratio = prx / prxDouble;
    assert(abs(ratio - 16.0) < 1e-4, 'Inverse 4th power law violated');
    
    % 3. Target kinematics propagation
    tgt = chirp.sim.TargetModel('Position', [0; 10; 0], 'Velocity', [2; 5; 0], 'Acceleration', [0; -1; 0]);
    tgtProp = tgt.propagate(2.0);
    % Expected pos: [0 + 2*2; 10 + 5*2 - 0.5*1*4; 0] = [4; 18; 0]
    assert(abs(tgtProp.Position(1) - 4.0) < 1e-6, 'Propagated x mismatch');
    assert(abs(tgtProp.Position(2) - 18.0) < 1e-6, 'Propagated y mismatch');
    
    % 4. Swerling models sampling
    s0 = chirp.sim.SwerlingModel.sample(0, 10.0, 100);
    assert(all(s0 == 10.0), 'Swerling 0 must be constant');
    
    s1 = chirp.sim.SwerlingModel.sample(1, 10.0, 100);
    assert(all(s1 == s1(1)), 'Swerling 1 must be constant across frame');
    
    s2 = chirp.sim.SwerlingModel.sample(2, 10.0, 100);
    assert(std(s2) > 0.5, 'Swerling 2 must fluctuate chirp-to-chirp');
    
    % 5. MIMO array geometry
    mimo = chirp.sim.MimoArray.automotive3Tx4Rx(lambda);
    assert(mimo.NumTx == 3, 'NumTx mismatch');
    assert(mimo.NumRx == 4, 'NumRx mismatch');
    assert(mimo.NumVirtual == 12, 'Virtual channels mismatch');
    
    sv = mimo.steeringVector(0, 0, lambda);
    assert(length(sv) == 12, 'Steering vector length mismatch');
    assert(abs(norm(sv) - sqrt(12)) < 1e-6, 'Steering vector norm mismatch');
    
    fprintf('  PASSED: test_radar_sim\n');
end
