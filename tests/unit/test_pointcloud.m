function test_pointcloud()
    % TEST_POINTCLOUD Unit tests for 4D radar point clouds and spatial clustering.
    
    fprintf('Running test_pointcloud...\n');
    
    % 1. PointCloud4D construction and cropping
    pos = [1.0, 2.0, 10.0; ...
           5.0, 6.0, 30.0; ...
           0.0, 0.5,  2.0];
    vel = [2.0, -5.0, 15.0];
    snr = [25.0, 20.0, 15.0];
    rcs = [1.0, 2.0, 0.5];
    cloud = chirp.pointcloud.PointCloud4D(pos, vel, snr, rcs, 0.0);
    assert(cloud.NumPoints == 3, 'Point count mismatch');
    
    cropped = cloud.cropBox(-5, 5, 0, 15, -1, 5);
    assert(cropped.NumPoints == 2, 'Crop box count mismatch');
    
    % 2. Rigid transformation
    R = chirp.utils.FastMath.rot3d(0, 0, 90); % 90 deg yaw
    t = [0; 0; 0];
    transformed = cloud.transform(R, t);
    % Point [1; 5; 0] rotated 90 deg around z becomes [-5; 1; 0]
    assert(abs(transformed.Positions(1, 1) - (-5.0)) < 1e-6, 'Rotated X mismatch');
    assert(abs(transformed.Positions(2, 1) - 1.0) < 1e-6, 'Rotated Y mismatch');
    
    % 3. SpatialKdTree radius and k-NN searches
    pts = [0.0, 1.0, 2.0, 10.0; ...
           0.0, 0.0, 0.0, 10.0];
    kd = chirp.pointcloud.SpatialKdTree(pts);
    neighbors = kd.queryRadius([0.0; 0.0], 1.5);
    % Points 1 and 2 are within 1.5 of [0, 0]
    assert(length(neighbors) == 2, 'KdTree radius neighbor count mismatch');
    
    [knnIdx, knnDist] = kd.queryKnn([0.0; 0.0], 3);
    assert(length(knnIdx) == 3, 'KdTree k-NN count mismatch');
    assert(knnDist(1) == 0.0, 'First neighbor must be exact match at dist 0');
    
    % 4. VelocityAwareDbscan clustering
    % Cluster 1: 8 points around (0, 10, 0) with v = 5 m/s
    c1Pos = [0.0, 0.1, -0.1, 0.05, -0.05, 0.08, -0.08, 0.02; ...
             10.0 + [0.0, 0.1, -0.1, 0.05, -0.05, 0.08, -0.08, 0.02]; ...
             zeros(1, 8)];
    c1Vel = 5.0 * ones(1, 8);
    
    % Cluster 2: 8 points around (0, 10, 0) with v = -10 m/s (Doppler separated obstacle!)
    c2Pos = c1Pos;
    c2Vel = -10.0 * ones(1, 8);
    
    bothPos = [c1Pos, c2Pos];
    bothVel = [c1Vel, c2Vel];
    cloudClust = chirp.pointcloud.PointCloud4D(bothPos, bothVel, ones(1, 16), ones(1, 16), 0);
    
    dbscan = chirp.pointcloud.VelocityAwareDbscan(1.5, 2.0, 3);
    [labels, nClusters] = dbscan.cluster(cloudClust);
    % Even though positions overlap at (0, 10), velocity awareness MUST separate them into 2 distinct clusters!
    fprintf('  DBSCAN separated clusters: %d (Expected: 2)\n', nClusters);
    assert(nClusters == 2, 'Doppler-aware DBSCAN failed to separate overlapping velocity modes');
    
    % 5. BoundingBoxEstimator
    bbox = chirp.pointcloud.BoundingBoxEstimator.fitCluster(c1Pos, c1Vel);
    assert(abs(bbox.centroid(2) - 10.0) < 0.5, 'Centroid Y mismatch');
    assert(abs(bbox.velocity - 5.0) < 0.5, 'Mean velocity mismatch');
    
    fprintf('  PASSED: test_pointcloud\n');
end
