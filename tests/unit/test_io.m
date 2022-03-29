function test_io()
    % TEST_IO Unit tests for hardware binary ingestion, codecs, and exporters.
    
    fprintf('Running test_io...\n');
    
    testDir = 'output/unit_test_io';
    if ~exist(testDir, 'dir'), mkdir(testDir); end
    cleanUpDir = onCleanup(@() rmdir(testDir, 's'));
    
    % 1. RadarCubeCodec bitwise round-trip and CRC32
    origCube = randn(32, 16, 4) + 1j * randn(32, 16, 4);
    rcubeFile = fullfile(testDir, 'test.rcube');
    chirp.io.RadarCubeCodec.save(rcubeFile, origCube);
    
    [loadedCube, dims] = chirp.io.RadarCubeCodec.load(rcubeFile);
    assert(all(dims == [32, 16, 4]), 'Dimensions mismatch');
    diffCube = max(abs(origCube(:) - loadedCube(:)));
    assert(diffCube < 1e-12, 'RadarCubeCodec precision loss');
    
    % Corrupt single byte to test CRC32 failure rejection
    fid = fopen(rcubeFile, 'r+b');
    fseek(fid, 30, 'bof');
    fwrite(fid, 255, 'uint8');
    fclose(fid);
    
    corruptCaught = false;
    try
        chirp.io.RadarCubeCodec.load(rcubeFile);
    catch err
        corruptCaught = true;
    end
    assert(corruptCaught, 'RadarCubeCodec failed to reject CRC32-corrupted file');
    
    % 2. DCA1000 EVM synthetic reader
    dca = chirp.io.Dca1000Reader('NumAdcSamples', 32, 'NumChirpsPerFrame', 16, 'NumRx', 4);
    dcaFile = fullfile(testDir, 'dca_test.bin');
    dca.writeSyntheticBinary(origCube, dcaFile);
    loadedDca = dca.readBinaryFile(dcaFile);
    assert(size(loadedDca, 1) == 32 && size(loadedDca, 2) == 16 && size(loadedDca, 3) == 4, 'DCA shape mismatch');
    
    % 3. PLY writer
    pts = [0.0; 10.0; 0.0];
    cloud = chirp.pointcloud.PointCloud4D(pts, 5.0, 20.0, 1.0, 0.0);
    plyFile = fullfile(testDir, 'test.ply');
    chirp.io.PlyWriter.writeAscii(plyFile, cloud);
    assert(exist(plyFile, 'file') == 2, 'PLY file was not generated');
    
    % 4. LAS writer
    lasFile = fullfile(testDir, 'test.las');
    chirp.io.LasWriter.write(lasFile, cloud);
    assert(exist(lasFile, 'file') == 2, 'LAS file was not generated');
    
    % 5. GeoJSON writer
    geoFile = fullfile(testDir, 'test.geojson');
    chirp.io.GeoJsonWriter.writePointCloud(geoFile, cloud);
    assert(exist(geoFile, 'file') == 2, 'GeoJSON file was not generated');
    
    fprintf('  PASSED: test_io\n');
end
