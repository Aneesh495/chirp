function test_mex_parity()
    % TEST_MEX_PARITY Verify exact numerical agreement between C MEX kernels and MATLAB code.
    
    fprintf('Running test_mex_parity...\n');
    
    % Check if MEX kernels are compiled
    if exist('chirp_cfar_2d_mex', 'file') ~= 3
        warning('chirp:mex:missing', 'MEX kernels not compiled in bin/mex. Skipping test_mex_parity.');
        return;
    end
    
    % 1. Parity test: 2D CA-CFAR
    nR = 64; nV = 32;
    rdmPwr = -log(max(rand(nR, nV), 1e-12));
    rdmPwr(30, 15) = 100.0; % Strong target
    
    % Pure MATLAB implementation
    caMat = chirp.cfar.CaCfar(1e-4, 6, 2, 3, 1);
    [maskMat, threshMat, nfMat] = caMat.detect(rdmPwr);
    
    % C MEX implementation
    [maskMex, threshMex, nfMex] = chirp_cfar_2d_mex(rdmPwr, 1e-4, 6, 2, 3, 1);
    
    % Verify bitwise identity of logical detection mask
    mismatches = sum(maskMat(:) ~= maskMex(:));
    assert(mismatches == 0, 'MEX vs MATLAB CFAR detection mask mismatch');
    
    % Verify numerical precision of noise floor matrix
    diffNf = max(abs(nfMat(:) - nfMex(:)));
    fprintf('  CFAR 2D Noise Floor Max Delta: %e\n', diffNf);
    assert(diffNf < 1e-10, 'MEX vs MATLAB CFAR noise floor numerical deviation');
    
    % 2. Parity test: 1D Range-FFT
    nAdc = 64; nChirps = 16; nVirt = 4;
    cube = randn(nAdc, nChirps, nVirt) + 1j * randn(nAdc, nChirps, nVirt);
    win = chirp.dsp.WindowFunction.create('hann', nAdc);
    
    % Native FFT
    winMat = repmat(win, [1, nChirps, nVirt]);
    fftNative = fft(cube .* winMat, nAdc, 1);
    
    % C MEX FFT
    fftMex = chirp_range_fft_mex(cube, win, nAdc);
    
    diffFft = max(abs(fftNative(:) - fftMex(:)));
    fprintf('  Range-FFT Max Delta vs Native FFT: %e\n', diffFft);
    assert(diffFft < 1e-10, 'MEX vs Native FFT numerical deviation');
    
    fprintf('  PASSED: test_mex_parity\n');
end
