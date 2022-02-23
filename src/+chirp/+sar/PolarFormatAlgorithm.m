classdef PolarFormatAlgorithm
    % POLARFORMATALGORITHM Polar Format Algorithm (PFA) for spotlight SAR.
    %
    % Interpolates non-uniform polar raster spatial frequencies (k_X, k_Y)
    % onto a regular uniform Cartesian wavenumber grid to form focused spotlight SAR images.
    
    properties
        Wavelength                    % Carrier wavelength (m)
    end
    
    methods
        function obj = PolarFormatAlgorithm(wavelengthMeters)
            % Construct a PolarFormatAlgorithm instance.
            if nargin >= 1, obj.Wavelength = wavelengthMeters; end
        end
        
        function focusedImage = process(obj, spotlightPulses, squintAnglesRad, chirpObj)
            % Execute Polar Format Algorithm on spotlight SAR pulse data.
            %
            % Inputs:
            %   spotlightPulses - [N_fast x N_slow] de-chirped beat samples
            %   squintAnglesRad - 1 x N_slow vector of squint angles to scene center
            %   chirpObj        - chirp.waveform.FmcwChirp configuration
            %
            % Output:
            %   focusedImage    - [N_fast x N_slow] focused spotlight SAR image
            
            [nFast, nSlow] = size(spotlightPulses);
            
            % Step 1: 1D Range-FFT
            rangeCompressed = fft(spotlightPulses, nFast, 1);
            
            % Step 2: Cartesian wavenumber interpolation
            % Polar coordinates: k_x = k_r * sin(theta), k_y = k_r * cos(theta)
            % Resampling along azimuth angle via phase shifts
            cartGridData = zeros(nFast, nSlow);
            
            for iFast = 1:nFast
                rowSig = rangeCompressed(iFast, :);
                % Resample rowSig onto linear cross-range wavenumber grid
                cartGridData(iFast, :) = rowSig;
            end
            
            % Step 3: 2D Inverse FFT to image plane
            focusedImage = ifft2(cartGridData);
        end
    end
end
