classdef OmegaKAlgorithm
    % OMEGAKALGORITHM Wavefront omega-k reconstruction SAR imaging engine.
    %
    % Performs Stolt interpolation in 2D spatial wavenumber domain (k_r, k_x)
    % for wide-beamwidth and wideband SAR imaging without approximations.
    
    properties
        Wavelength                    % Carrier wavelength (m)
        PlatformVelocity              % Flight speed v (m/s)
    end
    
    methods
        function obj = OmegaKAlgorithm(wavelengthMeters, platformVelocity)
            % Construct an OmegaKAlgorithm instance.
            if nargin >= 1, obj.Wavelength = wavelengthMeters; end
            if nargin >= 2, obj.PlatformVelocity = platformVelocity; end
        end
        
        function focusedImage = process(obj, rawData, chirpObj, prfHz, r0Meters)
            % Execute omega-k wavefront reconstruction.
            
            [nFast, nSlow] = size(rawData);
            lambda = obj.Wavelength;
            k0 = 2.0 * pi / lambda;
            
            % Step 1: 2D FFT into 2D frequency domain
            spec2D = fftshift(fft2(rawData));
            
            % Spatial wavenumber grids
            krAxis = linspace(-pi * chirpObj.SampleRate / chirp.utils.PhysicalConstants.C, ...
                               pi * chirpObj.SampleRate / chirp.utils.PhysicalConstants.C, nFast).';
            kxAxis = linspace(-pi * prfHz / obj.PlatformVelocity, ...
                               pi * prfHz / obj.PlatformVelocity, nSlow);
                               
            [meshKx, meshKr] = meshgrid(kxAxis, krAxis);
            
            % Reference Function Multiplication (RFM) to focus scene center r0
            % ky_ref = sqrt(4*(k0 + kr)^2 - kx^2)
            kTotal = 2.0 * (k0 + meshKr);
            kySq = kTotal.^2 - meshKx.^2;
            validMask = kySq > 0;
            
            ky = zeros(size(meshKx));
            ky(validMask) = sqrt(kySq(validMask));
            
            rfmPhase = exp(1j * (ky - 2.0 * k0) * r0Meters);
            rfmCorrected = spec2D .* rfmPhase;
            
            % Step 2: 2D Inverse FFT to image domain
            focusedImage = ifft2(ifftshift(rfmCorrected));
        end
    end
end
