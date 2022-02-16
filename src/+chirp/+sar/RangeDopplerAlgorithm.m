classdef RangeDopplerAlgorithm
    % RANGEDOPPLERALGORITHM Classical frequency-domain Range-Doppler SAR processor.
    %
    % Performs Range Cell Migration Correction (RCMC) and azimuth matched filtering
    % to reconstruct focused SAR imagery efficiently via Fast Fourier Transforms.
    
    properties
        Wavelength                    % Carrier wavelength (m)
        PlatformVelocity              % Flight speed v (m/s)
    end
    
    methods
        function obj = RangeDopplerAlgorithm(wavelengthMeters, platformVelocity)
            % Construct a RangeDopplerAlgorithm instance.
            if nargin >= 1, obj.Wavelength = wavelengthMeters; end
            if nargin >= 2, obj.PlatformVelocity = platformVelocity; end
        end
        
        function focusedImage = process(obj, rawEchoData, chirpObj, prfHz, r0Meters)
            % Execute Range-Doppler Algorithm on raw 2D SAR pulse data.
            %
            % Inputs:
            %   rawEchoData      - [N_fast x N_slow] raw de-chirped beat data
            %   chirpObj         - chirp.waveform.FmcwChirp configuration
            %   prfHz            - Pulse repetition frequency along azimuth (Hz)
            %   r0Meters         - Reference slant range to scene center (m)
            %
            % Output:
            %   focusedImage     - [N_fast x N_slow] 2D focused complex SAR image
            
            [nFast, nSlow] = size(rawEchoData);
            lambda = obj.Wavelength;
            v = obj.PlatformVelocity;
            
            % Step 1: Range compression (Range-FFT along dimension 1)
            rangeCompressed = fft(rawEchoData, nFast, 1);
            
            % Step 2: Azimuth FFT to transform into Range-Doppler domain (r, fa)
            rdData = fftshift(fft(rangeCompressed, nSlow, 2), 2);
            
            % Step 3: Range Cell Migration Correction (RCMC)
            % Azimuth frequency axis
            if mod(nSlow, 2) == 0
                faIndices = (-nSlow/2):(nSlow/2 - 1);
            else
                faIndices = (-(nSlow - 1)/2):((nSlow - 1)/2);
            end
            faAxis = faIndices * (prfHz / nSlow);
            
            % Range migration curve: delta_R(fa) = r0 * (1 / sqrt(1 - (lambda * fa / (2*v))^2) - 1)
            sinTheta = (lambda * faAxis) / (2.0 * v);
            validAngle = abs(sinTheta) < 0.95;
            
            rcmcData = zeros(nFast, nSlow);
            for iSlow = 1:nSlow
                if ~validAngle(iSlow)
                    continue;
                end
                deltaR = r0Meters * (1.0 / sqrt(1.0 - sinTheta(iSlow)^2) - 1.0);
                
                % Fast-time phase ramp to shift range in frequency domain
                binShift = deltaR / chirpObj.RangeResolution;
                shiftPhase = exp(-1j * 2.0 * pi * binShift * (0:(nFast - 1)).' / nFast);
                
                rcmcData(:, iSlow) = rdData(:, iSlow) .* shiftPhase;
            end
            
            % Step 4: Azimuth Matched Filtering
            % Phase reference: H_az(fa) = exp(j * 4*pi*r0 / lambda * sqrt(1 - (lambda*fa / (2*v))^2))
            azMatchedFilter = exp(1j * (4.0 * pi * r0Meters / lambda) * sqrt(max(1.0 - sinTheta.^2, 0.0)));
            filteredRD = rcmcData .* azMatchedFilter;
            
            % Step 5: Azimuth IFFT back to time domain
            focusedImage = ifft(ifftshift(filteredRD, 2), nSlow, 2);
        end
    end
end
