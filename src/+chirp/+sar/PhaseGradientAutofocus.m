classdef PhaseGradientAutofocus
    % PHASEGRADIENTAUTOFOCUS Non-parametric SAR autofocus motion error estimation.
    %
    % Corrects unmeasured platform flight perturbations and atmospheric phase jitter
    % by estimating phase error gradients across isolated prominent scatterers.
    
    properties
        NumIterations = 5;            % Maximum autofocus correction iterations
        WindowFraction = 0.5;         % Initial window width fraction
    end
    
    methods
        function obj = PhaseGradientAutofocus(numIterations, windowFraction)
            % Construct a PhaseGradientAutofocus instance.
            if nargin >= 1 && ~isempty(numIterations), obj.NumIterations = numIterations; end
            if nargin >= 2 && ~isempty(windowFraction), obj.WindowFraction = windowFraction; end
        end
        
        function [correctedImage, estimatedPhaseError] = correct(obj, blurredImage)
            % Execute iterative Phase Gradient Autofocus on 2D complex SAR image.
            %
            % Input:
            %   blurredImage        - [N_range x N_azimuth] un-focused SAR image
            % Output:
            %   correctedImage      - Focused sharp image
            %   estimatedPhaseError - Estimated 1 x N_azimuth phase correction profile
            
            [nR, nAz] = size(blurredImage);
            currentImg = blurredImage;
            totalPhaseError = zeros(1, nAz);
            
            for iter = 1:obj.NumIterations
                % Step 1: Range-by-range peak detection and circular shift centering
                shiftedImg = zeros(nR, nAz);
                for r = 1:nR
                    row = currentImg(r, :);
                    [~, pkCol] = max(abs(row));
                    shiftAmount = round(nAz / 2) - pkCol;
                    shiftedImg(r, :) = circshift(row, [0, shiftAmount]);
                end
                
                % Step 2: Windowing around centered prominent scatterer
                winWidth = max(8, round(nAz * (obj.WindowFraction^iter)));
                halfW = floor(winWidth / 2);
                midCol = round(nAz / 2);
                
                colMask = false(1, nAz);
                colMask(max(1, midCol - halfW):min(nAz, midCol + halfW)) = true;
                
                windowedData = shiftedImg;
                windowedData(:, ~colMask) = 0.0;
                
                % Step 3: Phase gradient estimation
                % Transform to range-Doppler domain
                G = fft(windowedData, nAz, 2);
                
                % Derivative along Doppler frequency via circular shift: G_dot
                G_dot = [diff(G, 1, 2), G(:, 1) - G(:, end)];
                
                % Linear unbiased estimator of phase gradient: phi_dot = Im(G* . G_dot) / |G|^2
                numerator = sum(imag(conj(G) .* G_dot), 1);
                denominator = sum(abs(G).^2, 1);
                
                phiDot = numerator ./ max(denominator, 1e-12);
                
                % Step 4: Integrate gradient to obtain phase error
                phaseError = cumsum(phiDot);
                % De-trend constant and linear terms
                phaseError = detrend(phaseError);
                
                totalPhaseError = totalPhaseError + phaseError;
                
                % Apply phase correction across columns of range-Doppler image
                imgRD = fft(currentImg, nAz, 2);
                corrRD = imgRD .* exp(-1j * phaseError);
                currentImg = ifft(corrRD, nAz, 2);
            end
            
            correctedImage = currentImg;
            estimatedPhaseError = totalPhaseError;
        end
    end
end
