classdef RootMusicEstimator
    % ROOTMUSICESTIMATOR Polynomial root-finding MUSIC for uniform linear arrays (ULA).
    %
    % Converts spatial noise subspace projection into an algebraic polynomial on the
    % complex z-plane. Finds roots closest to the unit circle and directly computes
    % arrival angles without spatial grid scanning.
    
    properties
        ElementSpacingMeters          % Inter-element spacing d along ULA (m)
        Wavelength                    % Carrier wavelength (m)
    end
    
    methods
        function obj = RootMusicEstimator(elementSpacingMeters, wavelengthMeters)
            % Construct RootMusicEstimator for a uniform linear array.
            if nargin >= 2
                obj.ElementSpacingMeters = elementSpacingMeters;
                obj.Wavelength = wavelengthMeters;
            end
        end
        
        function anglesDeg = estimate(obj, spatialSnapshots, numSignals)
            % Estimate angles of arrival directly via polynomial root finding.
            %
            % Inputs:
            %   spatialSnapshots - [M x N_snapshots] or [M x M] covariance matrix
            %   numSignals       - Number of source targets to resolve
            %
            % Output:
            %   anglesDeg        - 1 x numSignals estimated arrival angles in degrees
            
            if size(spatialSnapshots, 1) == size(spatialSnapshots, 2)
                R = spatialSnapshots;
            else
                nSnap = size(spatialSnapshots, 2);
                R = (spatialSnapshots * spatialSnapshots') / max(nSnap, 1);
            end
            
            % Forward-backward spatial smoothing for rank restoration
            R = chirp.utils.MatrixUtils.forward_backward_avg(R);
            
            M = size(R, 1);
            [~, En, ~] = chirp.utils.MatrixUtils.subspace_decompose(R, numSignals);
            
            % Compute matrix C = En * En'
            C = En * En';
            
            % Sum along diagonals of C to obtain polynomial coefficients
            % D(z) = sum_{m=-(M-1)}^{M-1} c_m * z^(-m)
            polyCoeffs = zeros(1, 2 * M - 1);
            for m = -(M - 1):(M - 1)
                idx = m + M;
                polyCoeffs(idx) = sum(diag(C, m));
            end
            
            % Find polynomial roots
            allRoots = roots(polyCoeffs);
            
            % Select roots strictly inside unit circle: |z| < 1
            insideMask = abs(allRoots) < 1.0;
            insideRoots = allRoots(insideMask);
            
            if isempty(insideRoots)
                anglesDeg = zeros(1, numSignals);
                return;
            end
            
            % Sort roots by proximity to the unit circle: 1 - |z|
            [~, sortIdx] = sort(abs(1.0 - abs(insideRoots)), 'ascend');
            selectedRoots = insideRoots(sortIdx);
            
            k = min(numSignals, length(selectedRoots));
            targetRoots = selectedRoots(1:k);
            
            % Angle extraction:
            % For steering vector a_m = exp(1j * m * 2*pi * d * sin(theta) / lambda),
            % the conjugate root phase yields: sin(theta) = -arg(z) * lambda / (2*pi * d)
            phases = angle(targetRoots);
            sinTheta = -(phases * obj.Wavelength) / (2.0 * pi * obj.ElementSpacingMeters);
            
            % Clamp to [-1, 1] for real arcsin
            sinTheta = max(-1.0, min(1.0, sinTheta));
            anglesDeg = asin(sinTheta) * (180.0 / pi);
            anglesDeg = sort(anglesDeg(:).', 'ascend');
        end
    end
end
