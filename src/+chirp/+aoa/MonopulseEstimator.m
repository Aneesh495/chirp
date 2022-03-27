classdef MonopulseEstimator
    % MONOPULSEESTIMATOR Single-pulse amplitude and phase comparison monopulse angle estimation.
    %
    % Forms simultaneous sum (Sigma) and difference (Delta) beam patterns to estimate
    % target off-boresight angle with high precision from a single snapshot.
    
    properties
        MonopulseSlope = 10.0;        % Error slope km for phase/amplitude ratio
        BeamwidthDeg = 15.0;          % Linear beamwidth envelope
    end
    
    methods
        function obj = MonopulseEstimator(monopulseSlope, beamwidthDeg)
            % Construct MonopulseEstimator instance.
            if nargin >= 1 && ~isempty(monopulseSlope)
                obj.MonopulseSlope = monopulseSlope;
            end
            if nargin >= 2 && ~isempty(beamwidthDeg)
                obj.BeamwidthDeg = beamwidthDeg;
            end
        end
        
        function offBoresightDeg = estimate(obj, channelA, channelB)
            % Estimate off-boresight angle from two squinted sub-array antenna channels.
            %
            % Inputs:
            %   channelA - Complex received signal at left sub-array
            %   channelB - Complex received signal at right sub-array
            %
            % Output:
            %   offBoresightDeg - Estimated angle in degrees relative to boresight
            
            sigmaBeam = channelA + channelB;
            deltaBeam = channelA - channelB;
            
            denom = abs(sigmaBeam).^2;
            if denom <= 1e-12
                offBoresightDeg = 0.0;
                return;
            end
            
            rReal = real(deltaBeam .* conj(sigmaBeam)) ./ denom;
            rImag = imag(deltaBeam .* conj(sigmaBeam)) ./ denom;
            
            % Select dominant component (amplitude monopulse = real, phase monopulse = imag)
            if abs(rReal) >= abs(rImag)
                ratio = rReal;
            else
                ratio = -rImag;
            end
            
            % Angle in degrees: for small angles theta_rad = ratio / km
            angleRad = ratio / obj.MonopulseSlope;
            offBoresightDeg = angleRad * (180.0 / pi);
            
            maxAngle = obj.BeamwidthDeg * 0.75;
            offBoresightDeg = max(-maxAngle, min(maxAngle, offBoresightDeg));
        end
    end
end
