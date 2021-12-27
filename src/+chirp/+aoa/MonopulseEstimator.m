classdef MonopulseEstimator
    % MONOPULSEESTIMATOR Single-pulse amplitude and phase comparison monopulse angle estimation.
    %
    % Forms simultaneous sum (Sigma) and difference (Delta) beam patterns to estimate
    % target off-boresight angle with high precision from a single snapshot:
    %   theta_est = Re(Delta / Sigma) / km
    
    properties
        MonopulseSlope = 1.6;         % Empirical/analytical error slope km (rad^(-1))
        BeamwidthDeg = 15.0;          % 3-dB beamwidth of individual sub-array beams
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
            
            % Sum and Difference channel synthesis
            sigmaBeam = channelA + channelB;
            deltaBeam = channelA - channelB;
            
            % Monopulse ratio: d_over_s = Re(Delta * conj(Sigma)) / |Sigma|^2
            denom = abs(sigmaBeam).^2;
            if denom <= 1e-12
                offBoresightDeg = 0.0;
                return;
            end
            
            ratio = real(deltaBeam .* conj(sigmaBeam)) ./ denom;
            
            % Angle in radians
            angleRad = ratio / obj.MonopulseSlope;
            offBoresightDeg = angleRad * (180.0 / pi);
            
            % Clamp within valid linear monopulse beamwidth envelope
            maxAngle = obj.BeamwidthDeg * 0.75;
            offBoresightDeg = max(-maxAngle, min(maxAngle, offBoresightDeg));
        end
    end
end
