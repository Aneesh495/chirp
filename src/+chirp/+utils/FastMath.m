classdef FastMath
    % FASTMATH High-performance numerical algorithms and signal interpolation routines.
    %
    % Provides sub-bin interpolation, phase wrapping, Euler 3D rotations,
    % Whittaker-Shannon sinc interpolation, and discrete geometric operations.
    
    methods (Static)
        function p = nextpow2_fast(n)
            % Compute the exponent p such that 2^p >= abs(n).
            if n == 0
                p = 0;
                return;
            end
            p = ceil(log2(abs(n)));
        end
        
        function ang = wrap_to_pi(angleRad)
            % Wrap angular values in radians to the interval [-pi, pi).
            ang = mod(angleRad + pi, 2.0 * pi) - pi;
        end
        
        function ang = wrap_to_2pi(angleRad)
            % Wrap angular values in radians to the interval [0, 2*pi).
            ang = mod(angleRad, 2.0 * pi);
        end
        
        function [delta, peakVal] = quadratic_peak(yLeft, yMid, yRight)
            % 3-point quadratic (parabolic) peak interpolation.
            %
            % Fits a parabola through three equidistant points (-1, yLeft), (0, yMid), (1, yRight).
            % Peak location offset delta is in [-0.5, 0.5] relative to center bin.
            %
            % Returns:
            %   delta   - Fractional bin offset from the middle bin
            %   peakVal - Estimated interpolated peak magnitude
            
            denom = 2.0 * (yLeft - 2.0 * yMid + yRight);
            if abs(denom) < 1e-12
                delta = 0.0;
                peakVal = yMid;
                return;
            end
            
            delta = (yLeft - yRight) / denom;
            delta = max(-0.5, min(0.5, delta));
            peakVal = yMid - 0.25 * (yLeft - yRight) * delta;
        end
        
        function [delta, peakVal] = gaussian_peak(yLeft, yMid, yRight)
            % 3-point Gaussian peak interpolation.
            %
            % Assumes peak follows a Gaussian bell curve. Takes natural logarithm of
            % positive magnitudes and applies quadratic interpolation on log scale.
            
            yLeftSafe = max(yLeft, eps);
            yMidSafe = max(yMid, eps);
            yRightSafe = max(yRight, eps);
            
            lLeft = log(yLeftSafe);
            lMid = log(yMidSafe);
            lRight = log(yRightSafe);
            
            [delta, ~] = chirp.utils.FastMath.quadratic_peak(lLeft, lMid, lRight);
            % Interpolated value on linear scale
            peakVal = yMidSafe * exp(-0.5 * (lLeft - lRight) * delta);
        end
        
        function yInterp = sinc_interp(x, tOrig, tNew)
            % Bandlimited Whittaker-Shannon sinc interpolation.
            %
            % Inputs:
            %   x      - 1D original signal vector
            %   tOrig  - Original sample time indices or coordinates
            %   tNew   - Desired target coordinates for resampling
            
            x = x(:).';
            tOrig = tOrig(:);
            tNew = tNew(:);
            
            dt = tOrig(2) - tOrig(1);
            if dt <= 0
                error('FastMath:InvalidGrid', 'Sampling coordinates must be strictly monotonically increasing.');
            end
            
            % Compute normalized distance matrix (N_new x N_orig)
            diffGrid = (tNew - tOrig.') ./ dt;
            
            % Vectorized sinc computation: sin(pi*u) / (pi*u)
            sincWeights = sin(pi * diffGrid) ./ (pi * diffGrid);
            sincWeights(diffGrid == 0) = 1.0;
            
            % Apply window to suppress truncation sidelobes
            numOrig = length(tOrig);
            win = 0.54 - 0.46 * cos(2.0 * pi * (0:numOrig-1).' / (numOrig - 1));
            sincWeights = sincWeights .* (win.');
            
            yInterp = (sincWeights * x.').';
        end
        
        function R = rot3d(rollDeg, pitchDeg, yawDeg)
            % Generate 3D directional cosine rotation matrix using Z-Y-X Euler angles.
            % Angles given in degrees:
            %   rollDeg  - Rotation about X axis (phi)
            %   pitchDeg - Rotation about Y axis (theta)
            %   yawDeg   - Rotation about Z axis (psi)
            
            phi = rollDeg * (pi / 180.0);
            theta = pitchDeg * (pi / 180.0);
            psi = yawDeg * (pi / 180.0);
            
            Rx = [1, 0, 0; 0, cos(phi), -sin(phi); 0, sin(phi), cos(phi)];
            Ry = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
            Rz = [cos(psi), -sin(psi), 0; sin(psi), cos(psi), 0; 0, 0, 1];
            
            R = Rz * Ry * Rx;
        end
    end
end
