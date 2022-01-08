classdef SubcellPeakInterpolator
    % SUBCELLPEAKINTERPOLATOR 2D parabolic and Gaussian sub-bin peak localization.
    %
    % Interpolates discrete integer peak indices (rIdx, vIdx) in Range-Doppler
    % matrices to fractional continuous coordinates (rEst, vEst) with sub-bin precision.
    
    methods (Static)
        function [rFine, vFine, peakMag] = interpolate2D(matrixPower, rIdx, vIdx, rangeAxis, velAxis, useGaussian)
            % Compute continuous coordinates for discrete detection peak.
            %
            % Inputs:
            %   matrixPower - [N_range x N_doppler] power matrix
            %   rIdx        - Discrete integer row index (range bin)
            %   vIdx        - Discrete integer column index (velocity bin)
            %   rangeAxis   - 1 x N_range physical range coordinate vector (m)
            %   velAxis     - 1 x N_doppler physical velocity coordinate vector (m/s)
            %   useGaussian - Boolean (true = Gaussian interpolation on log, false = parabolic)
            %
            % Returns:
            %   rFine       - Continuous physical range estimate (m)
            %   vFine       - Continuous physical radial velocity estimate (m/s)
            %   peakMag     - Interpolated peak magnitude
            
            if nargin < 6 || isempty(useGaussian), useGaussian = false; end
            
            [nR, nV] = size(matrixPower);
            
            % Boundary check
            if rIdx <= 1 || rIdx >= nR || vIdx <= 1 || vIdx >= nV
                rFine = rangeAxis(rIdx);
                vFine = velAxis(vIdx);
                peakMag = matrixPower(rIdx, vIdx);
                return;
            end
            
            % Range direction 3-point stencil
            yR_left = matrixPower(rIdx - 1, vIdx);
            yR_mid = matrixPower(rIdx, vIdx);
            yR_right = matrixPower(rIdx + 1, vIdx);
            
            % Velocity direction 3-point stencil
            yV_left = matrixPower(rIdx, vIdx - 1);
            yV_mid = matrixPower(rIdx, vIdx);
            yV_right = matrixPower(rIdx, vIdx + 1);
            
            if useGaussian
                [deltaR, ~] = chirp.utils.FastMath.gaussian_peak(yR_left, yR_mid, yR_right);
                [deltaV, peakMag] = chirp.utils.FastMath.gaussian_peak(yV_left, yV_mid, yV_right);
            else
                [deltaR, ~] = chirp.utils.FastMath.quadratic_peak(yR_left, yR_mid, yR_right);
                [deltaV, peakMag] = chirp.utils.FastMath.quadratic_peak(yV_left, yV_mid, yV_right);
            end
            
            % Physical spacing
            dr = rangeAxis(2) - rangeAxis(1);
            dv = velAxis(2) - velAxis(1);
            
            rFine = rangeAxis(rIdx) + deltaR * dr;
            vFine = velAxis(vIdx) + deltaV * dv;
        end
    end
end
