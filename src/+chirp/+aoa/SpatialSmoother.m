classdef SpatialSmoother
    % SPATIALSMOOTHER Sub-array spatial smoothing for rank restoration of coherent signals.
    %
    % De-correlates fully coherent multipath arrivals by partitioning an M-element ULA
    % into L overlapping sub-arrays of dimension p = M - L + 1 and averaging their covariances.
    
    properties
        SubarraySize                  % Number of elements in each sub-array p
    end
    
    methods
        function obj = SpatialSmoother(subarraySize)
            % Construct a SpatialSmoother with specified sub-array dimension.
            if nargin >= 1
                obj.SubarraySize = subarraySize;
            end
        end
        
        function Rsmooth = smoothCovariance(obj, fullSnapshots)
            % Compute forward spatially smoothed covariance matrix.
            %
            % Input:
            %   fullSnapshots - [M x N_snapshots] raw snapshot data
            % Output:
            %   Rsmooth       - [p x p] smoothed covariance matrix with restored rank
            
            [M, nSnap] = size(fullSnapshots);
            p = obj.SubarraySize;
            
            if isempty(p) || p <= 0 || p > M
                p = ceil(M * 0.75);
                obj.SubarraySize = p;
            end
            
            L = M - p + 1; % Number of overlapping sub-arrays
            Rsmooth = zeros(p, p);
            
            for l = 1:L
                subData = fullSnapshots(l:(l + p - 1), :);
                R_l = (subData * subData') / max(nSnap, 1);
                Rsmooth = Rsmooth + R_l;
            end
            
            Rsmooth = Rsmooth / L;
        end
    end
end
