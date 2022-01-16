classdef BoundingBoxEstimator
    % BOUNDINGBOXESTIMATOR 3D Oriented Bounding Box (OBB) fitting for radar clusters.
    %
    % Performs Principal Component Analysis (PCA) on clustered point sets to estimate
    % obstacle centroid, 3D dimensions (length, width, height), and heading yaw.
    
    methods (Static)
        function bbox = fitCluster(clusterPositions, clusterVelocities)
            % Compute 3D oriented bounding box for points in a cluster.
            %
            % Inputs:
            %   clusterPositions  - 3 x N matrix of Cartesian points [x; y; z]
            %   clusterVelocities - (Optional) 1 x N vector of radial velocities
            %
            % Returns:
            %   bbox struct with fields:
            %     centroid   - 3 x 1 center [xc; yc; zc] (m)
            %     dimensions - 3 x 1 [length; width; height] (m)
            %     yawDeg     - Principal heading orientation in degrees
            %     velocity   - Mean radial velocity (m/s)
            %     numPoints  - Integer point count N
            
            N = size(clusterPositions, 2);
            if N == 0
                bbox = struct('centroid', [0;0;0], 'dimensions', [0;0;0], 'yawDeg', 0, 'velocity', 0, 'numPoints', 0);
                return;
            end
            
            centroid = mean(clusterPositions, 2);
            
            if nargin >= 2 && ~isempty(clusterVelocities)
                velMean = mean(clusterVelocities);
            else
                velMean = 0.0;
            end
            
            if N < 3
                % Degenerate cluster: axis-aligned fallback
                minPos = min(clusterPositions, [], 2);
                maxPos = max(clusterPositions, [], 2);
                dims = max(maxPos - minPos, 0.2); % Minimum bounding extent
                yaw = 0.0;
            else
                % Center points in XY horizontal plane
                centeredXY = clusterPositions(1:2, :) - centroid(1:2);
                covXY = (centeredXY * centeredXY') / (N - 1);
                
                [V, D] = eig(covXY);
                [~, maxIdx] = max(diag(D));
                principalAxis = V(:, maxIdx);
                
                % Heading yaw angle in horizontal plane
                yaw = atan2(principalAxis(2), principalAxis(1)) * (180.0 / pi);
                
                % Rotate points into principal alignment frame
                cosYaw = cos(yaw * pi / 180.0);
                sinYaw = sin(yaw * pi / 180.0);
                R2d = [cosYaw, sinYaw; -sinYaw, cosYaw];
                
                alignedXY = R2d * centeredXY;
                
                len = max(alignedXY(1, :)) - min(alignedXY(1, :));
                wid = max(alignedXY(2, :)) - min(alignedXY(2, :));
                hgt = max(clusterPositions(3, :)) - min(clusterPositions(3, :));
                
                dims = [max(len, 0.4); max(wid, 0.4); max(hgt, 0.3)];
            end
            
            bbox = struct(...
                'centroid', centroid, ...
                'dimensions', dims, ...
                'yawDeg', yaw, ...
                'velocity', velMean, ...
                'numPoints', N ...
            );
        end
    end
end
