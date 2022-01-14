classdef OutlierFilter
    % OUTLIERFILTER Statistical Outlier Removal (SOR) for radar point clouds.
    %
    % Removes isolated noise returns, speckle, and multipath ghost spikes by
    % analyzing mean distance distributions to k-nearest spatial neighbors.
    
    properties
        NumNeighbors = 6;             % Number of nearest neighbors k to evaluate
        StdDevMultiplier = 1.8;       % Multiplier alpha for threshold (mean + alpha * std)
    end
    
    methods
        function obj = OutlierFilter(kNeighbors, stdMultiplier)
            % Construct an OutlierFilter instance.
            if nargin >= 1 && ~isempty(kNeighbors), obj.NumNeighbors = kNeighbors; end
            if nargin >= 2 && ~isempty(stdMultiplier), obj.StdDevMultiplier = stdMultiplier; end
        end
        
        function [filteredCloud, inlierMask] = filter(obj, pointCloud)
            % Filter point cloud and retain inlier points.
            %
            % Returns:
            %   filteredCloud - Inlier PointCloud4D instance
            %   inlierMask    - 1 x N boolean vector (true = inlier)
            
            N = pointCloud.NumPoints;
            k = obj.NumNeighbors;
            
            if N <= k
                filteredCloud = pointCloud;
                inlierMask = true(1, N);
                return;
            end
            
            kd = chirp.pointcloud.SpatialKdTree(pointCloud.Positions);
            meanDists = zeros(1, N);
            
            for i = 1:N
                % Query k+1 neighbors (includes point itself at dist 0)
                [~, dists] = kd.queryKnn(pointCloud.Positions(:, i), k + 1);
                % Exclude first distance (0.0)
                meanDists(i) = mean(dists(2:end));
            end
            
            mu = mean(meanDists);
            sigma = std(meanDists);
            thresh = mu + obj.StdDevMultiplier * sigma;
            
            inlierMask = (meanDists <= thresh);
            
            filteredCloud = chirp.pointcloud.PointCloud4D(...
                pointCloud.Positions(:, inlierMask), ...
                pointCloud.Velocities(inlierMask), ...
                pointCloud.SnrDb(inlierMask), ...
                pointCloud.Rcs(inlierMask), ...
                pointCloud.Timestamp ...
            );
        end
    end
end
