classdef VelocityAwareDbscan
    % VELOCITYAWAREDBSCAN Doppler-penalized spatial clustering for radar point clouds.
    %
    % Groups 4D radar returns into physical obstacle clusters using an ellipsoidal
    % metric that penalizes spatial distance and radial velocity discrepancies:
    %   d(p1, p2) = sqrt(||r1 - r2||^2 / eps_r^2 + (v1 - v2)^2 / eps_v^2) <= 1.0
    
    properties
        EpsilonRange = 2.0;           % Spatial neighborhood radius in meters
        EpsilonVelocity = 2.5;        % Doppler velocity neighborhood tolerance in m/s
        MinPoints = 3;                % Minimum points required to form a core cluster
    end
    
    methods
        function obj = VelocityAwareDbscan(epsRange, epsVel, minPts)
            % Construct a VelocityAwareDbscan clustering instance.
            if nargin >= 1 && ~isempty(epsRange), obj.EpsilonRange = epsRange; end
            if nargin >= 2 && ~isempty(epsVel), obj.EpsilonVelocity = epsVel; end
            if nargin >= 3 && ~isempty(minPts), obj.MinPoints = minPts; end
        end
        
        function [clusterLabels, numClusters] = cluster(obj, pointCloud)
            % Execute Doppler-aware DBSCAN clustering on PointCloud4D instance.
            %
            % Returns:
            %   clusterLabels - 1 x N integer array (1..numClusters, 0 = noise outlier)
            %   numClusters   - Total number of identified valid clusters
            
            N = pointCloud.NumPoints;
            if N == 0
                clusterLabels = [];
                numClusters = 0;
                return;
            end
            
            pos = pointCloud.Positions;
            vel = pointCloud.Velocities;
            
            % Normalized 4D coordinates for standard unit-radius search:
            % [x/eps_r; y/eps_r; z/eps_r; v/eps_v]
            normPts = [pos(1, :) / obj.EpsilonRange; ...
                       pos(2, :) / obj.EpsilonRange; ...
                       pos(3, :) / obj.EpsilonRange; ...
                       vel / obj.EpsilonVelocity];
                       
            kd = chirp.pointcloud.SpatialKdTree(normPts);
            
            visited = false(1, N);
            clusterLabels = zeros(1, N);
            currentCluster = 0;
            
            for i = 1:N
                if visited(i)
                    continue;
                end
                visited(i) = true;
                
                % Query neighbors within normalized distance 1.0
                neighbors = kd.queryRadius(normPts(:, i), 1.0);
                
                if length(neighbors) < obj.MinPoints
                    clusterLabels(i) = 0; % Tentative noise
                else
                    currentCluster = currentCluster + 1;
                    clusterLabels(i) = currentCluster;
                    
                    % Expand cluster
                    seedQueue = neighbors(neighbors ~= i);
                    qIdx = 1;
                    
                    while qIdx <= length(seedQueue)
                        curPt = seedQueue(qIdx);
                        qIdx = qIdx + 1;
                        
                        if ~visited(curPt)
                            visited(curPt) = true;
                            curNeighbors = kd.queryRadius(normPts(:, curPt), 1.0);
                            if length(curNeighbors) >= obj.MinPoints
                                % Add new unique points to queue
                                newPoints = curNeighbors(~visited(curNeighbors));
                                seedQueue = unique([seedQueue, newPoints]);
                            end
                        end
                        
                        if clusterLabels(curPt) == 0
                            clusterLabels(curPt) = currentCluster;
                        end
                    end
                end
            end
            
            numClusters = currentCluster;
        end
    end
end
