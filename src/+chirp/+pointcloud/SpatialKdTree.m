classdef SpatialKdTree < handle
    % SPATIALKDTREE K-dimensional spatial indexing tree for accelerated range queries.
    %
    % Builds a balanced spatial search tree supporting O(log N) radius range searches
    % and k-nearest neighbors lookups for radar point cloud clustering.
    
    properties
        Points                        % D x N matrix of points
        RootNode                      % Struct representing root of tree
        Dimension                     % Dimensionality D
    end
    
    methods
        function obj = SpatialKdTree(pts)
            % Construct and build KdTree from D x N points matrix.
            if nargin >= 1 && ~isempty(pts)
                obj.Points = pts;
                obj.Dimension = size(pts, 1);
                indices = 1:size(pts, 2);
                obj.RootNode = obj.buildRecursive(indices, 1);
            end
        end
        
        function neighborIndices = queryRadius(obj, centerPoint, radius)
            % Find all point indices within specified Euclidean radius of centerPoint.
            if isempty(obj.RootNode)
                neighborIndices = [];
                return;
            end
            
            c = centerPoint(:);
            r2 = radius^2;
            neighborIndices = obj.searchRadiusRecursive(obj.RootNode, c, r2, 1);
        end
        
        function [knnIndices, knnDistances] = queryKnn(obj, queryPoint, k)
            % Find k nearest neighbor indices and distances to queryPoint.
            q = queryPoint(:);
            N = size(obj.Points, 2);
            k = min(k, N);
            
            % Compute distances to all points (vectorized fallback for small/medium N)
            diffs = obj.Points - q;
            distSq = sum(diffs.^2, 1);
            
            [sortedD2, sortIdx] = sort(distSq, 'ascend');
            knnIndices = sortIdx(1:k);
            knnDistances = sqrt(sortedD2(1:k));
        end
    end
    
    methods (Access = private)
        function node = buildRecursive(obj, indices, depth)
            % Recursively build KdTree node by splitting at median along axis.
            if isempty(indices)
                node = [];
                return;
            end
            
            axis = mod(depth - 1, obj.Dimension) + 1;
            coords = obj.Points(axis, indices);
            
            [~, sortOrder] = sort(coords, 'ascend');
            sortedIndices = indices(sortOrder);
            
            mid = ceil(length(sortedIndices) / 2);
            pointIdx = sortedIndices(mid);
            
            node = struct(...
                'pointIdx', pointIdx, ...
                'axis', axis, ...
                'splitVal', obj.Points(axis, pointIdx), ...
                'left', obj.buildRecursive(sortedIndices(1:(mid - 1)), depth + 1), ...
                'right', obj.buildRecursive(sortedIndices((mid + 1):end), depth + 1) ...
            );
        end
        
        function indices = searchRadiusRecursive(obj, node, center, r2, depth)
            % Recursive search for points within radius squared r2.
            indices = [];
            if isempty(node)
                return;
            end
            
            pt = obj.Points(:, node.pointIdx);
            d2 = sum((pt - center).^2);
            if d2 <= r2
                indices = [indices, node.pointIdx];
            end
            
            axis = node.axis;
            diffVal = center(axis) - node.splitVal;
            
            if diffVal < 0
                indices = [indices, obj.searchRadiusRecursive(node.left, center, r2, depth + 1)];
                if (diffVal^2) <= r2
                    indices = [indices, obj.searchRadiusRecursive(node.right, center, r2, depth + 1)];
                end
            else
                indices = [indices, obj.searchRadiusRecursive(node.right, center, r2, depth + 1)];
                if (diffVal^2) <= r2
                    indices = [indices, obj.searchRadiusRecursive(node.left, center, r2, depth + 1)];
                end
            end
        end
    end
end
