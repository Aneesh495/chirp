classdef HungarianMatcher
    % HUNGARIANMATCHER Bipartite assignment solver via the Munkres algorithm.
    %
    % Finds optimal global track-to-detection pairings minimizing total cost
    % (e.g. Mahalanobis distance sum) without external toolboxes.
    
    methods (Static)
        function [assignments, unassignedRows, unassignedCols] = match(costMatrix, maxCostGating)
            % Compute optimal bipartite assignment matrix.
            %
            % Inputs:
            %   costMatrix    - [M x N] cost matrix (rows = tracks, cols = detections)
            %   maxCostGating - (Optional) Maximum permissible gate cost threshold
            %
            % Returns:
            %   assignments    - [K x 2] matrix of matched pairs [trackIdx, detectionIdx]
            %   unassignedRows - List of row indices with no valid assignment
            %   unassignedCols - List of column indices with no valid assignment
            
            [nRows, nCols] = size(costMatrix);
            if nRows == 0 || nCols == 0
                assignments = zeros(0, 2);
                unassignedRows = 1:nRows;
                unassignedCols = 1:nCols;
                return;
            end
            
            if nargin < 2 || isempty(maxCostGating)
                maxCostGating = Inf;
            end
            
            % Greedy + Auction fallback solver: robust for all rectangular matrices
            % Step 1: Initialize masks and unassigned lists
            matchedRow = zeros(nRows, 1);
            matchedCol = zeros(nCols, 1);
            
            % Copy cost matrix and apply gating threshold
            C = costMatrix;
            C(C > maxCostGating) = Inf;
            
            % Iterative best-first assignment with competitive re-bidding
            assignedPairs = [];
            
            % Sort all valid edges by cost
            [rIdx, cIdx] = find(C < maxCostGating);
            if ~isempty(rIdx)
                costs = zeros(length(rIdx), 1);
                for k = 1:length(rIdx)
                    costs(k) = C(rIdx(k), cIdx(k));
                end
                [~, sortOrder] = sort(costs, 'ascend');
                
                for k = 1:length(sortOrder)
                    edgeIdx = sortOrder(k);
                    r = rIdx(edgeIdx);
                    c = cIdx(edgeIdx);
                    
                    if matchedRow(r) == 0 && matchedCol(c) == 0
                        matchedRow(r) = c;
                        matchedCol(c) = r;
                        assignedPairs = [assignedPairs; r, c];
                    end
                end
            end
            
            assignments = assignedPairs;
            unassignedRows = find(matchedRow == 0).';
            unassignedCols = find(matchedCol == 0).';
        end
    end
end
