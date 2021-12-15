classdef MtiFilter
    % MTIFILTER Moving Target Indication (MTI) clutter rejection algorithms.
    %
    % Suppresses stationary ground and obstacle clutter in slow-time radar data
    % via delay-line cancellation, high-pass Butterworth IIR filtering, and
    % SVD principal clutter subspace projection.
    
    methods (Static)
        function filtered = twoPulseCanceler(data)
            % 2-pulse delay line canceler: y[m] = x[m] - x[m - 1] along dimension 2 (slow-time).
            % Frequency response: |H(f)| = 2 * |sin(pi * f * T)|
            
            dim = 2;
            nChirps = size(data, dim);
            if nChirps < 2
                filtered = data;
                return;
            end
            
            % Difference along dimension 2, zero-pad first column to preserve shape
            diffData = diff(data, 1, dim);
            firstSlice = zeros(size(data, 1), 1, size(data, 3));
            filtered = cat(dim, firstSlice, diffData);
        end
        
        function filtered = threePulseCanceler(data)
            % 3-pulse delay line canceler: y[m] = x[m] - 2*x[m - 1] + x[m - 2].
            % Frequency response: |H(f)| = 4 * sin^2(pi * f * T)
            
            dim = 2;
            nChirps = size(data, dim);
            if nChirps < 3
                filtered = data;
                return;
            end
            
            diffData = diff(data, 2, dim);
            padSlices = zeros(size(data, 1), 2, size(data, 3));
            filtered = cat(dim, padSlices, diffData);
        end
        
        function filtered = svdClutterFilter(data, numClutterSubspaces)
            % Subspace projection clutter filter using Singular Value Decomposition.
            %
            % Projects out the dominant zero-Doppler clutter modes spanning the
            % highest singular values across range-slowtime slices.
            
            if nargin < 2 || isempty(numClutterSubspaces)
                numClutterSubspaces = 1;
            end
            
            [nR, nC, nVirt] = size(data);
            filtered = zeros(nR, nC, nVirt);
            
            for iVirt = 1:nVirt
                X = squeeze(data(:, :, iVirt));
                [U, S, V] = svd(X, 'econ');
                
                % Zero out the first numClutterSubspaces singular values
                k = min(numClutterSubspaces, min(size(S)));
                S(1:k, 1:k) = 0.0;
                
                filtered(:, :, iVirt) = U * S * V';
            end
        end
    end
end
