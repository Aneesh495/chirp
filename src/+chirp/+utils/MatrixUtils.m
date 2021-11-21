classdef MatrixUtils
    % MATRIXUTILS Numerical linear algebra routines for radar spatial covariance matrices.
    %
    % Implements regularized inversion with diagonal loading, forward-backward
    % spatial smoothing, exchange matrix permutations, and robust subspace decomposition.
    
    methods (Static)
        function J = exchange_matrix(N)
            % Generate N x N antidiagonal exchange matrix (reflection operator).
            % J(i, j) = 1 when i + j = N + 1, and 0 otherwise.
            J = fliplr(eye(N));
        end
        
        function Rfb = forward_backward_avg(R)
            % Forward-backward spatial covariance smoothing:
            %   R_fb = 0.5 * (R + J * conj(R) * J)
            %
            % Restores rank in coherent signal environments (e.g. multipath reflections)
            % across uniform linear arrays (ULA).
            
            [M, N] = size(R);
            if M ~= N
                error('MatrixUtils:SquareRequired', 'Covariance matrix must be square.');
            end
            
            J = chirp.utils.MatrixUtils.exchange_matrix(M);
            Rfb = 0.5 * (R + J * conj(R) * J);
            % Ensure Hermitian symmetry
            Rfb = 0.5 * (Rfb + Rfb');
        end
        
        function Rinv = regularized_inv(R, loadingFactor)
            % Regularized matrix inversion with adaptive diagonal loading:
            %   R_reg = R + (loadingFactor * trace(R) / M) * I
            %
            % Used in Capon (MVDR) beamforming when covariance is ill-conditioned
            % or sample count is limited.
            
            if nargin < 2 || isempty(loadingFactor)
                loadingFactor = 1e-3;
            end
            
            [M, N] = size(R);
            if M ~= N
                error('MatrixUtils:SquareRequired', 'Covariance matrix must be square.');
            end
            
            tr = real(trace(R));
            if tr <= 0
                tr = 1.0;
            end
            
            delta = (loadingFactor * tr / M);
            Rreg = R + delta * eye(M);
            Rinv = inv(Rreg);
        end
        
        function [Es, En, lambda] = subspace_decompose(R, numSignals)
            % Eigendecomposition and subspace separation for MUSIC and ESPRIT algorithms.
            %
            % Inputs:
            %   R          - M x M spatial covariance matrix
            %   numSignals - Number of dominant signal sources (rank of signal subspace)
            %
            % Outputs:
            %   Es         - M x K signal subspace eigenvector matrix
            %   En         - M x (M-K) noise subspace eigenvector matrix
            %   lambda     - Sorted eigenvalues in descending order
            
            [M, N] = size(R);
            if M ~= N
                error('MatrixUtils:SquareRequired', 'Covariance matrix must be square.');
            end
            
            % Force exact Hermitian symmetry
            R = 0.5 * (R + R');
            
            [V, D] = eig(R);
            eigenvals = real(diag(D));
            
            % Sort eigenvalues in descending order
            [lambda, idx] = sort(eigenvals, 'descend');
            V = V(:, idx);
            
            % Bound number of signals
            numSignals = max(1, min(M - 1, round(numSignals)));
            
            Es = V(:, 1:numSignals);
            En = V(:, (numSignals + 1):M);
        end
        
        function fNorm = frobenius_norm(A)
            % Compute Frobenius norm of matrix A: sqrt(sum(abs(A(:)).^2))
            fNorm = sqrt(sum(abs(A(:)).^2));
        end
    end
end
