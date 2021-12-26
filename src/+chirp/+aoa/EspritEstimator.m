classdef EspritEstimator
    % ESPRITESTIMATOR Estimation of Signal Parameters via Rotational Invariance Techniques.
    %
    % Decomposes an array of identical doublet pairs into two displaced sub-arrays.
    % Solves for the rotational invariance operator Psi using Total Least Squares (TLS),
    % deriving arrival angles directly from the eigenvalues of Psi.
    
    properties
        DisplacementMeters            % Invariant baseline displacement between sub-arrays Delta (m)
        Wavelength                    % Carrier wavelength (m)
    end
    
    methods
        function obj = EspritEstimator(displacementMeters, wavelengthMeters)
            % Construct an ESPRIT estimator instance.
            if nargin >= 2
                obj.DisplacementMeters = displacementMeters;
                obj.Wavelength = wavelengthMeters;
            end
        end
        
        function anglesDeg = estimate(obj, spatialSnapshots, numSignals)
            % Estimate angles of arrival using Total Least Squares (TLS) ESPRIT.
            %
            % Inputs:
            %   spatialSnapshots - [M x N_snapshots] or [M x M] covariance matrix
            %   numSignals       - Number of source targets
            %
            % Output:
            %   anglesDeg        - 1 x numSignals estimated arrival angles in degrees
            
            if size(spatialSnapshots, 1) == size(spatialSnapshots, 2)
                R = spatialSnapshots;
            else
                nSnap = size(spatialSnapshots, 2);
                R = (spatialSnapshots * spatialSnapshots') / max(nSnap, 1);
            end
            
            % Forward-backward smoothing
            R = chirp.utils.MatrixUtils.forward_backward_avg(R);
            
            M = size(R, 1);
            [Es, ~, ~] = chirp.utils.MatrixUtils.subspace_decompose(R, numSignals);
            
            % Sub-array 1: elements 1 to M-1
            % Sub-array 2: elements 2 to M
            Es1 = Es(1:(M - 1), :);
            Es2 = Es(2:M, :);
            
            % Form concatenated matrix [Es1, Es2]
            E12 = [Es1, Es2];
            
            % SVD for Total Least Squares solution
            [~, ~, V] = svd(E12, 'econ');
            K = numSignals;
            
            % Partition V into K x K blocks
            V12 = V(1:K, (K + 1):(2 * K));
            V22 = V((K + 1):(2 * K), (K + 1):(2 * K));
            
            % Total Least Squares operator: Psi = -V12 * inv(V22)
            Psi = -V12 / V22;
            
            % Eigenvalues of Psi
            eigPsi = eig(Psi);
            
            % Phase shift: phi = arg(eigPsi) = 2*pi * Delta * sin(theta) / lambda
            phases = angle(eigPsi);
            sinTheta = (phases * obj.Wavelength) / (2.0 * pi * obj.DisplacementMeters);
            
            sinTheta = max(-1.0, min(1.0, sinTheta));
            anglesDeg = asin(sinTheta) * (180.0 / pi);
            anglesDeg = sort(anglesDeg(:).', 'ascend');
        end
    end
end
