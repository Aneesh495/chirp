classdef MusicEstimator
    % MUSICESTIMATOR MUltiple SIgnal Classification super-resolution spatial spectrum.
    %
    % Decomposes spatial covariance into signal and noise subspaces. Orthogonality
    % of steering vectors to noise subspace yields sharp spatial pseudo-spectrum peaks:
    %   P_music(theta) = 1 / (a(theta)' * En * En' * a(theta))
    
    properties
        MimoArray                     % chirp.sim.MimoArray instance
        Wavelength                    % Carrier wavelength (m)
        ForwardBackward = true;       % Apply forward-backward covariance averaging
    end
    
    methods
        function obj = MusicEstimator(arrayObj, wavelengthMeters, forwardBackward)
            % Construct a MusicEstimator instance.
            if nargin >= 2
                obj.MimoArray = arrayObj;
                obj.Wavelength = wavelengthMeters;
            end
            if nargin >= 3 && ~isempty(forwardBackward)
                obj.ForwardBackward = forwardBackward;
            end
        end
        
        function [spectrumDb, azGrid] = estimateAzimuth(obj, spatialSnapshots, numSignals, azGridDeg, elDeg)
            % Estimate 1D MUSIC pseudo-spectrum across azimuth grid in degrees.
            %
            % Inputs:
            %   spatialSnapshots - [N_virt x N_snapshots] or [N_virt x N_virt] covariance
            %   numSignals       - Estimated number of discrete target arrivals
            %   azGridDeg        - Scanning azimuth grid (default: -60:0.25:60)
            %   elDeg            - Fixed elevation angle (default: 0)
            
            if nargin < 4 || isempty(azGridDeg)
                azGridDeg = linspace(-60.0, 60.0, 481);
            end
            if nargin < 5 || isempty(elDeg)
                elDeg = 0.0;
            end
            
            if size(spatialSnapshots, 1) == size(spatialSnapshots, 2)
                R = spatialSnapshots;
            else
                nSnap = size(spatialSnapshots, 2);
                R = (spatialSnapshots * spatialSnapshots') / max(nSnap, 1);
            end
            
            % Forward-backward spatial smoothing if requested
            if obj.ForwardBackward
                R = chirp.utils.MatrixUtils.forward_backward_avg(R);
            end
            
            % Eigendecomposition and noise subspace isolation
            [~, En, ~] = chirp.utils.MatrixUtils.subspace_decompose(R, numSignals);
            
            % Noise projection matrix: P_n = En * En'
            Pn = En * En';
            
            nGrid = length(azGridDeg);
            spectrum = zeros(1, nGrid);
            
            for k = 1:nGrid
                az = azGridDeg(k);
                sv = obj.MimoArray.steeringVector(az, elDeg, obj.Wavelength);
                denom = real(sv' * Pn * sv);
                spectrum(k) = 1.0 / max(denom, 1e-12);
            end
            
            spectrumDb = 10.0 * log10(max(spectrum, eps));
            azGrid = azGridDeg;
        end
    end
end
