classdef CaponBeamformer
    % CAPONBEAMFORMER Minimum Variance Distortionless Response (MVDR) spatial estimator.
    %
    % Minimizes total output power subject to unity distortionless gain in the look direction:
    %   P_capon(theta) = 1 / (a(theta)' * R^(-1) * a(theta))
    
    properties
        MimoArray                     % chirp.sim.MimoArray instance
        Wavelength                    % Carrier wavelength (m)
        LoadingFactor = 1e-3;         % Diagonal loading regularization factor
    end
    
    methods
        function obj = CaponBeamformer(arrayObj, wavelengthMeters, loadingFactor)
            % Construct a CaponBeamformer instance.
            if nargin >= 2
                obj.MimoArray = arrayObj;
                obj.Wavelength = wavelengthMeters;
            end
            if nargin >= 3 && ~isempty(loadingFactor)
                obj.LoadingFactor = loadingFactor;
            end
        end
        
        function [spectrumDb, azGrid] = estimateAzimuth(obj, spatialSnapshots, azGridDeg, elDeg)
            % Estimate 1D Capon MVDR spatial pseudo-spectrum in decibels.
            
            if nargin < 3 || isempty(azGridDeg)
                azGridDeg = linspace(-60.0, 60.0, 241);
            end
            if nargin < 4 || isempty(elDeg)
                elDeg = 0.0;
            end
            
            if size(spatialSnapshots, 1) == size(spatialSnapshots, 2)
                R = spatialSnapshots;
            else
                nSnap = size(spatialSnapshots, 2);
                R = (spatialSnapshots * spatialSnapshots') / max(nSnap, 1);
            end
            
            % Regularized matrix inversion with diagonal loading
            Rinv = chirp.utils.MatrixUtils.regularized_inv(R, obj.LoadingFactor);
            
            nGrid = length(azGridDeg);
            spectrum = zeros(1, nGrid);
            
            for k = 1:nGrid
                az = azGridDeg(k);
                sv = obj.MimoArray.steeringVector(az, elDeg, obj.Wavelength);
                denom = real(sv' * Rinv * sv);
                spectrum(k) = 1.0 / max(denom, 1e-12);
            end
            
            spectrumDb = 10.0 * log10(max(spectrum, eps));
            azGrid = azGridDeg;
        end
    end
end
