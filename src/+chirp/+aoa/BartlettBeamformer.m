classdef BartlettBeamformer
    % BARTLETTBEAMFORMER Conventional Fourier spatial beamformer across virtual antenna arrays.
    %
    % Computes spatial power spectrum by steering phase weights across the array:
    %   P(theta, phi) = (a(theta, phi)' * R * a(theta, phi)) / (a' * a)
    
    properties
        MimoArray                     % chirp.sim.MimoArray instance
        Wavelength                    % Carrier wavelength (m)
    end
    
    methods
        function obj = BartlettBeamformer(arrayObj, wavelengthMeters)
            % Construct a BartlettBeamformer instance.
            if nargin >= 2
                obj.MimoArray = arrayObj;
                obj.Wavelength = wavelengthMeters;
            end
        end
        
        function [spectrumDb, azGrid] = estimateAzimuth(obj, spatialSnapshots, azGridDeg, elDeg)
            % Estimate 1D azimuth spatial spectrum across specified grid in degrees.
            %
            % Inputs:
            %   spatialSnapshots - [N_virt x N_snapshots] complex matrix or [N_virt x N_virt] covariance
            %   azGridDeg        - Vector of azimuth angles to scan (default: -60:0.5:60)
            %   elDeg            - Elevation angle in degrees (default: 0)
            %
            % Returns:
            %   spectrumDb       - 1 x N_grid spatial power in decibels
            %   azGrid           - 1 x N_grid evaluation grid
            
            if nargin < 3 || isempty(azGridDeg)
                azGridDeg = linspace(-60.0, 60.0, 241);
            end
            if nargin < 4 || isempty(elDeg)
                elDeg = 0.0;
            end
            
            % Compute or retrieve spatial covariance matrix R
            if size(spatialSnapshots, 1) == size(spatialSnapshots, 2)
                R = spatialSnapshots;
            else
                nSnap = size(spatialSnapshots, 2);
                R = (spatialSnapshots * spatialSnapshots') / max(nSnap, 1);
            end
            
            nGrid = length(azGridDeg);
            spectrum = zeros(1, nGrid);
            
            for k = 1:nGrid
                az = azGridDeg(k);
                sv = obj.MimoArray.steeringVector(az, elDeg, obj.Wavelength);
                normFactor = real(sv' * sv);
                spectrum(k) = real(sv' * R * sv) / max(normFactor, eps);
            end
            
            spectrumDb = 10.0 * log10(max(spectrum, eps));
            azGrid = azGridDeg;
        end
    end
end
