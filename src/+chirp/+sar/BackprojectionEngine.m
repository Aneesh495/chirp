classdef BackprojectionEngine
    % BACKPROJECTIONENGINE Exact time-domain Backprojection Algorithm (BPA) for SAR.
    %
    % Reconstructs focused 2D/3D radar images from range-compressed pulse streams
    % for arbitrary flight paths, wide beamwidths, and irregular topographical terrain.
    
    properties
        Wavelength                    % Carrier wavelength (m)
        UseMex = true;                % Flag to call native C MEX kernel if compiled
    end
    
    methods
        function obj = BackprojectionEngine(wavelengthMeters, useMex)
            % Construct a BackprojectionEngine instance.
            if nargin >= 1, obj.Wavelength = wavelengthMeters; end
            if nargin >= 2, obj.UseMex = useMex; end
        end
        
        function sarImage = reconstruct2D(obj, rangeProfiles, rangeAxis, platformPositions, ...
                                          gridX, gridY, groundZ)
            % Execute 2D Backprojection Algorithm to form complex image I(x, y).
            %
            % Inputs:
            %   rangeProfiles     - [N_range x N_pulses] complex range-compressed data
            %   rangeAxis         - 1 x N_range vector of range coordinates (m)
            %   platformPositions - [3 x N_pulses] matrix of platform positions [x; y; z]
            %   gridX             - 1 x N_x vector of cross-range image coordinates (m)
            %   gridY             - 1 x N_y vector of ground-range image coordinates (m)
            %   groundZ           - (Optional) Ground elevation scalar or heightmap (default: 0)
            %
            % Output:
            %   sarImage          - [N_y x N_x] complex focused radar reflectivity image
            
            if nargin < 7 || isempty(groundZ)
                groundZ = 0.0;
            end
            
            [nRange, nPulses] = size(rangeProfiles);
            nX = length(gridX);
            nY = length(gridY);
            
            % Check if compiled MEX kernel is available and requested
            if obj.UseMex && exist('chirp_backprojection_mex', 'file') == 3
                sarImage = chirp_backprojection_mex(rangeProfiles, rangeAxis, platformPositions, ...
                                                   gridX, gridY, groundZ, obj.Wavelength);
                return;
            end
            
            % Pure vectorized MATLAB implementation
            [meshX, meshY] = meshgrid(gridX, gridY);
            sarImage = zeros(nY, nX);
            
            c = chirp.utils.PhysicalConstants.C;
            dr = rangeAxis(2) - rangeAxis(1);
            r0 = rangeAxis(1);
            lambda = obj.Wavelength;
            kWave = 4.0 * pi / lambda;
            
            % Loop across synthetic aperture pulses
            for u = 1:nPulses
                posU = platformPositions(:, u);
                xU = posU(1); yU = posU(2); zU = posU(3);
                
                % Compute Euclidean slant range from platform to all grid pixels
                dx = meshX - xU;
                dy = meshY - yU;
                dz = groundZ - zU;
                slantR = sqrt(dx.^2 + dy.^2 + dz.^2);
                
                % Convert slant range to continuous range bin indices
                binIdx = (slantR - r0) / dr + 1.0;
                
                % Bounds check: only interpolate pixels within valid range profile window
                validMask = (binIdx >= 1.0) & (binIdx <= (nRange - 1.0));
                
                profileU = rangeProfiles(:, u);
                
                % Linear interpolation across range profile
                idxFloor = floor(binIdx(validMask));
                frac = binIdx(validMask) - idxFloor;
                
                s1 = profileU(idxFloor);
                s2 = profileU(idxFloor + 1);
                interpSig = (1.0 - frac) .* s1 + frac .* s2;
                
                % Matched filter phase compensation: exp(j * 4*pi*R / lambda)
                phaseComp = exp(1j * kWave * slantR(validMask));
                
                % Accumulate into complex image
                currentImage = zeros(nY, nX);
                currentImage(validMask) = interpSig .* phaseComp;
                sarImage = sarImage + currentImage;
            end
        end
    end
end
