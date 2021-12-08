classdef ClutterGenerator
    % CLUTTERGENERATOR Radar surface and volume clutter returns.
    %
    % Simulates distributed ground, road surface, and sea clutter reflections
    % using Weibull, K-distribution, and constant sigma-0 area backscattering models.
    
    properties
        ClutterType = 'weibull';       % Clutter distribution: 'constant', 'weibull', or 'k-distribution'
        Sigma0Dbm2 = -20.0;            % Area reflectivity sigma_0 in dB(m^2/m^2)
        WeibullShape = 1.5;            % Weibull shape parameter (b)
        KDistShape = 2.0;              % K-distribution shape parameter (nu)
        MinRange = 2.0;                % Minimum clutter generation range (m)
        MaxRange = 80.0;               % Maximum clutter generation range (m)
        NumPatches = 150;              % Number of discrete discrete clutter scatterers
    end
    
    methods
        function obj = ClutterGenerator(varargin)
            % Construct a ClutterGenerator with optional name-value overrides.
            if nargin > 0
                for k = 1:2:nargin
                    propName = varargin{k};
                    propVal = varargin{k+1};
                    if isprop(obj, propName)
                        obj.(propName) = propVal;
                    else
                        error('ClutterGenerator:UnknownProperty', 'Property %s not found.', propName);
                    end
                end
            end
        end
        
        function clutterTargets = generatePatches(obj, sensorHeight)
            % Generate synthetic discrete clutter scatterer patches on ground plane (z = 0).
            %
            % Inputs:
            %   sensorHeight - Radar height above ground in meters (default: 1.0 m)
            %
            % Returns:
            %   clutterTargets - Array of chirp.sim.TargetModel objects
            
            if nargin < 2 || isempty(sensorHeight)
                sensorHeight = 1.0;
            end
            
            N = obj.NumPatches;
            clutterTargets = chirp.sim.TargetModel.empty(0, N);
            
            % Random distribution of ranges and azimuth angles on ground
            r = obj.MinRange + (obj.MaxRange - obj.MinRange) * sqrt(rand(1, N));
            azDeg = -60.0 + 120.0 * rand(1, N);
            
            % Area of each annular ring patch
            dr = (obj.MaxRange - obj.MinRange) / sqrt(N);
            dAzRad = (120.0 * pi / 180.0) / N;
            patchArea = r .* dr .* dAzRad;
            
            meanRcsBase = chirp.utils.PhysicalConstants.db2pow(obj.Sigma0Dbm2) .* patchArea;
            
            switch lower(obj.ClutterType)
                case 'constant'
                    rcsVals = meanRcsBase;
                case 'weibull'
                    % Weibull random variable generation
                    u = max(rand(1, N), 1e-12);
                    rcsVals = meanRcsBase .* ((-log(u)).^(1.0 / obj.WeibullShape));
                case 'k-distribution'
                    % K-distribution: product of gamma (underlying texture) and exponential (speckle)
                    nu = obj.KDistShape;
                    % Gamma random variable with shape nu and scale 1/nu
                    % Approximation via sum of exponentials when nu is small
                    tex = (-log(max(rand(1, N), 1e-12)) + -log(max(rand(1, N), 1e-12))) / 2.0;
                    speckle = -log(max(rand(1, N), 1e-12));
                    rcsVals = meanRcsBase .* tex .* speckle;
                otherwise
                    rcsVals = meanRcsBase;
            end
            
            for k = 1:N
                % Ground plane coordinates (z = -sensorHeight relative to radar)
                azRad = azDeg(k) * (pi / 180.0);
                x = r(k) * sin(azRad);
                y = r(k) * cos(azRad);
                z = -sensorHeight;
                
                clutterTargets(k) = chirp.sim.TargetModel(...
                    'Position', [x; y; z], ...
                    'Velocity', [0.0; 0.0; 0.0], ...
                    'Rcs', max(rcsVals(k), 1e-6), ...
                    'SwerlingType', 0, ...
                    'TargetId', 1000 + k ...
                );
            end
        end
    end
end
