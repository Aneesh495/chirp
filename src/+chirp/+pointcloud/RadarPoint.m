classdef RadarPoint
    % RADARPOINT 4D kinematic radar point measurement value object.
    %
    % Encapsulates 3D Cartesian coordinates, radial Doppler velocity,
    % SNR, radar cross section, spherical coordinates, and timestamp.
    
    properties
        X = 0.0;                      % Lateral cross-range coordinate (m)
        Y = 0.0;                      % Longitudinal range coordinate (m)
        Z = 0.0;                      % Vertical height coordinate (m)
        Velocity = 0.0;               % Radial line-of-sight velocity (m/s)
        SnrDb = 0.0;                  % Detection Signal-to-Noise Ratio (dB)
        Rcs = 1.0;                    % Estimated Radar Cross Section (m^2)
        AzimuthDeg = 0.0;             % Azimuth angle (deg)
        ElevationDeg = 0.0;           % Elevation angle (deg)
        Range = 0.0;                  % Slant range (m)
        Timestamp = 0.0;              % Measurement timestamp (s)
    end
    
    methods
        function obj = RadarPoint(varargin)
            % Construct RadarPoint with optional property overrides.
            if nargin > 0
                for k = 1:2:nargin
                    propName = varargin{k};
                    propVal = varargin{k+1};
                    if isprop(obj, propName)
                        obj.(propName) = propVal;
                    else
                        error('RadarPoint:UnknownProperty', 'Property %s not found.', propName);
                    end
                end
            end
        end
        
        function pos = getPosition(obj)
            % Returns 3x1 Cartesian position vector [x; y; z].
            pos = [obj.X; obj.Y; obj.Z];
        end
    end
end
