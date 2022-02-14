classdef SarGeometry
    % SARGEOMETRY Flight trajectory kinematics and ground coordinate mapping for SAR.
    %
    % Models platform motion for Stripmap, Spotlight, and Circular Synthetic Aperture Radar (SAR),
    % computing instant slant ranges to arbitrary ground grid pixels.
    
    properties
        FlightMode = 'stripmap';      % SAR mode: 'stripmap', 'spotlight', or 'circular'
        FlightAltitude = 1000.0;      % Platform altitude h above ground (m)
        PlatformVelocity = 100.0;     % Flight speed v along track (m/s)
        SceneCenter = [0.0; 5000.0; 0.0]; % 3x1 target scene center coordinate [x; y; z] (m)
        ApertureLength = 200.0;       % Total synthetic aperture baseline length L (m)
    end
    
    methods
        function obj = SarGeometry(varargin)
            % Construct a SarGeometry instance.
            if nargin > 0
                for k = 1:2:nargin
                    propName = varargin{k};
                    propVal = varargin{k+1};
                    if isprop(obj, propName)
                        obj.(propName) = propVal;
                    else
                        error('SarGeometry:UnknownProperty', 'Property %s not found.', propName);
                    end
                end
            end
        end
        
        function pos = platformPosition(obj, uApertureCoord)
            % Returns 3x1 platform position [x; y; z] at synthetic aperture coordinate u.
            %
            % uApertureCoord is in [-L/2, L/2] along cross-range X axis.
            
            switch lower(obj.FlightMode)
                case 'stripmap'
                    % Straight-and-level flight along X axis at constant altitude z = h
                    pos = [uApertureCoord; 0.0; obj.FlightAltitude];
                    
                case 'spotlight'
                    % Linear flight with boresight continuously steered toward SceneCenter
                    pos = [uApertureCoord; 0.0; obj.FlightAltitude];
                    
                case 'circular'
                    % Circular orbit of radius R around SceneCenter in horizontal plane
                    radius = norm(obj.SceneCenter(1:2));
                    ang = (uApertureCoord / obj.ApertureLength) * (2.0 * pi);
                    pos = [radius * sin(ang); radius * cos(ang); obj.FlightAltitude];
                    
                otherwise
                    pos = [uApertureCoord; 0.0; obj.FlightAltitude];
            end
        end
        
        function r = slantRange(obj, platformPos, pixelCoord)
            % Compute Euclidean distance between platform and a ground pixel.
            diffVec = platformPos(:) - pixelCoord(:);
            r = sqrt(sum(diffVec.^2));
        end
    end
end
