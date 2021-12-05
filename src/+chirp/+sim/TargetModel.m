classdef TargetModel
    % TARGETMODEL Kinematic target state and electromagnetic radar cross section.
    %
    % Represents a dynamic 3D physical scatterer characterized by position,
    % velocity, acceleration, mean radar cross section, and Swerling index.
    
    properties
        Position = [0.0; 20.0; 0.0];      % 3x1 Cartesian position [x; y; z] in meters
        Velocity = [0.0; 0.0; 0.0];       % 3x1 Cartesian velocity [vx; vy; vz] in m/s
        Acceleration = [0.0; 0.0; 0.0];   % 3x1 Cartesian acceleration [ax; ay; az] in m/s^2
        Rcs = 1.0;                        % Mean radar cross section in square meters (m^2)
        SwerlingType = 0;                 % Swerling model index (0, 1, 2, 3, or 4)
        TargetId = 1;                     % Unique integer tracking identifier
    end
    
    properties (Dependent)
        Range                 % Slant range from sensor origin (m)
        AzimuthDeg            % Azimuth angle relative to sensor boresight (+y) in degrees
        ElevationDeg          % Elevation angle relative to horizontal plane in degrees
        RadialVelocity        % Line-of-sight velocity (m/s, positive = receding)
    end
    
    methods
        function obj = TargetModel(varargin)
            % Construct a TargetModel instance with optional name-value parameters.
            if nargin > 0
                for k = 1:2:nargin
                    propName = varargin{k};
                    propVal = varargin{k+1};
                    if isprop(obj, propName)
                        obj.(propName) = propVal;
                    else
                        error('TargetModel:UnknownProperty', 'Property %s not found in TargetModel.', propName);
                    end
                end
            end
            obj.Position = obj.Position(:);
            obj.Velocity = obj.Velocity(:);
            obj.Acceleration = obj.Acceleration(:);
        end
        
        function r = get.Range(obj)
            r = sqrt(sum(obj.Position.^2));
        end
        
        function az = get.AzimuthDeg(obj)
            % Boresight along +y axis; x is lateral cross-range
            az = atan2(obj.Position(1), obj.Position(2)) * (180.0 / pi);
        end
        
        function el = get.ElevationDeg(obj)
            horizDist = sqrt(obj.Position(1)^2 + obj.Position(2)^2);
            el = atan2(obj.Position(3), horizDist) * (180.0 / pi);
        end
        
        function vr = get.RadialVelocity(obj)
            % Radial velocity vr = (r . v) / ||r||
            r = obj.Range;
            if r <= 1e-6
                vr = 0.0;
            else
                vr = sum(obj.Position .* obj.Velocity) / r;
            end
        end
        
        function updatedTarget = propagate(obj, dtSeconds)
            % Propagate kinematic state forward by dtSeconds using Newtonian mechanics.
            %   p(t + dt) = p(t) + v(t)*dt + 0.5*a(t)*dt^2
            %   v(t + dt) = v(t) + a(t)*dt
            
            updatedTarget = obj;
            p = obj.Position + obj.Velocity * dtSeconds + 0.5 * obj.Acceleration * (dtSeconds^2);
            v = obj.Velocity + obj.Acceleration * dtSeconds;
            
            updatedTarget.Position = p;
            updatedTarget.Velocity = v;
        end
    end
end
